## Builds the binary payload behind the "Data browser" article.
##
## Reads the per-dataset RDS files from DATA_DIR plus the shared PoPS table and
## writes a static bundle to pkgdown/assets/databrowser/, which pkgdown copies to
## the site root. Re-run whenever a dataset changes, then rebuild the site.
##
##   Rscript data-raw/build_databrowser_data.R            # everything
##   Rscript data-raw/build_databrowser_data.R dlpfc liver  # only matching ids
##
## Visium datasets are stored per spot with exact raw counts. CosMx datasets hold
## 10^5-10^6 cells, which is far past both the repo budget and what a ~350px panel
## can resolve, so their cells are aggregated onto a spatial grid of about
## TARGET_BINS bins and each bin carries the MEAN raw count of its cells.
##
## Layout written:
##   databrowser/manifest.json          dataset list + disease list
##   databrowser/pops.bin               float32[nPopsGenes*nDiseases], gene-major
##   databrowser/pops_genes.json        gene order for pops.bin
##   databrowser/<id>/meta.json         genes, cell types, per-gene max/scale, pops map
##   databrowser/<id>/coords.bin        float32[nSpots*2], x,y interleaved
##   databrowser/<id>/celltype.bin      float32[nSpots*nCellTypes], spot-major
##   databrowser/<id>/esv.bin           float32[nCellTypes*nGenes], cell-type-major
##   databrowser/<id>/expr_NNN.bin      uint16[chunk*nSpots], gene-major
##                                      value = stored * exprScale[gene]

suppressMessages(library(Matrix))

DATA_DIR <- "/Users/leeh16/Documents/spacelink/GitHub_repo_data_browser"
OUT <- "pkgdown/assets/databrowser"
TARGET_BINS <- 20000L
BYTES_PER_CHUNK <- 700000    # aim for a ~0.7 MB per-gene fetch

datasets <- list(
  list(id = "visium_human_dlpfc",         label = "Visium human DLPFC",
       prefix = "Visium_human_DLPFC",         bin = FALSE),
  list(id = "visium_human_liver",         label = "Visium human liver",
       prefix = "Visium_human_liver",         bin = FALSE),
  list(id = "visium_human_lymph_node",    label = "Visium human lymph node",
       prefix = "Visium_human_lymph_node",    bin = FALSE),
  list(id = "cosmx_human_frontal_cortex", label = "CosMx human frontal cortex",
       prefix = "CosMx_human_frontal_cortex", bin = TRUE),
  list(id = "cosmx_human_liver",          label = "CosMx human liver",
       prefix = "CosMx_human_liver",          bin = TRUE),
  list(id = "cosmx_human_lymph",          label = "CosMx human lymph node",
       prefix = "CosMx_human_lymph",          bin = TRUE)
)

sel <- commandArgs(trailingOnly = TRUE)
if (length(sel)) {
  datasets <- Filter(function(d) any(grepl(paste(sel, collapse = "|"), d$id)), datasets)
  message("Selected: ", paste(vapply(datasets, `[[`, "", "id"), collapse = ", "))
}

# ---- io helpers -------------------------------------------------------------

write_f32 <- function(x, path) {
  con <- file(path, "wb"); on.exit(close(con))
  writeBin(as.double(x), con, size = 4L, endian = "little")
}

write_u16 <- function(x, path) {
  stopifnot(all(is.finite(x)), min(x) >= 0, max(x) <= 65535)
  con <- file(path, "wb"); on.exit(close(con))
  writeBin(as.integer(x), con, size = 2L, endian = "little")
}

json <- function(x, path) {
  writeLines(jsonlite::toJSON(x, auto_unbox = TRUE, digits = 8, na = "null"), path)
}

# ---- ESV gene-name reconciliation -------------------------------------------
#
# The ESV tables were round-tripped through a spreadsheet, which mangled some
# symbols: SEPT2 -> "2-Sep", LSM3 -> "LSM<nbsp>3.00", COP1 -> "COP 1", and -
# via scientific notation - CYP2E1 -> "CYP 20.00". Reversing those blind is
# ambiguous, so instead each real symbol from the counts matrix is mangled
# forward into the keys it could have become and matched against the leftovers.

esv_key <- function(x) {
  x <- gsub(" ", " ", x)
  x <- sub("\\.0+$", "", x)
  toupper(gsub("[^A-Za-z0-9]", "", x))
}

count_keys <- function(g) {
  up <- toupper(g)
  keys <- gsub("[^A-Z0-9]", "", up)
  m <- regmatches(up, regexec(
    "^(SEPT|SEP|MARCH|MARC|MAR|DEC|APR|FEB|OCT|NOV|JUN|JUL|AUG|JAN|MAY)([0-9]+)$", up))[[1]]
  if (length(m) == 3L) keys <- c(keys, paste0(m[3], substr(m[2], 1, 3)))
  e <- regmatches(up, regexec("^([A-Z]+)([0-9]+)E([0-9]+)$", up))[[1]]
  if (length(e) == 4L) {
    v <- as.numeric(e[3]) * 10^as.numeric(e[4])
    if (is.finite(v)) keys <- c(keys, paste0(e[2], format(v, scientific = FALSE, trim = TRUE)))
  }
  unique(keys)
}

## Returns esv_names with mangled entries renamed to their real symbols.
repair_names <- function(esv_names, counts_genes) {
  lost <- setdiff(counts_genes, esv_names)
  if (!length(lost)) return(esv_names)
  spare <- setdiff(esv_names, counts_genes)
  sk <- esv_key(spare)
  fixed <- character(0)
  for (g in lost) {
    hit <- spare[match(count_keys(g), sk)]
    hit <- hit[!is.na(hit)]
    if (length(hit)) fixed[hit[1]] <- g
  }
  if (anyDuplicated(unname(fixed))) stop("gene-name repair produced duplicate targets")
  hit <- match(esv_names, names(fixed))
  esv_names[!is.na(hit)] <- unname(fixed)[hit[!is.na(hit)]]
  esv_names
}

## The cell-type matrices and the ESV tables occasionally disagree on spelling
## (the liver sets label a column "Erthyroid.cells" where ESV says
## "Erythroid.cells"). Align on the ESV spelling, but only for a near-identical
## unambiguous match, so a genuinely different cell type is never silently
## folded into another.
repair_cell_types <- function(ct_cols, esv_types) {
  lost <- setdiff(ct_cols, esv_types)
  if (!length(lost)) return(ct_cols)
  spare <- setdiff(esv_types, c("whole", ct_cols))
  for (nm in lost) {
    d <- as.vector(adist(nm, spare, ignore.case = TRUE))
    near <- spare[d <= 2 & d == min(d)]
    if (length(near) != 1L) {
      stop("cell type '", nm, "' has no unambiguous ESV counterpart; candidates: ",
           if (length(near)) paste(near, collapse = ", ") else "none")
    }
    message("  cell type '", nm, "' -> '", near, "' (ESV spelling)")
    ct_cols[ct_cols == nm] <- near
    spare <- setdiff(spare, near)
  }
  ct_cols
}

# ---- spatial binning --------------------------------------------------------

## Picks the finest square-ish grid whose occupied-bin count reaches `target`.
choose_grid <- function(xy, target) {
  rx <- range(xy[, 1]); ry <- range(xy[, 2])
  dx <- diff(rx); dy <- diff(ry)
  ar <- if (dy > 0) dx / dy else 1
  assign_bins <- function(nx) {
    ny <- max(1L, as.integer(round(nx / ar)))
    ix <- pmin(nx, 1L + as.integer(floor((xy[, 1] - rx[1]) / dx * nx)))
    iy <- pmin(ny, 1L + as.integer(floor((xy[, 2] - ry[1]) / dy * ny)))
    list(nx = nx, ny = ny, key = (iy - 1L) * nx + ix)
  }
  lo <- 10L; hi <- 2000L
  while (lo < hi) {
    mid <- (lo + hi) %/% 2L
    if (length(unique(assign_bins(mid)$key)) < target) lo <- mid + 1L else hi <- mid
  }
  g <- assign_bins(lo)
  u <- sort(unique(g$key))
  bx <- ((u - 1L) %% g$nx) + 1L
  by <- ((u - 1L) %/% g$nx) + 1L
  list(bin = match(g$key, u), nBins = length(u),
       cx = rx[1] + (bx - 0.5) * dx / g$nx,
       cy = ry[1] + (by - 0.5) * dy / g$ny,
       nx = g$nx, ny = g$ny)
}

# ---- per-dataset build ------------------------------------------------------

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(OUT, "pops"), recursive = TRUE)   # superseded by pops.bin
info <- list()

for (ds in datasets) {
  message("== ", ds$label, " ==")
  f <- function(suffix) file.path(DATA_DIR, paste0(ds$prefix, suffix))

  counts <- readRDS(f("_counts.rds"))
  # Guarded so an already-CSC matrix is never copied - the largest dataset holds
  # ~5e8 non-zeros and a spare copy alone exhausts the vector limit.
  if (!is(counts, "CsparseMatrix")) counts <- as(as.matrix(counts), "CsparseMatrix")
  coords <- as.matrix(readRDS(f("_spatial_coords.rds")))[, 1:2, drop = FALSE]
  ct <- as.matrix(readRDS(f("_cell_type_data.rds")))
  esv <- readRDS(f("_ESV_results.rds"))

  # Some datasets carry a few spots in the counts matrix that the cell-type
  # table lacks; keep the intersection so every panel plots the same positions.
  cells <- Reduce(intersect, list(colnames(counts), rownames(coords), rownames(ct)))
  cells <- colnames(counts)[colnames(counts) %in% cells]   # keep counts order
  dropped <- ncol(counts) - length(cells)
  if (dropped > 0) {
    message("  dropping ", dropped, " of ", ncol(counts),
            " spot(s) absent from the coordinate or cell-type table")
  }
  stopifnot(length(cells) > 0)
  if (dropped > 0) counts <- counts[, cells, drop = FALSE]
  coords <- coords[cells, , drop = FALSE]
  ct <- ct[cells, , drop = FALSE]
  stopifnot(all(counts@x == round(counts@x)), min(counts@x) >= 0)

  # Row order is applied later, to the (much smaller) binned matrix, so the
  # full-size counts matrix is never reordered in place.
  genes <- sort(rownames(counts))
  gene_row <- match(genes, rownames(counts))
  colnames(ct) <- repair_cell_types(colnames(ct), unique(esv$Celltype))
  cell_types <- c("whole", colnames(ct))
  stopifnot(setequal(unique(esv$Celltype), cell_types))

  # ESV matrix, cell-type-major, aligned to `genes`
  esv_mat <- matrix(NA_real_, nrow = length(cell_types), ncol = length(genes),
                    dimnames = list(cell_types, genes))
  # A handful of (gene, cell type) pairs legitimately have no ESV - the score
  # was not computed there - and those stay NA. A large shortfall instead means
  # the names failed to reconcile, which is a bug worth stopping for.
  for (cty in cell_types) {
    s <- esv[esv$Celltype == cty, ]
    s$GeneName <- repair_names(s$GeneName, genes)
    hit <- match(genes, s$GeneName)
    n_missing <- sum(is.na(hit))
    if (n_missing > max(50L, 0.01 * length(genes))) {
      stop(ds$id, ": ESV for '", cty, "' is missing ", n_missing,
           " of ", length(genes), " genes - name reconciliation likely failed; e.g. ",
           paste(head(genes[is.na(hit)], 5), collapse = ", "))
    }
    if (n_missing > 0) {
      message("  ESV for '", cty, "': ", n_missing, " gene(s) without a score")
    }
    esv_mat[cty, ] <- s$ESV[hit]
  }
  rm(esv); gc(FALSE)

  # Aggregate CosMx cells onto a grid; Visium spots are kept as they are.
  if (ds$bin) {
    g <- choose_grid(coords, TARGET_BINS)
    message("  binning ", length(cells), " cells -> ", g$nBins,
            " bins on a ", g$nx, "x", g$ny, " grid")
    n_pos <- g$nBins
    per_bin <- tabulate(g$bin, n_pos)
    # averaging operator: cells x bins, each column summing to 1
    B <- sparseMatrix(i = seq_along(g$bin), j = g$bin, x = 1 / per_bin[g$bin],
                      dims = c(length(g$bin), n_pos))
    out_coords <- cbind(g$cx, g$cy)
    out_ct <- as.matrix(Matrix::crossprod(B, ct))
    # One sparse product for the whole matrix; chunked row-slicing of a CSC
    # matrix this large would rescan every non-zero on each chunk.
    message("  aggregating expression ...")
    expr <- as.matrix(counts %*% B)[gene_row, , drop = FALSE]
    rm(B); rm(counts); gc(FALSE)
  } else {
    n_pos <- length(cells)
    out_coords <- coords
    out_ct <- ct
    expr <- counts[gene_row, , drop = FALSE]
    rm(counts); gc(FALSE)
  }
  rm(ct, coords); gc(FALSE)

  dir <- file.path(OUT, ds$id)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  write_f32(as.vector(t(out_coords)), file.path(dir, "coords.bin"))
  write_f32(as.vector(t(out_ct)), file.path(dir, "celltype.bin"))
  write_f32(as.vector(t(esv_mat)), file.path(dir, "esv.bin"))

  n_genes <- length(genes)
  chunk <- max(1L, min(64L, as.integer(floor(BYTES_PER_CHUNK / (n_pos * 2)))))
  n_chunks <- as.integer(ceiling(n_genes / chunk))
  gmax <- numeric(n_genes)
  gscale <- rep(1, n_genes)

  for (k in seq_len(n_chunks)) {
    lo <- (k - 1L) * chunk + 1L
    hi <- min(k * chunk, n_genes)
    block <- as.matrix(expr[lo:hi, , drop = FALSE])
    mx <- apply(block, 1L, max)
    gmax[lo:hi] <- mx
    if (ds$bin) {
      # bin means are fractional: carry them as uint16 scaled to each gene's max
      sc <- ifelse(mx > 0, mx / 65535, 1)
      gscale[lo:hi] <- sc
      block <- round(block / sc)
    }
    # t() so the column-major dump lands gene-major: all positions of gene 1, then gene 2
    write_u16(t(block), file.path(dir, sprintf("expr_%03d.bin", k - 1L)))
    rm(block)
    if (k %% 20 == 0) gc(FALSE)
  }

  info[[ds$id]] <- list(
    id = ds$id, label = ds$label, binned = ds$bin,
    nSpots = n_pos, nGenes = n_genes, nSource = length(cells),
    chunkSize = chunk, nChunks = n_chunks,
    genes = genes, cellTypes = cell_types,
    exprMax = gmax, exprScale = gscale
  )
  message("  wrote ", n_chunks, " expression chunk(s), ", n_pos, " positions, ",
          n_genes, " genes")
  rm(expr, esv_mat, out_ct, out_coords); gc(FALSE)
}

# ---- shared PoPS table ------------------------------------------------------

message("== PoPS ==")
pops <- as.data.frame(readRDS(file.path(DATA_DIR, "pops_score.rds")))
stopifnot(!anyDuplicated(pops$GeneName))
diseases <- setdiff(colnames(pops), "GeneName")

gene_union <- sort(unique(unlist(lapply(info, `[[`, "genes"), use.names = FALSE)))
pops_genes <- sort(intersect(gene_union, pops$GeneName))
message("  ", length(pops_genes), " of ", length(gene_union), " gene(s) have PoPS scores")

pops_sub <- pops[match(pops_genes, pops$GeneName), diseases, drop = FALSE]
# t() so the column-major dump lands gene-major: all diseases of gene 1, then gene 2
write_f32(as.vector(t(as.matrix(pops_sub))), file.path(OUT, "pops.bin"))
json(pops_genes, file.path(OUT, "pops_genes.json"))

for (nm in names(info)) {
  d <- info[[nm]]
  json(list(
    id = d$id, label = d$label, binned = d$binned,
    nSpots = d$nSpots, nGenes = d$nGenes, nSource = d$nSource,
    chunkSize = d$chunkSize, nChunks = d$nChunks,
    genes = d$genes, cellTypes = d$cellTypes,
    exprMax = round(d$exprMax, 4), exprScale = d$exprScale,
    popsRow = match(d$genes, pops_genes) - 1L     # -1L -> NA -> null in JSON
  ), file.path(OUT, d$id, "meta.json"))
}

# Merge with whatever the manifest already lists, so rebuilding a subset does
# not drop the datasets that were not rebuilt this run.
manifest_path <- file.path(OUT, "manifest.json")
entries <- list()
if (file.exists(manifest_path)) {
  prev <- jsonlite::fromJSON(manifest_path, simplifyDataFrame = FALSE)$datasets
  for (e in prev) if (dir.exists(file.path(OUT, e$id))) entries[[e$id]] <- e
}
for (nm in names(info)) {
  d <- info[[nm]]
  entries[[d$id]] <- list(id = d$id, label = d$label, nSpots = d$nSpots,
                          nGenes = d$nGenes, binned = d$binned)
}
all_ids <- vapply(datasets, `[[`, "", "id")
known <- c(all_ids, setdiff(names(entries), all_ids))
manifest_ds <- unname(entries[intersect(known, names(entries))])
json(list(datasets = manifest_ds, diseases = diseases), manifest_path)

message("Done -> ", OUT)
