## Builds the binary payload behind the "Data browser" article.
##
## Reads the per-dataset RDS files listed in `datasets` below plus the shared
## PoPS table, and writes a static bundle to pkgdown/assets/databrowser/, which
## pkgdown copies to the site root. Re-run whenever a dataset is added or the
## source RDS files change, then rebuild the site.
##
##   Rscript data-raw/build_databrowser_data.R
##
## Layout written:
##   databrowser/manifest.json          dataset list + disease list
##   databrowser/pops/<k>.bin           float32[nPopsGenes], one file per disease
##   databrowser/pops_genes.json        gene order for the pops/*.bin files
##   databrowser/<id>/meta.json         genes, cell types, per-gene max, pops row map
##   databrowser/<id>/coords.bin        float32[nSpots*2], x,y interleaved
##   databrowser/<id>/celltype.bin      float32[nSpots*nCellTypes], spot-major
##   databrowser/<id>/esv.bin           float32[nCellTypes*nGenes], cell-type-major
##   databrowser/<id>/expr_NNN.bin      uint8[chunk*nSpots], gene-major, per-gene scaled

suppressMessages(library(Matrix))

DL <- "~/Downloads"
OUT <- "pkgdown/assets/databrowser"
CHUNK <- 64L

datasets <- list(
  list(
    id     = "visium_human_dlpfc",
    label  = "Visium human DLPFC",
    counts = file.path(DL, "Visium_human_DLPFC_counts.rds"),
    coords = file.path(DL, "Visium_human_DLPFC_spatial_coords.rds"),
    esv    = file.path(DL, "Visium_human_DLPFC_ESV_results.rds"),
    ct     = file.path(DL, "Visium_human_DLPFC_cell_type_data.rds")
  )
)
pops_file <- file.path(DL, "pops_score.rds")

# ---- helpers ----------------------------------------------------------------

write_f32 <- function(x, path) {
  con <- file(path, "wb"); on.exit(close(con))
  writeBin(as.double(x), con, size = 4L, endian = "little")
}

write_u8 <- function(x, path) {
  con <- file(path, "wb"); on.exit(close(con))
  writeBin(as.raw(x), con)
}

json <- function(x, path) {
  writeLines(jsonlite::toJSON(x, auto_unbox = TRUE, digits = 8, na = "null"), path)
}

## The ESV tables were round-tripped through a spreadsheet, which mangled some
## gene symbols (SEPT2 -> "2-Sep", LSM3 -> "LSM<nbsp>3.00"). Both forms are
## reversible; `stopifnot` below refuses to guess if a future dataset differs.
unmangle <- function(x) {
  x <- sub("^([A-Za-z]+) ([0-9]+)\\.00$", "\\1\\2", x)
  mon <- c(Jan = "JAN", Feb = "FEB", Mar = "MARCH", Apr = "APR", May = "MAY",
           Jun = "JUN", Jul = "JUL", Aug = "AUG", Sep = "SEPT", Oct = "OCT",
           Nov = "NOV", Dec = "DEC")
  m <- regmatches(x, regexec("^([0-9]+)-([A-Za-z]{3})$", x))
  for (i in seq_along(x)) {
    if (length(m[[i]]) == 3L && !is.na(mon[m[[i]][3]])) {
      x[i] <- paste0(mon[[m[[i]][3]]], m[[i]][2])
    }
  }
  x
}

# ---- shared PoPS table ------------------------------------------------------

message("Reading PoPS table ...")
pops <- as.data.frame(readRDS(pops_file))
stopifnot(!anyDuplicated(pops$GeneName))
diseases <- setdiff(colnames(pops), "GeneName")

dir.create(file.path(OUT, "pops"), recursive = TRUE, showWarnings = FALSE)

# ---- per-dataset ------------------------------------------------------------

manifest_ds <- list()
pops_gene_union <- character()

prepared <- lapply(datasets, function(ds) {
  message("Reading ", ds$label, " ...")
  counts <- readRDS(ds$counts)
  coords <- readRDS(ds$coords)
  esv <- readRDS(ds$esv)
  ct <- readRDS(ds$ct)

  spots <- colnames(counts)
  stopifnot(identical(rownames(coords), spots), identical(rownames(ct), spots))

  genes <- sort(rownames(counts))
  counts <- counts[genes, , drop = FALSE]

  esv$GeneName <- unmangle(esv$GeneName)
  cell_types <- c("whole", colnames(ct))
  stopifnot(setequal(unique(esv$Celltype), cell_types))

  # ESV as cell-type-major matrix aligned to `genes`; NA where a pair is absent.
  esv_mat <- matrix(NA_real_, nrow = length(cell_types), ncol = length(genes),
                    dimnames = list(cell_types, genes))
  for (cty in cell_types) {
    s <- esv[esv$Celltype == cty, ]
    hit <- match(genes, s$GeneName)
    if (anyNA(hit)) {
      stop("ESV table for cell type '", cty, "' is missing ",
           sum(is.na(hit)), " gene(s) present in the counts matrix, e.g. ",
           paste(head(genes[is.na(hit)], 5), collapse = ", "))
    }
    esv_mat[cty, ] <- s$ESV[hit]
  }

  list(ds = ds, counts = counts, coords = coords, ct = ct,
       genes = genes, spots = spots, cell_types = cell_types, esv_mat = esv_mat)
})

for (p in prepared) pops_gene_union <- union(pops_gene_union, p$genes)
pops_gene_union <- sort(intersect(pops_gene_union, pops$GeneName))
message("PoPS covers ", length(pops_gene_union), " gene(s) across all datasets")

pops_sub <- pops[match(pops_gene_union, pops$GeneName), , drop = FALSE]
for (k in seq_along(diseases)) {
  write_f32(pops_sub[[diseases[k]]], file.path(OUT, "pops", sprintf("%d.bin", k - 1L)))
}
json(pops_gene_union, file.path(OUT, "pops_genes.json"))

for (p in prepared) {
  ds <- p$ds
  dir <- file.path(OUT, ds$id)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  message("Writing ", ds$id, " ...")

  n_spots <- length(p$spots)
  n_genes <- length(p$genes)

  write_f32(as.vector(t(as.matrix(p$coords[, 1:2]))), file.path(dir, "coords.bin"))
  write_f32(as.vector(t(as.matrix(p$ct))), file.path(dir, "celltype.bin"))
  write_f32(as.vector(t(p$esv_mat)), file.path(dir, "esv.bin"))

  # Expression, quantised to uint8 against each gene's own maximum. The colour
  # ramp has ~7 stops, so 8-bit is well past what the plot can resolve.
  gmax <- numeric(n_genes)
  n_chunks <- ceiling(n_genes / CHUNK)
  for (k in seq_len(n_chunks)) {
    lo <- (k - 1L) * CHUNK + 1L
    hi <- min(k * CHUNK, n_genes)
    block <- as.matrix(p$counts[lo:hi, , drop = FALSE])
    mx <- apply(block, 1L, max)
    gmax[lo:hi] <- mx
    scaled <- block / ifelse(mx > 0, mx, 1)
    # t() so the column-major dump lands gene-major: all spots of gene 1, then gene 2
    write_u8(round(pmin(pmax(t(scaled), 0), 1) * 255),
             file.path(dir, sprintf("expr_%03d.bin", k - 1L)))
  }

  json(list(
    id = ds$id,
    label = ds$label,
    nSpots = n_spots,
    nGenes = n_genes,
    chunkSize = CHUNK,
    nChunks = n_chunks,
    genes = p$genes,
    cellTypes = p$cell_types,
    exprMax = round(gmax, 5),
    popsRow = match(p$genes, pops_gene_union) - 1L   # -1L -> NA -> null in JSON
  ), file.path(dir, "meta.json"))

  manifest_ds[[length(manifest_ds) + 1L]] <- list(
    id = ds$id, label = ds$label, nSpots = n_spots, nGenes = n_genes
  )
}

json(list(datasets = manifest_ds, diseases = diseases),
     file.path(OUT, "manifest.json"))

message("Done -> ", OUT)
