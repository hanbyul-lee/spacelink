# Shared synthetic data for all tests
# Ensures the package (and its compiled C++ code) is loaded when running
# testthat::test_file() directly without devtools::test().
if (!isNamespaceLoaded("spacelink")) {
  pkg_root <- normalizePath(file.path(dirname(getwd()), ".."), mustWork = FALSE)
  if (file.exists(file.path(pkg_root, "DESCRIPTION"))) {
    pkgload::load_all(pkg_root, quiet = TRUE)
  } else {
    library(spacelink)
  }
}

make_test_data <- function(n_genes = 5, n_spots = 30, seed = 42) {
  set.seed(seed)

  # Genes x Spots count matrix
  counts <- matrix(
    rpois(n_genes * n_spots, lambda = 10),
    nrow = n_genes, ncol = n_spots
  )
  rownames(counts) <- paste0("gene", seq_len(n_genes))

  # Spots x 2 spatial coordinates on a grid with slight jitter
  grid_n <- ceiling(sqrt(n_spots))
  coords <- expand.grid(x = seq_len(grid_n), y = seq_len(grid_n))[seq_len(n_spots), ]
  coords <- as.matrix(coords) + matrix(rnorm(n_spots * 2, sd = 0.1), ncol = 2)

  list(counts = counts, coords = coords)
}

make_ct_data <- function(n_genes = 5, n_spots = 30, n_celltypes = 3, seed = 42) {
  base <- make_test_data(n_genes, n_spots, seed)

  set.seed(seed + 1)
  props <- matrix(runif(n_spots * n_celltypes), nrow = n_spots, ncol = n_celltypes)
  props <- props / rowSums(props)
  colnames(props) <- paste0("Type", LETTERS[seq_len(n_celltypes)])

  c(base, list(props = props))
}
