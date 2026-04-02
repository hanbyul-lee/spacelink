test_that("spacelink returns a data.frame with correct dimensions (full mode)", {
  d   <- make_test_data(n_genes = 4, n_spots = 25)
  res <- spacelink(d$counts, d$coords, n_lengthscales = 3, n_workers = 1)

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 4)
  expect_true(all(c("tau.sq", "ESV", "pval", "padj", "ESV_adj") %in% colnames(res)))
})

test_that("spacelink preserves gene names as rownames", {
  d   <- make_test_data(n_genes = 3, n_spots = 25)
  res <- spacelink(d$counts, d$coords, n_lengthscales = 3, n_workers = 1)
  expect_equal(rownames(res), rownames(d$counts))
})

test_that("spacelink returns valid p-values in [0, 1]", {
  d   <- make_test_data(n_genes = 4, n_spots = 25)
  res <- spacelink(d$counts, d$coords, n_lengthscales = 3, n_workers = 1)

  expect_true(all(res$pval  >= 0 & res$pval  <= 1))
  expect_true(all(res$padj  >= 0 & res$padj  <= 1))
})

test_that("spacelink ESV is in [0, 1]", {
  d   <- make_test_data(n_genes = 4, n_spots = 25)
  res <- spacelink(d$counts, d$coords, n_lengthscales = 3, n_workers = 1)

  expect_true(all(res$ESV     >= 0 & res$ESV     <= 1))
  expect_true(all(res$ESV_adj >= 0 & res$ESV_adj <= 1))
})

test_that("spacelink works with Gaussian kernel", {
  d   <- make_test_data(n_genes = 3, n_spots = 25)
  res <- spacelink(d$counts, d$coords, kernel = "Gaussian",
                   n_lengthscales = 3, n_workers = 1)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 3)
})

test_that("spacelink works in lite mode", {
  d   <- make_test_data(n_genes = 3, n_spots = 30)
  res <- spacelink(d$counts, d$coords, lite = TRUE,
                   n_lengthscales = 3, n_workers = 1)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 3)
})

test_that("spacelink works with covariates", {
  d   <- make_test_data(n_genes = 3, n_spots = 30)
  set.seed(10)
  cov <- matrix(rnorm(30), ncol = 1)
  res <- spacelink(d$counts, d$coords, covariates = cov,
                   n_lengthscales = 3, n_workers = 1)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 3)
})

test_that("spacelink works with sparse matrix input", {
  skip_if_not_installed("Matrix")
  d      <- make_test_data(n_genes = 3, n_spots = 25)
  sparse <- Matrix::Matrix(d$counts, sparse = TRUE)
  res    <- spacelink(sparse, d$coords, n_lengthscales = 3, n_workers = 1)
  expect_s3_class(res, "data.frame")
})

test_that("spacelink works with data.frame inputs", {
  d   <- make_test_data(n_genes = 3, n_spots = 25)
  res <- spacelink(as.data.frame(d$counts), as.data.frame(d$coords),
                   n_lengthscales = 3, n_workers = 1)
  expect_s3_class(res, "data.frame")
})

test_that("spacelink handles unnamed genes gracefully", {
  d   <- make_test_data(n_genes = 3, n_spots = 25)
  rownames(d$counts) <- NULL
  res <- spacelink(d$counts, d$coords, n_lengthscales = 3, n_workers = 1)
  expect_s3_class(res, "data.frame")
})

test_that("spacelink errors on mismatched dimensions", {
  d <- make_test_data(n_genes = 3, n_spots = 25)
  expect_error(spacelink(d$counts, d$coords[1:10, ], n_lengthscales = 3))
})

test_that("spacelink errors on invalid kernel", {
  d <- make_test_data()
  expect_error(spacelink(d$counts, d$coords, kernel = "Matern"))
})
