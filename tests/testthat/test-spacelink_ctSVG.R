test_that("spacelink_ctSVG returns a data.frame with correct dimensions", {
  d   <- make_ct_data(n_genes = 3, n_spots = 30, n_celltypes = 3)
  res <- spacelink_ctSVG(d$counts, d$coords, d$props,
                         focal_cell_type = "TypeA",
                         n_lengthscales = 3, n_workers = 1)

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 3)
  expect_true(all(c("time", "pval", "padj", "ESV", "ESV_adj") %in% colnames(res)))
})

test_that("spacelink_ctSVG preserves gene names as rownames", {
  d   <- make_ct_data(n_genes = 3, n_spots = 30)
  res <- spacelink_ctSVG(d$counts, d$coords, d$props,
                         focal_cell_type = "TypeA",
                         n_lengthscales = 3, n_workers = 1)
  expect_equal(rownames(res), rownames(d$counts))
})

test_that("spacelink_ctSVG returns valid p-values in [0, 1]", {
  d   <- make_ct_data(n_genes = 3, n_spots = 30)
  res <- spacelink_ctSVG(d$counts, d$coords, d$props,
                         focal_cell_type = "TypeA",
                         n_lengthscales = 3, n_workers = 1)
  expect_true(all(res$pval  >= 0 & res$pval  <= 1))
  expect_true(all(res$padj  >= 0 & res$padj  <= 1))
})

test_that("spacelink_ctSVG ESV is in [0, 1]", {
  d   <- make_ct_data(n_genes = 3, n_spots = 30)
  res <- spacelink_ctSVG(d$counts, d$coords, d$props,
                         focal_cell_type = "TypeA",
                         n_lengthscales = 3, n_workers = 1)
  expect_true(all(res$ESV     >= 0 & res$ESV     <= 1))
  expect_true(all(res$ESV_adj >= 0 & res$ESV_adj <= 1))
})

test_that("spacelink_ctSVG ESV_adj is 0 for non-significant genes", {
  d   <- make_ct_data(n_genes = 3, n_spots = 30)
  res <- spacelink_ctSVG(d$counts, d$coords, d$props,
                         focal_cell_type = "TypeA",
                         n_lengthscales = 3, n_workers = 1)
  expect_true(all(res$ESV_adj[res$padj > 0.05] == 0))
})

test_that("spacelink_ctSVG works with covariates", {
  d   <- make_ct_data(n_genes = 3, n_spots = 30)
  set.seed(99)
  cov <- matrix(rnorm(30), ncol = 1)
  res <- spacelink_ctSVG(d$counts, d$coords, d$props,
                         focal_cell_type = "TypeA",
                         covariates = cov,
                         n_lengthscales = 3, n_workers = 1)
  expect_s3_class(res, "data.frame")
})

test_that("spacelink_ctSVG falls back to spacelink when no mixing exists", {
  # All spots belong 100% to one cell type → falls back to spacelink()
  d       <- make_test_data(n_genes = 3, n_spots = 20)
  n_spots <- ncol(d$counts)
  props   <- matrix(0, nrow = n_spots, ncol = 2)
  props[, 1] <- 1
  colnames(props) <- c("TypeA", "TypeB")

  res <- spacelink_ctSVG(d$counts, d$coords, props,
                         focal_cell_type = "TypeA",
                         n_lengthscales = 3, n_workers = 1)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 3)
  # Both paths return the same 5 columns
  expect_equal(sort(colnames(res)), sort(c("time", "pval", "padj", "ESV", "ESV_adj")))
})

test_that("spacelink_ctSVG removes spots with zero total proportion", {
  d <- make_ct_data(n_genes = 3, n_spots = 30)
  # Zero out first 5 spots
  d$props[1:5, ] <- 0
  expect_no_error(
    spacelink_ctSVG(d$counts, d$coords, d$props,
                    focal_cell_type = "TypeA",
                    n_lengthscales = 3, n_workers = 1)
  )
})

test_that("spacelink_ctSVG errors when focal_cell_type not in colnames", {
  d <- make_ct_data()
  expect_error(
    spacelink_ctSVG(d$counts, d$coords, d$props,
                    focal_cell_type = "NonExistent"),
    "'focal_cell_type'"
  )
})

test_that("spacelink_ctSVG errors on mismatched spot counts", {
  d <- make_ct_data(n_spots = 30)
  expect_error(
    spacelink_ctSVG(d$counts, d$coords, d$props[1:20, ],
                    focal_cell_type = "TypeA"),
    "must equal"
  )
})

test_that("spacelink_ctSVG errors when c1 or c2 out of [0,1]", {
  d <- make_ct_data()
  expect_error(spacelink_ctSVG(d$counts, d$coords, d$props,
                               focal_cell_type = "TypeA", c1 = -0.1),
               "'c1'")
  expect_error(spacelink_ctSVG(d$counts, d$coords, d$props,
                               focal_cell_type = "TypeA", c2 = 1.5),
               "'c2'")
})
