test_that(".validate_inputs accepts valid inputs", {
  d <- make_test_data()
  expect_no_error(
    spacelink:::.validate_inputs(d$counts, d$coords, NULL,
                                 lite = FALSE, grid_size = NULL,
                                 kernel = "Exponential",
                                 n_lengthscales = 5, n_workers = 1)
  )
})

test_that(".validate_inputs converts data.frame inputs to matrix", {
  d <- make_test_data()
  res <- spacelink:::.validate_inputs(
    as.data.frame(d$counts), as.data.frame(d$coords), NULL,
    lite = FALSE, grid_size = NULL, kernel = "Exponential",
    n_lengthscales = 5, n_workers = 1
  )
  expect_true(is.matrix(res$normalized_counts))
  expect_true(is.matrix(res$spatial_coords))
})

test_that(".validate_inputs rejects mismatched dimensions", {
  d <- make_test_data(n_spots = 30)
  bad_coords <- d$coords[1:20, ]
  expect_error(
    spacelink:::.validate_inputs(d$counts, bad_coords, NULL,
                                 FALSE, NULL, "Exponential", 5, 1),
    "must equal"
  )
})

test_that(".validate_inputs rejects non-numeric normalized_counts", {
  d <- make_test_data()
  char_counts <- matrix(as.character(d$counts), nrow = nrow(d$counts))
  expect_error(
    spacelink:::.validate_inputs(char_counts, d$coords, NULL,
                                 FALSE, NULL, "Exponential", 5, 1)
  )
})

test_that(".validate_inputs rejects invalid kernel", {
  d <- make_test_data()
  expect_error(
    spacelink:::.validate_inputs(d$counts, d$coords, NULL,
                                 FALSE, NULL, "Matern", 5, 1),
    "'kernel' must be"
  )
})

test_that(".validate_inputs rejects singular covariates", {
  d <- make_test_data()
  # Two identical columns → singular
  cov_mat <- cbind(rep(1, ncol(d$counts)), rep(1, ncol(d$counts)))
  expect_error(
    spacelink:::.validate_inputs(d$counts, d$coords, cov_mat,
                                 FALSE, NULL, "Exponential", 5, 1),
    "singular"
  )
})

test_that(".make_kernel returns correct values", {
  exp_k  <- spacelink:::.make_kernel("Exponential")
  gaus_k <- spacelink:::.make_kernel("Gaussian")

  # At distance 0, both kernels equal 1
  expect_equal(exp_k(1, 0),  1)
  expect_equal(gaus_k(1, 0), 1)

  # Exponential: exp(-phi * d)
  expect_equal(exp_k(0.5, 2), exp(-1))

  # Gaussian: exp(-(phi*d)^2 / 2)
  expect_equal(gaus_k(1, 2), exp(-2))

  # Kernel values are in (0, 1] for positive distances
  d_vals <- c(0.1, 1, 5, 10)
  expect_true(all(exp_k(1, d_vals)  > 0 & exp_k(1, d_vals)  <= 1))
  expect_true(all(gaus_k(1, d_vals) > 0 & gaus_k(1, d_vals) <= 1))
})

test_that(".make_kernel rejects unknown kernel name", {
  expect_error(spacelink:::.make_kernel("Matern"))
})

test_that(".make_residual_fn returns zero-mean residuals (no covariates)", {
  fn <- spacelink:::.make_residual_fn(NULL)
  y  <- rnorm(50, mean = 5)
  r  <- fn(y)
  expect_equal(mean(r), 0, tolerance = 1e-10)
})

test_that(".make_residual_fn projects out covariates", {
  set.seed(1)
  n   <- 50
  cov <- matrix(rnorm(n), ncol = 1)
  fn  <- spacelink:::.make_residual_fn(cov)
  y   <- rnorm(n)
  r   <- fn(y)

  # Residuals should be orthogonal to intercept and covariate
  expect_equal(sum(r),          0,  tolerance = 1e-8)
  expect_equal(sum(cov[, 1] * r), 0, tolerance = 1e-8)
})
