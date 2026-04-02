#' @importFrom RANN nn2
#' @importFrom fields rdist
#' @importFrom Matrix Matrix tril drop0
#' @importFrom methods as
#' @importFrom pracma logspace
#' @importFrom RcppML nnls
#' @importFrom stats dist pchisq
NULL

.validate_inputs <- function(normalized_counts, spatial_coords, covariates,
                              lite, grid_size, kernel, n_lengthscales, n_workers) {
  if (inherits(normalized_counts, "sparseMatrix")) {
    if (!is.numeric(normalized_counts@x))
      stop("'normalized_counts' contains non-numeric values.", call. = FALSE)
  } else if (is.data.frame(normalized_counts)) {
    if (!all(vapply(normalized_counts, is.numeric, logical(1))))
      stop("'normalized_counts' contains non-numeric columns.", call. = FALSE)
    normalized_counts <- as.matrix(normalized_counts)
  } else if (!is.matrix(normalized_counts) || !is.numeric(normalized_counts)) {
    stop("'normalized_counts' must be a numeric matrix, data.frame, or sparse matrix.", call. = FALSE)
  }

  if (is.data.frame(spatial_coords)) {
    if (!all(vapply(spatial_coords, is.numeric, logical(1))))
      stop("'spatial_coords' contains non-numeric columns.", call. = FALSE)
    if (ncol(spatial_coords) != 2)
      stop("'spatial_coords' must have 2 columns (x, y).", call. = FALSE)
    spatial_coords <- as.matrix(spatial_coords)
  } else if (!is.matrix(spatial_coords) || !is.numeric(spatial_coords)) {
    stop("'spatial_coords' must be a numeric matrix or data.frame.", call. = FALSE)
  }

  if (ncol(normalized_counts) != nrow(spatial_coords))
    stop(sprintf(
      "ncol(normalized_counts) (%d) must equal nrow(spatial_coords) (%d).",
      ncol(normalized_counts), nrow(spatial_coords)
    ), call. = FALSE)

  if (!is.null(covariates)) {
    if (is.vector(covariates) || is.factor(covariates))
      covariates <- matrix(covariates, ncol = 1)
    if (is.data.frame(covariates))
      covariates <- model.matrix(~ . - 1, data = covariates)
    if (!is.matrix(covariates) || !is.numeric(covariates))
      stop("'covariates' must be a numeric vector, matrix, or data.frame.", call. = FALSE)
    if (nrow(covariates) != nrow(spatial_coords))
      stop(sprintf(
        "nrow(covariates) (%d) must equal nrow(spatial_coords) (%d).",
        nrow(covariates), nrow(spatial_coords)
      ), call. = FALSE)
    if (rcond(crossprod(cbind(1, covariates))) < .Machine$double.eps)
      stop("'covariates' (with intercept) is singular or near-singular.", call. = FALSE)
  }

  if (!is.logical(lite))
    stop("'lite' must be logical.", call. = FALSE)
  if (!is.null(grid_size) && (!is.numeric(grid_size) || grid_size <= 0))
    stop("'grid_size' must be NULL or a positive number.", call. = FALSE)
  if (!kernel %in% c("Exponential", "Gaussian"))
    stop("'kernel' must be 'Exponential' or 'Gaussian'.", call. = FALSE)
  if (!is.numeric(n_lengthscales) || n_lengthscales < 1)
    stop("'n_lengthscales' must be a positive integer.", call. = FALSE)
  if (!is.numeric(n_workers) || n_workers < 1)
    stop("'n_workers' must be a positive integer.", call. = FALSE)

  list(normalized_counts = normalized_counts, spatial_coords = spatial_coords,
       covariates = covariates)
}

.make_residual_fn <- function(covariates) {
  if (is.null(covariates)) {
    function(y) y - mean(y)
  } else {
    X       <- cbind(1, covariates)
    XtX_inv <- solve(crossprod(X))
    M       <- diag(nrow(X)) - X %*% XtX_inv %*% t(X)
    function(y) as.numeric(M %*% y)
  }
}

.make_kernel <- function(kernel) {
  switch(kernel,
    "Exponential" = function(phi, D) exp(-phi * D),
    "Gaussian"    = function(phi, D) exp(-(phi * D)^2 / 2),
    stop("Unsupported kernel: choose 'Exponential' or 'Gaussian'", call. = FALSE)
  )
}

.build_dist_context <- function(spatial_coords, lite, grid_size) {
  N <- nrow(spatial_coords)

  if (!lite) {
    D <- dist(spatial_coords)
    return(list(type = "full", D = D, N = N, min_dist = min(D), max_dist = max(D)))
  }

  # Distance range: max from convex hull, min from nearest neighbor
  ch       <- chull(spatial_coords)
  max_dist <- max(dist(spatial_coords[ch, ]))
  nn       <- nn2(spatial_coords, k = 2)
  min_dist <- min(nn$nn.dists[, -1])
  rm(ch, nn)

  if (is.null(grid_size)) grid_size <- min_dist * 5

  # Assign spots to grid cells
  ox <- min(spatial_coords[, 1]); oy <- min(spatial_coords[, 2])
  gx <- floor((spatial_coords[, 1] - ox) / grid_size)
  gy <- floor((spatial_coords[, 2] - oy) / grid_size)
  grid_coords <- cbind(
    grid_x   = gx, grid_y   = gy,
    center_x = ox + (gx + 0.5) * grid_size,
    center_y = oy + (gy + 0.5) * grid_size
  )
  rownames(grid_coords) <- seq_len(N)

  # Group spot indices by grid cell
  grid_idx <- split(seq_len(N), list(gx, gy))
  grid_idx <- grid_idx[sapply(grid_idx, length) > 0]
  n_vec    <- sapply(grid_idx, length)

  # Unique grid centers (sorted: primary gy, secondary gx)
  grid_centers <- unique(grid_coords)
  grid_centers <- grid_centers[order(grid_centers[, "grid_x"]), ]
  grid_centers <- grid_centers[order(grid_centers[, "grid_y"]), ]

  # Sparse lower-triangle distance matrix between grid centers
  G        <- nrow(grid_centers)
  D_grid   <- dist(grid_centers[, c("center_x", "center_y")])
  D_sparse <- as(tril(Matrix(1, G, G, sparse = TRUE), k = -1), "TsparseMatrix")
  D_sparse@x <- as.vector(D_grid)
  rm(D_grid)

  # Adjacent grid pairs: store exact spot-to-spot distances
  is_near     <- D_sparse@x <= grid_size * sqrt(2)
  near_i      <- D_sparse@i[is_near] + 1L
  near_j      <- D_sparse@j[is_near] + 1L
  near_D_list <- lapply(seq_along(near_i), function(k)
    rdist(
      matrix(spatial_coords[grid_idx[[near_i[k]]], ], ncol = 2),
      matrix(spatial_coords[grid_idx[[near_j[k]]], ], ncol = 2)
    )
  )

  # Keep only far grid pairs in D_sparse
  D_sparse@x[is_near] <- 0
  D_sparse <- as(drop0(D_sparse), "TsparseMatrix")

  # Within-grid spot distances
  diag_D_list <- lapply(grid_idx, function(idx) dist(spatial_coords[idx, ]))

  list(
    type        = "lite",
    grid_idx    = grid_idx,
    D_sparse    = D_sparse,
    near_i      = near_i,
    near_j      = near_j,
    near_D_list = near_D_list,
    diag_D_list = diag_D_list,
    n_vec       = n_vec,
    N           = N,
    min_dist    = min_dist,
    max_dist    = max_dist
  )
}

.compute_AtA_Atx <- function(y, phi_seq, kernel_fun, dist_ctx) {
  y        <- as.vector(y)
  sum.y.sq <- sum(y^2)
  N        <- dist_ctx$N
  P        <- length(phi_seq)

  AtA <- matrix(0, P, P)
  Atx <- numeric(P)

  if (dist_ctx$type == "full") {
    D   <- dist_ctx$D
    yyt <- unlist(lapply(seq_len(N - 1), function(j) y[(j + 1):N] * y[j]))

    for (i in 1:P) {
      for (j in 1:P) {
        AtA[i, j] <- 2 * sum(kernel_fun(phi_seq[i] + phi_seq[j], D)) + N
      }
      Atx[i] <- 2 * sum(kernel_fun(phi_seq[i], D) * yyt) + sum.y.sq
    }

  } else {
    D_sparse    <- dist_ctx$D_sparse
    near_i      <- dist_ctx$near_i
    near_j      <- dist_ctx$near_j
    near_D_list <- dist_ctx$near_D_list
    diag_D_list <- dist_ctx$diag_D_list
    n_vec       <- dist_ctx$n_vec
    grid_idx    <- dist_ctx$grid_idx
    y_sum_vec   <- sapply(grid_idx, function(x) sum(y[x]))

    for (i in 1:P) {
      for (j in i:P) {
        phi       <- phi_seq[i] + phi_seq[j]
        AtA[i, j] <-
          2 * sum(kernel_fun(phi, D_sparse@x) * n_vec[D_sparse@i + 1] * n_vec[D_sparse@j + 1]) +
          2 * sum(sapply(near_D_list, function(x) sum(kernel_fun(phi, x)))) +
          2 * sum(sapply(diag_D_list, function(x) sum(kernel_fun(phi, x)))) + N
      }
      phi    <- phi_seq[i]
      Atx[i] <-
        2 * sum(kernel_fun(phi, D_sparse@x) * y_sum_vec[D_sparse@i + 1] * y_sum_vec[D_sparse@j + 1]) +
        2 * sum(sapply(seq_along(near_D_list), function(k)
          sum((y[grid_idx[[near_i[k]]]] %*% kernel_fun(phi, near_D_list[[k]])) *
                y[grid_idx[[near_j[k]]]]))) +
        2 * sum(sapply(seq_along(diag_D_list), function(k) {
          idx <- grid_idx[[k]]
          if (length(idx) == 1) return(0)
          pairs <- unlist(lapply(seq_len(length(idx) - 1), function(ii)
            y[idx[ii]] * y[idx[(ii + 1):length(idx)]]))
          sum(pairs * kernel_fun(phi, diag_D_list[[k]]))
        })) +
        sum.y.sq
    }
  }

  AtA <- AtA + t(AtA) - diag(diag(AtA))
  AtA <- rbind(rep(N, P), AtA)
  AtA <- cbind(rep(N, P + 1), AtA)
  Atx <- c(sum.y.sq, Atx)

  list(AtA = AtA, Atx = Atx)
}

.select_lengthscales <- function(y, dist_ctx, n_lengthscales, kernel_fun) {
  phi_seq <- 1 / logspace(log10(dist_ctx$min_dist), log10(dist_ctx$max_dist),
                          2 * n_lengthscales)

  mom   <- .compute_AtA_Atx(y, phi_seq, kernel_fun, dist_ctx)
  theta <- nnls(mom$AtA, as.matrix(mom$Atx))

  mme_res    <- data.frame(t(theta))
  colnames(mme_res) <- c("tau.sq", paste0("sigma.sq", seq_along(phi_seq)))
  sigma_vals <- mme_res[, grep("sigma.sq", colnames(mme_res)), drop = FALSE]

  phi_idx <- apply(sigma_vals > 1e-5, 1, function(x) {
    idx <- which(x)
    n   <- length(phi_seq)
    if (length(idx) == 0) {
      c(n - 1, n)
    } else if (length(idx) == 1) {
      if      (idx == 1) c(1, 2)
      else if (idx == n) c(n - 1, n)
      else               c(idx - 1, idx + 1)
    } else {
      c(idx[1], idx[length(idx)])
    }
  })

  t(apply(phi_idx, 2, function(x)
    logspace(log10(phi_seq[x[1]]), log10(phi_seq[x[2]]), n_lengthscales)))
}

.score_test <- function(N, sum.y.sq, ytSy_vec, tr_Sigma_sq_vec) {
  vapply(seq_along(ytSy_vec), function(i) {
    score  <- N / tr_Sigma_sq_vec[i] * ytSy_vec[i]
    df_val <- N^2 / tr_Sigma_sq_vec[i]
    exp(pchisq(score * N / sum.y.sq, df_val, lower.tail = FALSE, log.p = TRUE))
  }, numeric(1))
}

.aggregate_results <- function(out, n_lengthscales, gene_names) {
  results  <- do.call("rbind", lapply(out, function(x) x$mme_res))
  pval_mat <- do.call("rbind", lapply(out, function(x) x$pval_vec))
  colnames(pval_mat) <- paste0("pval", 1:n_lengthscales)
  results  <- cbind(results, pval_mat)

  suppressWarnings(results$pval <- apply(pval_mat, 1, ACAT))
  results$padj    <- p.adjust(results$pval, method = "BH")
  results$ESV_adj <- results$ESV
  results$ESV_adj[results$padj > 0.05] <- 0

  if (!is.null(gene_names)) rownames(results) <- make.unique(gene_names)
  results
}

.validate_inputs_ctSVG <- function(normalized_counts, spatial_coords,
                                    cell_type_proportions, focal_cell_type,
                                    covariates, lite, grid_size, kernel,
                                    n_lengthscales, c1, c2, n_workers) {
  if (inherits(normalized_counts, "sparseMatrix")) {
    if (!is.numeric(normalized_counts@x))
      stop("'normalized_counts' contains non-numeric values.", call. = FALSE)
  } else if (is.data.frame(normalized_counts)) {
    if (!all(vapply(normalized_counts, is.numeric, logical(1))))
      stop("'normalized_counts' contains non-numeric columns.", call. = FALSE)
    normalized_counts <- as.matrix(normalized_counts)
  } else if (!is.matrix(normalized_counts) || !is.numeric(normalized_counts)) {
    stop("'normalized_counts' must be a numeric matrix, data.frame, or sparse matrix.", call. = FALSE)
  }

  if (is.data.frame(spatial_coords)) {
    if (!all(vapply(spatial_coords, is.numeric, logical(1))))
      stop("'spatial_coords' contains non-numeric columns.", call. = FALSE)
    if (ncol(spatial_coords) != 2)
      stop("'spatial_coords' must have 2 columns (x, y).", call. = FALSE)
    spatial_coords <- as.matrix(spatial_coords)
  } else if (!is.matrix(spatial_coords) || !is.numeric(spatial_coords)) {
    stop("'spatial_coords' must be a numeric matrix or data.frame.", call. = FALSE)
  }

  if (is.data.frame(cell_type_proportions)) {
    if (!all(vapply(cell_type_proportions, is.numeric, logical(1))))
      stop("'cell_type_proportions' contains non-numeric columns.", call. = FALSE)
    cell_type_proportions <- as.matrix(cell_type_proportions)
  } else if (!is.matrix(cell_type_proportions) || !is.numeric(cell_type_proportions)) {
    stop("'cell_type_proportions' must be a numeric matrix or data.frame.", call. = FALSE)
  }

  if (ncol(normalized_counts) != nrow(spatial_coords))
    stop(sprintf(
      "ncol(normalized_counts) (%d) must equal nrow(spatial_coords) (%d).",
      ncol(normalized_counts), nrow(spatial_coords)
    ), call. = FALSE)

  if (nrow(cell_type_proportions) != nrow(spatial_coords))
    stop(sprintf(
      "nrow(cell_type_proportions) (%d) must equal nrow(spatial_coords) (%d).",
      nrow(cell_type_proportions), nrow(spatial_coords)
    ), call. = FALSE)

  if (!(focal_cell_type %in% colnames(cell_type_proportions)))
    stop("'focal_cell_type' must be one of column names of 'cell_type_proportions'.", call. = FALSE)

  if (!is.null(covariates)) {
    if (is.vector(covariates) || is.factor(covariates))
      covariates <- matrix(covariates, ncol = 1)
    if (is.data.frame(covariates))
      covariates <- model.matrix(~ . - 1, data = covariates)
    if (!is.matrix(covariates) || !is.numeric(covariates))
      stop("'covariates' must be a numeric vector, matrix, or data.frame.", call. = FALSE)
    if (nrow(covariates) != nrow(spatial_coords))
      stop(sprintf(
        "nrow(covariates) (%d) must equal nrow(spatial_coords) (%d).",
        nrow(covariates), nrow(spatial_coords)
      ), call. = FALSE)
    if (rcond(crossprod(cbind(1, covariates))) < .Machine$double.eps)
      stop("'covariates' (with intercept) is singular or near-singular.", call. = FALSE)
  }

  if (!is.logical(lite))
    stop("'lite' must be logical.", call. = FALSE)
  if (!is.null(grid_size) && (!is.numeric(grid_size) || grid_size <= 0))
    stop("'grid_size' must be NULL or a positive number.", call. = FALSE)
  if (!kernel %in% c("Exponential", "Gaussian"))
    stop("'kernel' must be 'Exponential' or 'Gaussian'.", call. = FALSE)
  if (!is.numeric(n_lengthscales) || n_lengthscales < 1)
    stop("'n_lengthscales' must be a positive integer.", call. = FALSE)
  if (!is.numeric(c1) || c1 > 1 || c1 < 0)
    stop("'c1' must be a number between 0 and 1.", call. = FALSE)
  if (!is.numeric(c2) || c2 > 1 || c2 < 0)
    stop("'c2' must be a number between 0 and 1.", call. = FALSE)
  if (!is.numeric(n_workers) || n_workers < 1)
    stop("'n_workers' must be a positive integer.", call. = FALSE)

  list(normalized_counts = normalized_counts, spatial_coords = spatial_coords,
       cell_type_proportions = cell_type_proportions, covariates = covariates)
}

.score_test_ctSVG <- function(ct_idx, null_res, y, D, kernel_fun,
                               cell_type_proportions, coloc, phi_list) {
  N        <- nrow(D)
  null_res <- as.numeric(null_res)
  null_cov <- null_res[1] * diag(N)

  rest_ct_idx <- (1:ncol(cell_type_proportions))[-ct_idx]
  temp_idx    <- 1
  for (i in rest_ct_idx) {
    if (i %in% coloc) {
      Pi.sq <- cell_type_proportions[, i] %*% t(cell_type_proportions[, i])
      for (phi in phi_list[[i]]) {
        temp_idx <- temp_idx + 1
        if (null_res[temp_idx] > 0)
          null_cov <- null_cov + null_res[temp_idx] * kernel_fun(phi, D) * Pi.sq
      }
    }
  }
  if (sum(diag(null_cov) == 0) > 0)
    diag(null_cov)[diag(null_cov) == 0] <- 1e-10

  null_cov_inv <- tryCatch(invM_arma(null_cov), error = function(e) solve(null_cov))
  rm(null_cov)

  n_kernels       <- length(phi_list[[ct_idx]])
  tr.Sigma.seq    <- numeric(n_kernels)
  tr.Sigma.sq.seq <- numeric(n_kernels)
  ytSy.seq        <- numeric(n_kernels)
  Pi.sq <- cell_type_proportions[, ct_idx] %*% t(cell_type_proportions[, ct_idx])

  for (i in seq_len(n_kernels)) {
    Sigma.prod         <- (kernel_fun(phi_list[[ct_idx]][i], D) * Pi.sq) %*% null_cov_inv
    ytSy.seq[i]        <- vMMv_arma(y, null_cov_inv, Sigma.prod, y)
    tr.Sigma.seq[i]    <- sum(diag(Sigma.prod))
    tr.Sigma.sq.seq[i] <- sum(t(Sigma.prod) * Sigma.prod)
  }
  rm(Sigma.prod, Pi.sq, D, null_cov_inv)

  pval_vec <- vapply(seq_len(n_kernels), function(i) {
    score  <- (tr.Sigma.seq[i] / tr.Sigma.sq.seq[i]) * ytSy.seq[i]
    df_val <- (tr.Sigma.seq[i]^2) / tr.Sigma.sq.seq[i]
    exp(pchisq(score, df_val, lower.tail = FALSE, log.p = TRUE))
  }, numeric(1))

  list(pval_vec = pval_vec)
}
