#' Cell Type-Specific Spatial Variability Analysis
#'
#' Detects cell type-specific spatially variable genes (ctSVGs) using a
#' score test with a data-driven colocalization gating strategy to correct
#' for spatial colocalization bias.
#'
#' @param normalized_counts Normalized count matrix (genes x spots). Accepts a
#'   numeric matrix, data.frame, or sparse matrix (\code{sparseMatrix}).
#' @param spatial_coords Spatial coordinate matrix (spots x 2, x and y columns).
#'   Accepts a numeric matrix or data.frame.
#' @param cell_type_proportions Cell type proportion matrix (spots x cell types,
#'   e.g., estimated by RCTD). Accepts a numeric matrix or data.frame.
#' @param focal_cell_type Name of the cell type to be tested. Must match a
#'   column name of \code{cell_type_proportions}.
#' @param covariates Optional covariate matrix (spots x covariates, without
#'   intercept). Default is \code{NULL}.
#' @param lite Logical. If \code{TRUE}, uses a grid-based approximation.
#'   Default is \code{FALSE}.
#' @param grid_size Grid cell size for lite mode. If \code{NULL}, set
#'   automatically. Default is \code{NULL}.
#' @param kernel Spatial kernel type. One of \code{"Exponential"} (default) or
#'   \code{"Gaussian"}.
#' @param n_lengthscales Number of length scales (kernels). Default is 5.
#' @param c1 Gating parameter controlling the proportion threshold for
#'   determining if the focal cell type dominates. Default is 0.1.
#' @param c2 Gating parameter controlling the colocalization threshold.
#'   Default is 0.1.
#' @param n_workers Number of workers for parallelization. Default is 1.
#'
#' @return A data frame (genes x results) containing:
#' \describe{
#'   \item{time}{Computation time per gene (seconds)}
#'   \item{pval}{Combined p-value (ACAT) for cell type-specific spatial patterns}
#'   \item{padj}{Benjamini-Hochberg adjusted p-value}
#'   \item{ESV}{Cell type-specific Effective Spatial Variability score}
#'   \item{ESV_adj}{ESV adjusted for multiple testing (0 if padj > 0.05)}
#' }
#' If all spots are pure (no mixing), falls back to \code{\link{spacelink}}
#' output (same columns) for the focal cell type spots.
#'
#' @examples
#' library(spacelink)
#' set.seed(123)
#' n_spots    <- 100
#' n_genes    <- 10
#' n_celltypes <- 2
#'
#' coords <- expand.grid(
#'   x = seq(0, 10, length.out = sqrt(n_spots)),
#'   y = seq(0, 10, length.out = sqrt(n_spots))
#' )
#'
#' cell_type_props <- matrix(nrow = n_spots, ncol = n_celltypes)
#' cell_type_props[1:(n_spots / 2), 1]            <- runif(n_spots / 2, 0, 0.5)
#' cell_type_props[(n_spots / 2 + 1):n_spots, 1]  <- runif(n_spots / 2, 0.5, 1)
#' cell_type_props[, 2] <- 1 - cell_type_props[, 1]
#' colnames(cell_type_props) <- paste0("Cell_type_", 1:n_celltypes)
#'
#' expr_data <- matrix(nrow = n_genes, ncol = n_spots)
#' rownames(expr_data) <- paste0("Gene_", 1:n_genes)
#' D <- as.matrix(dist(coords)); K <- exp(-D)
#' Sigma <- chol(K) %*% diag(cell_type_props[, "Cell_type_1"])
#' for (i in 1:(n_genes / 2))
#'   expr_data[i, ] <- rnorm(n_spots) + i * rnorm(n_spots) %*% Sigma
#' Sigma <- chol(K) %*% diag(cell_type_props[, "Cell_type_2"])
#' for (i in (n_genes / 2 + 1):n_genes)
#'   expr_data[i, ] <- rnorm(n_spots) + rnorm(n_spots) %*% Sigma
#'
#' results <- spacelink_ctSVG(
#'   normalized_counts     = expr_data,
#'   spatial_coords        = coords,
#'   cell_type_proportions = cell_type_props,
#'   focal_cell_type       = "Cell_type_1"
#' )
#' print(results[, c("pval", "ESV")])
#'
#' @importFrom ACAT ACAT
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom gaston lmm.aireml
#' @importFrom stats p.adjust
#' @export
spacelink_ctSVG <- function(normalized_counts,
                            spatial_coords,
                            cell_type_proportions,
                            focal_cell_type,
                            covariates = NULL,
                            lite = FALSE,
                            grid_size = NULL,
                            kernel = "Exponential",
                            n_lengthscales = 5,
                            c1 = 0.1,
                            c2 = 0.1,
                            n_workers = 1) {
  validated <- .validate_inputs_ctSVG(normalized_counts, spatial_coords,
                                      cell_type_proportions, focal_cell_type,
                                      covariates, lite, grid_size, kernel,
                                      n_lengthscales, c1, c2, n_workers)
  normalized_counts     <- validated$normalized_counts
  spatial_coords        <- validated$spatial_coords
  cell_type_proportions <- validated$cell_type_proportions
  covariates            <- validated$covariates
  rm(validated)

  cell_type_proportions[is.na(cell_type_proportions)] <- 0
  zero_spot_inds <- (apply(cell_type_proportions, 1, sum) == 0)
  if (sum(zero_spot_inds) > 0) {
    normalized_counts <- normalized_counts[, !zero_spot_inds]
    if (is.vector(normalized_counts)) normalized_counts <- matrix(normalized_counts,nrow=1)
  }
  spatial_coords        <- spatial_coords[!zero_spot_inds, ]
  cell_type_proportions <- cell_type_proportions[!zero_spot_inds, ]
  cell_type_proportions <- as.matrix(cell_type_proportions / apply(cell_type_proportions, 1, sum))
  if (!is.null(covariates)) covariates <- covariates[!zero_spot_inds, ]

  if (sum((cell_type_proportions > 1e-7) * (cell_type_proportions < 1 - (1e-7))) == 0) {
    ct_idx      <- which(colnames(cell_type_proportions) == focal_cell_type)
    focal_spots <- (cell_type_proportions[, ct_idx] >= 1 - (1e-7))
    normalized_counts <- normalized_counts[, focal_spots]
    if (is.vector(normalized_counts)) normalized_counts <- matrix(normalized_counts,nrow=1)
    spatial_coords    <- spatial_coords[focal_spots, ]
    if (!is.null(covariates)) covariates <- covariates[focal_spots, ]
    results <- spacelink(normalized_counts, spatial_coords, covariates,
                         lite, grid_size, kernel, n_lengthscales, n_workers)
    results <- results[, c("time", "pval", "padj", "ESV", "ESV_adj")]
    return(results)
  }

  residual_fn <- .make_residual_fn(covariates)
  kernel_fun  <- .make_kernel(kernel)
  dist_ctx    <- .build_dist_context(spatial_coords, lite, grid_size)

  ct_idx <- which(colnames(cell_type_proportions) == focal_cell_type)
  coloc  <- c()
  for (i in 1:ncol(cell_type_proportions)) {
    if (i == ct_idx) next
    q <- sum((cell_type_proportions[, i] > cell_type_proportions[, ct_idx]) &
               (cell_type_proportions[, ct_idx] > c1)) /
      sum(cell_type_proportions[, ct_idx] > c1)
    if (q > c2) coloc <- c(coloc, i)
  }

  out <- bplapply(1:nrow(normalized_counts), function(gene_idx) {
    runtime <- system.time({
      y <- as.numeric(normalized_counts[gene_idx, ])

      phi_seq <- .select_lengthscales(
        residual_fn(as.numeric(normalized_counts[gene_idx, ])),
        dist_ctx, n_lengthscales, kernel_fun
      )
      colnames(phi_seq) <- paste0("phi", 1:n_lengthscales)

      phi_list <- NULL
      for (i in 1:ncol(cell_type_proportions)) {
        phi_list[[i]] <- if (i == ct_idx) phi_seq else median(phi_seq)
      }

      # Null model REML estimation
      D          <- as.matrix(dist(spatial_coords))
      Sigma.list <- NULL
      rest_ct_idx <- (1:ncol(cell_type_proportions))[-ct_idx]
      temp_idx   <- 1
      for (i in rest_ct_idx) {
        if (i %in% coloc) {
          Pi.sq <- cell_type_proportions[, i] %*% t(cell_type_proportions[, i])
          for (phi in phi_list[[i]]) {
            Sigma.list[[temp_idx]] <- kernel_fun(phi, D) * Pi.sq
            temp_idx <- temp_idx + 1
          }
        }
      }
      if (is.null(Sigma.list)) {
        local_cov <- if (is.null(covariates)) cell_type_proportions else
          cbind(covariates, cell_type_proportions)
        mu_hat   <- local_cov %*% solve(t(local_cov) %*% local_cov, t(local_cov) %*% y)
        null_res <- mean((y - mu_hat)^2)
      } else {
        reml_model <- lmm.aireml(
          Y       = y,
          X       = cbind(rep(1, nrow(cell_type_proportions)), covariates, cell_type_proportions),
          K       = Sigma.list,
          min_s2  = 0,
          verbose = FALSE
        )
        mu_hat   <- cbind(rep(1, nrow(cell_type_proportions)), covariates, cell_type_proportions) %*%
          reml_model$BLUP_beta
        null_res <- c(reml_model$sigma2, reml_model$tau)
      }
      rm(Sigma.list)

      test_res <- .score_test_ctSVG(ct_idx, null_res, y - mu_hat, D,
                                    kernel_fun, cell_type_proportions, coloc, phi_list)

      N          <- nrow(spatial_coords)
      weight_vec <- sapply(phi_list[[ct_idx]],
                           function(x) sqrt(max(0, 1 - (N / sum(exp(-2 * x * D))))))

      # Full model REML estimation
      Sigma.list <- NULL
      temp_idx   <- 1
      Pi.sq      <- cell_type_proportions[, ct_idx] %*% t(cell_type_proportions[, ct_idx])
      for (phi in phi_list[[ct_idx]]) {
        Sigma.list[[temp_idx]] <- kernel_fun(phi, D) * Pi.sq
        temp_idx <- temp_idx + 1
      }
      for (i in rest_ct_idx) {
        if (i %in% coloc) {
          Pi.sq <- cell_type_proportions[, i] %*% t(cell_type_proportions[, i])
          for (phi in phi_list[[i]]) {
            Sigma.list[[temp_idx]] <- kernel_fun(phi, D) * Pi.sq
            temp_idx <- temp_idx + 1
          }
        }
      }
      reml_model <- lmm.aireml(
        Y       = y,
        X       = cbind(rep(1, nrow(cell_type_proportions)), covariates, cell_type_proportions),
        K       = Sigma.list,
        min_s2  = 0,
        verbose = FALSE
      )

      numerator   <- sum(reml_model$tau[1:length(phi_list[[ct_idx]])] * weight_vec)
      denominator <- sum(reml_model$tau[1:length(phi_list[[ct_idx]])]) + reml_model$sigma2
      ESV         <- if (denominator < 1e-5) 0 else numerator / denominator
    })
    list(time = runtime[["elapsed"]], pval_vec = test_res$pval_vec, ESV = ESV)
  }, BPPARAM = MulticoreParam(workers = n_workers))

  results  <- data.frame(time = do.call("rbind", lapply(out, function(x) x$time)))
  pval_mat <- do.call("rbind", lapply(out, function(x) x$pval_vec))
  suppressWarnings({ results$pval <- apply(pval_mat, 1, ACAT) })
  results$padj    <- p.adjust(results$pval, method = "BH")
  results$ESV     <- do.call("rbind", lapply(out, function(x) x$ESV))
  results$ESV_adj <- results$ESV
  results$ESV_adj[results$padj > 0.05] <- 0

  if (!is.null(rownames(normalized_counts)))
    rownames(results) <- make.unique(rownames(normalized_counts))

  results
}
