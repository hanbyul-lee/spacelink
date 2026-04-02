#' Global Spatial Variability Analysis
#'
#' Conducts hypothesis testing to identify spatially variable genes (SVGs) using an adaptive multi-kernel approach,
#' and provides Effective Spatial Variability (ESV) scores (ranging from 0 to 1), where higher values indicate greater spatial variability in gene expression.
#'
#' @param normalized_counts Normalized count matrix (genes x spots). Accepts a
#'   numeric matrix, data.frame, or sparse matrix (\code{sparseMatrix}).
#' @param spatial_coords Two-dimensional spatial coordinate matrix (spots x 2).
#'   Accepts a numeric matrix or data.frame.
#' @param covariates Optional covariate matrix (spots x covariates, without
#'   intercept). Accepts a numeric vector, matrix, or data.frame. Default is
#'   \code{NULL}.
#' @param lite Logical. If \code{TRUE}, uses a grid-based approximation to
#'   speed up computation for large datasets. Default is \code{FALSE}.
#' @param grid_size Grid cell size for lite mode. If \code{NULL}, set
#'   automatically to 5 times the minimum nearest-neighbor distance.
#'   Default is \code{NULL}.
#' @param kernel Spatial kernel type. One of \code{"Exponential"} (default) or
#'   \code{"Gaussian"}.
#' @param n_lengthscales Number of length scales (kernels). Default is 5.
#' @param n_workers Number of workers for parallelization. Default is 1.
#'
#' @return A data frame (genes x results) containing:
#' \describe{
#'   \item{tau.sq}{Nugget variance component}
#'   \item{sigma.sq1, sigma.sq2, ...}{Spatial variance components for each kernel}
#'   \item{phi1, phi2, ...}{Inverse length scales for each kernel}
#'   \item{ESV}{Effective Spatial Variability score}
#'   \item{pval1, pval2, ...}{Score test p-values for each kernel}
#'   \item{pval}{Combined p-value (ACAT)}
#'   \item{padj}{Benjamini-Hochberg adjusted p-value}
#'   \item{ESV_adj}{ESV adjusted for non-SVGs (0 if padj > 0.05)}
#'   \item{time}{Computation time per gene (seconds)}
#' }
#'
#' @examples
#' library(spacelink)
#' set.seed(123)
#' n_spots <- 100
#' n_genes <- 10
#'
#' coords <- expand.grid(
#'   x = seq(0, 10, length.out = sqrt(n_spots)),
#'   y = seq(0, 10, length.out = sqrt(n_spots))
#' )
#'
#' expr_data <- matrix(nrow = n_genes, ncol = n_spots)
#' rownames(expr_data) <- paste0("Gene_", 1:n_genes)
#' D <- as.matrix(dist(coords))
#' K <- exp(-D / 2); chol_K <- chol(K)
#' for (i in 1:(n_genes / 2)) {
#'   expr_data[i, ] <- rnorm(n_spots) + i * rnorm(n_spots) %*% chol_K
#' }
#' for (i in (n_genes / 2 + 1):n_genes) {
#'   expr_data[i, ] <- rnorm(n_spots)
#' }
#'
#' results <- spacelink(normalized_counts = expr_data, spatial_coords = coords)
#' print(results[, c("ESV", "pval", "padj")])
#'
#' @importFrom ACAT ACAT
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom RcppML nnls
#' @importFrom stats p.adjust
#' @export
spacelink <- function(normalized_counts,
                      spatial_coords,
                      covariates = NULL,
                      lite = FALSE,
                      grid_size = NULL,
                      kernel = "Exponential",
                      n_lengthscales = 5,
                      n_workers = 1) {
  validated <- .validate_inputs(normalized_counts, spatial_coords, covariates,
                                lite, grid_size, kernel, n_lengthscales, n_workers)
  normalized_counts <- validated$normalized_counts
  spatial_coords    <- validated$spatial_coords
  covariates        <- validated$covariates
  rm(validated)

  residual_fn <- .make_residual_fn(covariates)
  N           <- nrow(spatial_coords)
  kernel_fun  <- .make_kernel(kernel)
  dist_ctx    <- .build_dist_context(spatial_coords, lite, grid_size)

  out <- bplapply(1:nrow(normalized_counts), function(gene_idx) {
    runtime <- system.time({
      y        <- residual_fn(as.numeric(normalized_counts[gene_idx, ]))
      sum.y.sq <- sum(y^2)

      phi_seq <- .select_lengthscales(y, dist_ctx, n_lengthscales, kernel_fun)
      colnames(phi_seq) <- paste0("phi", 1:n_lengthscales)

      mom <- .compute_AtA_Atx(y, phi_seq, kernel_fun, dist_ctx)
      AtA <- mom$AtA
      Atx <- mom$Atx

      theta   <- nnls(AtA, as.matrix(Atx))
      mme_res <- data.frame(t(theta))
      colnames(mme_res) <- c("tau.sq", paste0("sigma.sq", 1:n_lengthscales))
      mme_res <- cbind(mme_res, phi_seq)

      sigma.sq_vec <- unlist(mme_res[, grep("sigma", colnames(mme_res))])
      tau.sq       <- unlist(mme_res[, grep("tau",   colnames(mme_res))])
      weight_vec   <- sqrt((diag(AtA)[-1] - N) / diag(AtA)[-1])
      mme_res$ESV  <- if (sum(sigma.sq_vec) + tau.sq == 0) 0 else
        sum(sigma.sq_vec * weight_vec) / (sum(sigma.sq_vec) + tau.sq)

      pval_vec <- .score_test(N, sum.y.sq, Atx[-1], diag(AtA)[-1])
    })
    mme_res$time <- runtime[["elapsed"]]
    list(mme_res = mme_res, pval_vec = pval_vec)
  }, BPPARAM = MulticoreParam(workers = n_workers))

  .aggregate_results(out, n_lengthscales, rownames(normalized_counts))
}
