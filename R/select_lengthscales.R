#' Select Optimal Length Scales for Spatial Analysis
#'
#' Automatically selects optimal length scales (kernel bandwidths) for spatial
#' gene expression analysis using a two-step approach based on non-negative
#' least squares fitting.
#'
#' @param Y Normalized count matrix. Rows and columns indicate genes and spots, respectively.
#' @param spatial_coords Spatial coordinates matrix. Rows indicate spots.
#' @param n_lengthscales Number of length-scales (i.e., number of kernels). Default is 5.
#' @param M Controls the minimum length-scale. Minimum length-scale is set to be
#'   the minimum distance multiplied by M. Default is 1.
#' @param is.inverse Should the inverse of lengthscale (phi) be returned? Default is TRUE.
#' @param n_workers Number of workers for parallelization. Default is 1.
#'
#' @details
#' This function implements a two-step approach for selecting optimal length scales:
#' 1. First, it fits models across a broad range of length scales
#' 2. Then, for each gene, it identifies the range where variance components are non-zero
#' 3. Finally, it selects gene-specific length scales within those ranges
#'
#' The function uses \code{\link{run_nnls}} internally to fit non-negative least
#' squares models for each gene across different length scales.
#'
#' @return A list containing:
#' \describe{
#'   \item{phi_mat}{Matrix of selected length scale parameters (phi = 1/lengthscale)
#'                  for each gene. Each row corresponds to a gene.}
#'   \item{time}{Vector of computation times for each gene.}
#' }
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' n_spots <- 100
#' n_genes <- 20
#' coords <- matrix(runif(n_spots * 2), ncol = 2)
#' expr_data <- matrix(rnorm(n_genes * n_spots), nrow = n_genes)
#'
#' # Select optimal length scales
#' length_scales <- select_lengthscales(
#'   Y = expr_data,
#'   spatial_coords = coords,
#'   n_lengthscales = 5
#' )
#'
#' # View selected scales for first few genes
#' head(length_scales$phi_mat)
#' }
#'
#' @seealso \code{\link{run_nnls}}, \code{\link{spacelink_global}}
#' @keywords internal
select_lengthscales <- function(Y, spatial_coords, n_lengthscales = 5, M = 1, is.inverse = TRUE, n_workers = 1){
  D <- as.matrix(dist(spatial_coords))
  diag(D) <- NA
  D <- D[!is.na(D)]
  phi_seq <- 1/pracma::logspace(log10(M*min(D)), log10(max(D)), 2*n_lengthscales)

  Y <- Y - rowMeans(Y)
  D <- as.matrix(dist(spatial_coords))

  out <- bplapply(1:nrow(Y), function(gene_idx) {
    runtime <- system.time({
      y <- as.numeric(Y[gene_idx,])
      res <- run_nnls(y, spatial_coords, phi_seq = phi_seq)
    })
    res$time <- runtime[["elapsed"]]
    return(res)
  }, BPPARAM = MulticoreParam(workers = n_workers))

  nnls_res <- do.call("rbind", out)

  phi_idx <- apply(nnls_res[,grep("sigma.sq", colnames(nnls_res))] > 1e-5, 1,
                   function(x){idx <- which(x);
                   if(length(idx)==0){
                     c(length(phi_seq)-1, length(phi_seq))
                   }else if(length(idx)==1){
                     if(idx[1]==1){
                       c(idx[1], idx[1]+1)
                     }else if(idx[1]==length(phi_seq)){
                       c(idx[1]-1, idx[1])
                     }else{
                       c(idx[1]-1, idx[1]+1)
                     }
                   }else{c(idx[1], idx[length(idx)])}})

  phi_mat <- t(apply(phi_idx, 2,
                     function(x)pracma::logspace(log10(phi_seq[x[1]]), log10(phi_seq[x[2]]), n_lengthscales)))

  if(is.inverse){
    return(list(phi_mat=phi_mat, time=nnls_res$time))
  }else{
    return(list(phi_mat=1/phi_mat, time=nnls_res$time))
  }
}
