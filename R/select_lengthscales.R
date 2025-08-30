#################################################
#' @param Y Normalized count matrix. Rows and columns indicate genes and spots, respectively.
#' @param spatial_coords Spatial coordinates matrix. Rows indicate spots.
#' @param n_lengthscales Number of length-scales (i.e., number of kernels)
#' @param M It controls the minimum length-scale. Minimum length-scale is set to be the minimum distance multiplied by M.
#' @param is.inverse Is an inverse of lengthscale returned?
#' @param n_workers Number of workers for parallelization.
#################################################

select_lengthscales <- function(Y, spatial_coords, n_lengthscales = 5, M = 1, is.inverse = TRUE, n_workers = 1){
  D <- as.matrix(dist(spatial_coords))
  diag(D) <- NA
  D <- D[!is.na(D)]
  phi_seq <- 1/logspace(log10(M*min(D)), log10(max(D)), 2*n_lengthscales)
  
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
                     function(x)logspace(log10(phi_seq[x[1]]), log10(phi_seq[x[2]]), n_lengthscales)))
  
  if(is.inverse){
    return(list(phi_mat=phi_mat, time=nnls_res$time))
  }else{
    return(list(phi_mat=1/phi_mat, time=nnls_res$time))
  }
}