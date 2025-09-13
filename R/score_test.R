#' Score Test for Spatial Gene Expression
#'
#' Performs score tests to assess spatial gene expression patterns across
#' multiple length scales.
#'
#' @param y Gene expression vector.
#' @param spatial_coords Spatial coordinates matrix.
#' @param phi_seq Sequence of length-scale parameters.
#'
#' @return List containing p-values for each length scale.
#' @importFrom stats dist pchisq
#' @keywords internal
score_test <- function(y, spatial_coords, phi_seq){
  if(is.data.frame(phi_seq)){
    phi_seq <- as.matrix(phi_seq)
  }
  n_kernels <- length(phi_seq)

  sum.y.sq <- sum(y^2)

  if(sum.y.sq == 0){
    return(list(pval_vec=rep(0, n_kernels)))
  }

  N <- nrow(spatial_coords)
  D <- as.matrix(dist(spatial_coords))

  tr.Sigma.seq <- rep(N, n_kernels)
  tr.Sigma.sq.seq <- rep(NA, n_kernels)
  ytSy.seq <- rep(NA, n_kernels)
  for(i in 1:length(phi_seq)){
    Sigma <- exp(-phi_seq[i]*D)
    tr.Sigma.sq.seq[i] <- norm(Sigma, "F")^2
    ytSy.seq[i] <- vMv_arma(y, Sigma, y)
  }
  rm(Sigma, D)

  pval_vec <- rep(NA, n_kernels)
  for(i in 1:n_kernels){
    score <- (tr.Sigma.seq[i]/tr.Sigma.sq.seq[i])*ytSy.seq[i]
    df_val <- (tr.Sigma.seq[i]^2)/tr.Sigma.sq.seq[i]
    pval_vec[i] <- exp(pchisq(score*N/sum.y.sq, df_val, lower.tail=FALSE, log.p=TRUE))
  }
  return(list(pval_vec=pval_vec))
}
