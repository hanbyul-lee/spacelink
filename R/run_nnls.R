#' Non-Negative Least Squares for Spatial Kernels
#'
#' Performs non-negative least squares fitting for spatial kernel regression.
#'
#' @param y Numeric vector of gene expression values.
#' @param spatial_coords Matrix of spatial coordinates.
#' @param phi_seq Numeric vector of length-scale parameters.
#' @param is.sparse Logical indicating whether to use sparse regularization.
#' @param lambda Regularization parameter for L1 constraint.
#' @param is.obj Logical indicating whether to return objective function value.
#'
#' @return Data frame with variance components and optional objective value.
#'
#' @importFrom RcppML nnls
#' @export
# lambda: tuning parameter for L1 constraint, used only if is.sparse = TRUE
run_nnls <- function(y, spatial_coords, phi_seq, is.sparse = FALSE, lambda = 0, is.obj = FALSE){

  if(is.null(phi_seq)){
    stop("phi_seq is required.\n")
  }

  # Construct matrix A
  N <- nrow(spatial_coords)
  A <- as.vector(diag(N))
  D <- as.matrix(dist(spatial_coords))
  for(phi in phi_seq){
    K <- exp(-phi*D)
    A <- cbind(A, as.vector(K))
  }

  # Run nnls
  x <- as.vector(y %*% t(y))
  if(is.sparse){
    ntheta <- ncol(A)
    AtA <- t(A) %*% A
    x.sparse <- x - 0.5*lambda*(A %*% solve(AtA, c(0, rep(1,ntheta-1))))
    theta <- RcppML::nnls(crossprod(A), crossprod(A,x.sparse))
  }else{
    theta <- RcppML::nnls(crossprod(A), crossprod(A,x))
  }

  res <- data.frame(t(theta))
  colnames(res) <- c("tau.sq", paste0("sigma.sq",1:length(phi_seq)))

  if(is.obj){
    # Calculate objective
    obj <- sum((A %*% theta - x)^2)
    if(is.sparse){
      obj <- obj + lambda*sum(theta[-1])
    }
    res$obj <- obj
  }

  return(res)
}
