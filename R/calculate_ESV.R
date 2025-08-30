#' Calculate Effective Spatial Variability
#'
#' Calculates the Effective Spatial Variability (ESV) score based on variance
#' components and spatial weights.
#'
#' @param sigma.sq_vec Vector of spatial variance components.
#' @param tau.sq Nugget variance component.
#' @param spatial_coords Spatial coordinates matrix.
#' @param phi_seq Sequence of length-scale parameters.
#'
#' @return Numeric value representing the ESV score (between 0 and 1).
#' @keywords internal
calculate_ESV <- function(sigma.sq_vec, tau.sq, spatial_coords, phi_seq){
  N <- nrow(spatial_coords)
  D <- as.matrix(dist(spatial_coords))

  weight_vec <- sapply(phi_seq, function(x)sqrt(1 - (N/sum(exp(-2*x*D))))); rm(D)

  if(sum(sigma.sq_vec)+tau.sq == 0){
    return(0)
  }else{
    return(sum(sigma.sq_vec*weight_vec)/(sum(sigma.sq_vec)+tau.sq))
  }
}
