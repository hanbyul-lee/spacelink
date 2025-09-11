#' Cell Type-Specific Spatial Variability Analysis
#'
#' Performs cell type-specific spatial gene variability analysis
#'
#' @param normalized_counts Normalized count matrix. Rows and columns indicate genes and spots, respectively.
#' @param spatial_coords Spatial coordinates matrix. Rows indicate spots.
#' @param cell_type_proportions Cell-type data matrix. Cell-type proportions (e.g., RCTD) or cell-type abundance (e.g., Cell2location).
#' @param focal_cell_type Name of cell-type to be tested.
#' @param covariates Spatial features (without intercept)
#' @param global_spacelink_results Spacelink global SVG test results from \code{\link{spacelink_global}}.
#' @param lengthscales If NULL, length-scales are calculated by using our two-step approach.
#' @param n_lengthscales Number of length-scales (i.e., number of kernels)
#' @param M It controls the minimum length-scale. Minimum length-scale is set to be the minimum distance multiplied by M.
#' @param calculate_ESV If TRUE, ct-ESV is returned.
#' @param c1 Gating parameter 1. Default is 0.25.
#' @param c2 Gating parameter 2. Default is 0.2.
#' @param n_workers Number of workers for parallelization. Default is 1.
#'
#' @return A data frame containing:
#' \describe{
#'   \item{time}{Computation time}
#'   \item{pval}{Combined p-value for cell type-specific spatial patterns}
#'   \item{ESV}{ct-ESV for cell type-specific spatial patterns}
#' }
#'
#' @examples
#' \dontrun{
#' # First run global analysis
#' global_results <- spacelink_global(normalized_counts, spatial_coords)
#'
#' # Then cell type-specific analysis
#' ct_results <- spacelink_cell_type(
#'   normalized_counts = normalized_counts,
#'   spatial_coords = spatial_coords,
#'   cell_proportions = cell_proportions,
#'   focal_cell_type = 'cell type name',
#'   global_spacelink_results = global_results
#' )
#' }
#'
#' @export
spacelink_cell_type <- function(normalized_counts, spatial_coords, cell_type_proportions, focal_cell_type, covariates = NULL, 
                                global_spacelink_results = NULL, lengthscales = NULL, n_lengthscales = 5, M = 1,
                                calculate_ESV = TRUE,
                                c1 = 0.25, c2 = 0.2, n_workers = 1){
  
  Y <- normalized_counts
  X <- covariates
  
  if(!is.null(global_spacelink_results)){
    if(nrow(Y)!=nrow(global_results)){
      stop("The row dimensions of normalized_counts and global_spacelink_results should be the same.")
    }
    phi_mat = global_results[,grep("phi",colnames(global_results))]
    n_lengthscales <- ncol(phi_mat)
  }else{
    if(!is.null(lengthscales)){
      if(nrow(lengthscales)!=nrow(Y)){
        stop("The row dimensions of normalized counts and lengthscales should be the same.")
      }
      phi_mat <- 1/lengthscales
      n_lengthscales <- ncol(phi_mat)
    }else{
      out <- select_lengthscales(Y, spatial_coords, n_lengthscales, M = M, is.inverse = TRUE, n_workers)
      phi_mat <- out$phi_mat
    }
  }
  colnames(phi_mat) <- paste0("phi", 1:n_lengthscales)
  phi_mat <- as.data.frame(phi_mat)
  
  ct_idx <- which(colnames(cell_type_proportions) == focal_cell_type)
  cell_type_prop = as.matrix(cell_type_proportions/apply(cell_type_proportions,1,sum))
  if(sum((cell_type_prop > 1e-7)*(cell_type_prop < 1-(1e-7)))==0){
    stop("Use spacelink_global with cell-type indicators (i.e., cell resolution data).")
  }
  if(sum(cell_type_prop[,ct_idx] > 0.5) / sum(cell_type_prop[,ct_idx] > 0.05) > c1){
    coloc <- c()
  }else{
    coloc <- which(colSums(cell_type_prop[cell_type_prop[,ct_idx] > 0.05,]-cell_type_prop[cell_type_prop[,ct_idx] > 0.05,ct_idx] > 0.05)/colSums(cell_type_prop > 0.05) > c2)
  }
  
  out <- bplapply(1:nrow(Y), function(gene_idx) {
    runtime <- system.time({
      y <- Y[gene_idx,]
      phi_seq <- as.numeric(phi_mat[gene_idx,])
      
      phi_list <- NULL
      for(i in 1:ncol(cell_type_prop)){
        ct_D <- as.matrix(dist(spatial_coords[cell_type_prop[,i]>0.05,]))
        phi_lower <- quantile(1/ct_D[upper.tri(ct_D)], 0.1)
        if(i == ct_idx){
          phi_list[[i]] <- phi_seq[phi_seq > phi_lower]
        }else{
          phi_list[[i]] <- median(phi_seq[phi_seq > phi_lower])
        }
      }
      suppressWarnings(rm(ct_D, phi_lower))
      
      if(length(phi_list[[ct_idx]])==0){
        if(calculate_ESV){
          return(list(time=0, pval_vec=1, ESV=0))
        }else{
          return(list(time=0, pval_vec=1))
        }
      }
      
      # REML Estimation
      D <- as.matrix(dist(spatial_coords))
      Sigma.list <- NULL
      rest_ct_idx <- 1:ncol(cell_type_prop)
      rest_ct_idx <- rest_ct_idx[-ct_idx]
      temp_idx <- 1
      for(i in rest_ct_idx){
        if(i %in% coloc){
          Pi.sq <- (cell_type_prop[,i]%*%t(cell_type_prop[,i]))
          for(phi in phi_list[[i]]){
            Sigma.list[[temp_idx]] <- exp(-phi*D)*Pi.sq
            temp_idx <- temp_idx + 1
          }
        }
      }
      if(is.null(Sigma.list)){
        if(is.null(X)){
          X <- cell_type_prop
        }else{
          X <- cbind(X, cell_type_prop)
        }
        mu_hat <- X %*% solve(t(X) %*% X, t(X) %*% y)
        null_res <- mean((y - mu_hat)^2)
      }else{
        reml_model <- gaston::lmm.aireml(Y = y, X = cbind(rep(1,nrow(cell_type_prop)),X,cell_type_prop), 
                                         K = Sigma.list, min_s2 = 0, verbose = FALSE)
        null_res <- c(reml_model$sigma2, reml_model$tau)
        mu_hat <- cbind(rep(1,nrow(cell_type_prop)),cell_type_prop) %*% reml_model$BLUP_beta
      }
      rm(Sigma.list)
      
      test_res <- score_test_cell_type(ct_idx, null_res, y-mu_hat, spatial_coords, cell_type_prop, phi_list = phi_list, c1 = c1, c2 = c2)
      
      if(calculate_ESV){
        N <- nrow(spatial_coords)
        weight_vec <- sapply(phi_list[[ct_idx]], function(x)sqrt(1 - (N/sum(exp(-2*x*D)))))
        
        # REML Estimation
        Sigma.list <- NULL
        temp_idx <- 1
        Pi.sq <- (cell_type_prop[,ct_idx]%*%t(cell_type_prop[,ct_idx]))
        for(phi in phi_list[[ct_idx]]){
          Sigma.list[[temp_idx]] <- exp(-phi*D)*Pi.sq
          temp_idx <- temp_idx + 1
        }
        rest_ct_idx <- 1:ncol(cell_type_prop)
        rest_ct_idx <- rest_ct_idx[-ct_idx]
        for(i in rest_ct_idx){
          if(i %in% coloc){
            Pi.sq <- (cell_type_prop[,i]%*%t(cell_type_prop[,i]))
            for(phi in phi_list[[i]]){
              Sigma.list[[temp_idx]] <- exp(-phi*D)*Pi.sq
              temp_idx <- temp_idx + 1
            }
          }
        }
        reml_model <- gaston::lmm.aireml(Y = y, X = cbind(rep(1,nrow(cell_type_prop)),X,cell_type_prop), 
                                         K = Sigma.list, min_s2 = 0, verbose = FALSE)
        
        numerator <- sum(reml_model$tau[1:length(phi_list[[ct_idx]])]*weight_vec)
        denominator <- sum(reml_model$tau[1:length(phi_list[[ct_idx]])]) + reml_model$sigma2
        
        if(denominator < 1e-5){
          ESV <- 0
        }else{
          ESV <- numerator/denominator
        }
      }
    })
    if(calculate_ESV){
      return(list(time=runtime[["elapsed"]], pval_vec=test_res$pval_vec, ESV=ESV))
    }else{
      return(list(time=runtime[["elapsed"]], pval_vec=test_res$pval_vec))
    }
  }, BPPARAM = MulticoreParam(workers = n_workers))
  
  results <- data.frame(time=do.call("rbind", lapply(out, function(x)x$time)))
  
  pval_mat <- do.call("rbind", lapply(out, function(x)x$pval_vec))
  suppressWarnings({results$pval <- apply(pval_mat,1,ACAT)})
  
  if(calculate_ESV){
    ESV <- do.call("rbind", lapply(out, function(x)x$ESV))
    results$ESV <- ESV
  }

  if(!is.null(rownames(Y))){
    rownames(results) <- make.unique(rownames(Y))
  }

  return(results)
}
