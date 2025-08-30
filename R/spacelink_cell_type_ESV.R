####################################
#' @param global_test_res Spacelink global SVG test results.
#' @param ct_idx Index of cell-type to be tested.
#' @param Y Normalized count matrix. Rows and columns indicate genes and spots, respectively.
#' @param spatial_coords Spatial coordinates matrix. Rows indicate spots.
#' @param cell_type_data Cell-type data matrix. Cell-type proportions (e.g., RCTD) or cell-type abundance (e.g., Cell2location).
#' @param c1 Gating parameter 1.
#' @param c2 Gating parameter 2.
#' @param is_phi_reduce If TRUE, kernel bandwidths are restricted to be smaller than the maximum distance among cells within target cell type.
#' @param n_workers Number of workers for parallelization.
####################################

spacelink_cell_type_ESV <- function(global_test_res, ct_idx, Y, spatial_coords, cell_type_data, c1 = 0.25, c2 = 0.2,
                                     is_phi_reduce = TRUE, n_workers = 1){
  
  cell_type_prop = cell_type_data/apply(cell_type_data,1,sum)
  
  if(sum((cell_type_prop > 1e-7)*(cell_type_prop < 1-(1e-7)))==0){
    stop("Use spacelink_global with cell-type indicators (i.e., cell resolution data).")
  }
  
  phi_mat <- global_test_res[,grep("phi",colnames(global_test_res))]
  rm(global_test_res)
  
  out <- bplapply(1:nrow(Y), function(gene_idx) {
    runtime <- system.time({
      y <- Y[gene_idx,]
      phi_seq <- as.numeric(phi_mat[gene_idx,])
      
      phi_list <- NULL
      for(i in 1:ncol(cell_type_data)){
        if(is_phi_reduce){
          ct_D <- as.matrix(dist(spatial_coords[cell_type_prop[,i]>0.05,]))
          phi_lower <- quantile(1/ct_D[upper.tri(ct_D)], 0.1)
          if(i == ct_idx){
            phi_list[[i]] <- phi_seq[phi_seq > phi_lower]
          }else{
            phi_list[[i]] <- median(phi_seq[phi_seq > phi_lower])
          }
        }else{
          if(i == ct_idx){
            phi_list[[i]] <- phi_seq
          }else{
            phi_list[[i]] <- median(phi_seq)
          }
        }
      }
      suppressWarnings(rm(ct_D, phi_lower))
      
      ########################
      if(length(phi_list[[ct_idx]])==0){
        return(list(time=0, ESV=0))
      }
      ########################
      
      if(sum(cell_type_prop[,ct_idx] > 0.5) / sum(cell_type_prop[,ct_idx] > 0.05) > c1){
        temp_ct_idx <- c()
      }else{
        temp_ct_idx <- which(colSums(cell_type_prop[cell_type_prop[,ct_idx] > 0.05,]
                                     -cell_type_prop[cell_type_prop[,ct_idx] > 0.05,ct_idx] > 0.05)
                             /colSums(cell_type_prop > 0.05) > c2)
      }
      
      N <- nrow(spatial_coords)
      D <- as.matrix(dist(spatial_coords))
      
      weight_vec <- sapply(phi_list[[ct_idx]], function(x)sqrt(1 - (N/sum(exp(-2*x*D)))))
      
      # REML Estimation
      Sigma.list <- NULL
      temp_idx <- 1
      Pi.sq <- (cell_type_data[,ct_idx]%*%t(cell_type_data[,ct_idx]))
      for(phi in phi_list[[ct_idx]]){
        Sigma.list[[temp_idx]] <- exp(-phi*D)*Pi.sq
        temp_idx <- temp_idx + 1
      }
      rest_ct_idx <- 1:ncol(cell_type_data)
      rest_ct_idx <- rest_ct_idx[-ct_idx]
      for(i in rest_ct_idx){
        if(i %in% temp_ct_idx){
          Pi.sq <- (cell_type_data[,i]%*%t(cell_type_data[,i]))
          for(phi in phi_list[[i]]){
            Sigma.list[[temp_idx]] <- exp(-phi*D)*Pi.sq
            temp_idx <- temp_idx + 1
          }
        }
      }
      reml_model <- gaston::lmm.aireml(Y = y, X = cbind(rep(1,nrow(cell_type_data)),cell_type_data), 
                                       K = Sigma.list, min_s2 = 0, verbose = FALSE)
      
      numerator <- sum(reml_model$tau[1:length(phi_list[[ct_idx]])]*weight_vec)
      denominator <- sum(reml_model$tau[1:length(phi_list[[ct_idx]])]) + reml_model$sigma2

      if(denominator < 1e-5){
        ESV <- 0
      }else{
        ESV <- numerator/denominator
      }
    })
    return(list(time=runtime[["elapsed"]], ESV=ESV))
  }, BPPARAM = MulticoreParam(workers = n_workers))
  
  results <- data.frame(time=do.call("rbind", lapply(out, function(x)x$time)))
  
  ESV <- do.call("rbind", lapply(out, function(x)x$ESV))
  results$ESV <- ESV
  
  if(!is.null(rownames(Y))){
    rownames(results) <- make.unique(rownames(Y))
  }
  
  return(results)
}