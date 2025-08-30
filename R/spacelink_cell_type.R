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


spacelink_cell_type <- function(global_test_res, ct_idx, Y, spatial_coords, cell_type_data, c1 = 0.25, c2 = 0.2,
                                is_phi_reduce = TRUE, n_workers = 1){

  cell_type_prop = cell_type_data/apply(cell_type_data,1,sum)

  if(sum((cell_type_prop > 1e-7)*(cell_type_prop < 1-(1e-7)))==0){
    stop("Use spacelink_global with cell-type indicators (i.e., cell resolution data).")
  }

  phi_mat <- global_test_res[,grep("phi",colnames(global_test_res))]
  rm(global_test_res)

  out <- bplapply(1:nrow(Y), function(gene_idx) {
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

    if(length(phi_list[[ct_idx]])==0){
      return(list(time=0, pval_vec=1))
    }

    runtime <- system.time({
      if(sum(cell_type_prop[,ct_idx] > 0.5) / sum(cell_type_prop[,ct_idx] > 0.05) > c1){
        temp_ct_idx <- c()
      }else{
        temp_ct_idx <- which(colSums(cell_type_prop[cell_type_prop[,ct_idx] > 0.05,]
                                     -cell_type_prop[cell_type_prop[,ct_idx] > 0.05,ct_idx] > 0.05)
                             /colSums(cell_type_prop > 0.05) > c2)
      }

      # REML Estimation
      D <- as.matrix(dist(spatial_coords))
      Sigma.list <- NULL
      rest_ct_idx <- 1:ncol(cell_type_data)
      rest_ct_idx <- rest_ct_idx[-ct_idx]
      temp_idx <- 1
      for(i in rest_ct_idx){
        if(i %in% temp_ct_idx){
          Pi.sq <- (cell_type_data[,i]%*%t(cell_type_data[,i]))
          for(phi in phi_list[[i]]){
            Sigma.list[[temp_idx]] <- exp(-phi*D)*Pi.sq
            temp_idx <- temp_idx + 1
          }
        }
      }
      if(is.null(Sigma.list)){
        X <- cell_type_data
        mu_hat <- X %*% solve(t(X) %*% X, t(X) %*% y)
        null_res <- mean((y - mu_hat)^2)
      }else{
        reml_model <- gaston::lmm.aireml(Y = y, X = cbind(rep(1,nrow(cell_type_data)),cell_type_data),
                                         K = Sigma.list, min_s2 = 0, verbose = FALSE)
        null_res <- c(reml_model$sigma2, reml_model$tau)
        mu_hat <- cbind(rep(1,nrow(cell_type_data)),cell_type_data) %*% reml_model$BLUP_beta
      }
      rm(D, Sigma.list)

      test_res <- score_test_cell_type(ct_idx, null_res, y-mu_hat, spatial_coords, cell_type_data, phi_list = phi_list, c1 = c1, c2 = c2)
    })

    return(list(time=runtime[["elapsed"]], pval_vec=test_res$pval_vec))
  }, BPPARAM = MulticoreParam(workers = n_workers))

  results <- data.frame(time=do.call("rbind", lapply(out, function(x)x$time)))

  pval_mat <- do.call("rbind", lapply(out, function(x)x$pval_vec))
  suppressWarnings({results$pval <- apply(pval_mat,1,ACAT)})

  if(!is.null(rownames(Y))){
    rownames(results) <- make.unique(rownames(Y))
  }

  return(results)
}
