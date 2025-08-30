score_test_cell_type <- function(ct_idx, null_res, y, spatial_coords, cell_type_data, phi_list, c1, c2){
  N <- nrow(spatial_coords)
  D <- as.matrix(dist(spatial_coords))
  
  cell_type_prop = cell_type_data/apply(cell_type_data,1,sum)
  
  if(sum(cell_type_prop[,ct_idx] > 0.5) / sum(cell_type_prop[,ct_idx] > 0.05) > c1){
    temp_ct_idx <- c()
  }else{
    temp_ct_idx <- which(colSums(cell_type_prop[cell_type_prop[,ct_idx] > 0.05,]
                                 -cell_type_prop[cell_type_prop[,ct_idx] > 0.05,ct_idx] > 0.05)
                         /colSums(cell_type_prop > 0.05) > c2)
  }

  null_res <- as.numeric(null_res)
  null_cov <- null_res[1]*diag(N)
  rest_ct_idx <- 1:ncol(cell_type_data)
  rest_ct_idx <- rest_ct_idx[-ct_idx]
  temp_idx <- 1
  for(i in rest_ct_idx){
    if(i %in% temp_ct_idx){
      Pi.sq <- cell_type_data[,i]%*%t(cell_type_data[,i])
      for(phi in phi_list[[i]]){
        temp_idx = temp_idx + 1
        if(null_res[temp_idx] > 0){
          null_cov <- null_cov + null_res[temp_idx]*exp(-phi*D)*Pi.sq
        }
      }
    }
  }
  if(sum(diag(null_cov) == 0)>0){diag(null_cov)[diag(null_cov) == 0] <- 1e-10}
  null_cov_inv <- invM_arma(null_cov)
  rm(null_cov)
  
  n_kernels <- length(phi_list[[ct_idx]])
  
  tr.Sigma.seq <- rep(NA, n_kernels)
  tr.Sigma.sq.seq <- rep(NA, n_kernels)
  ytSy.seq <- rep(NA, n_kernels)
  Pi.sq <- cell_type_data[,ct_idx]%*%t(cell_type_data[,ct_idx])
  for(i in 1:n_kernels){
    Sigma.prod <- (exp(-phi_list[[ct_idx]][i]*D)*Pi.sq) %*% null_cov_inv
    ytSy.seq[i] <- vMMv_arma(y, null_cov_inv, Sigma.prod, y)
    tr.Sigma.seq[i] <- sum(diag(Sigma.prod))
    tr.Sigma.sq.seq[i] <- sum(t(Sigma.prod)*Sigma.prod)
  }
  rm(Sigma.prod, Pi.sq, D, null_cov_inv)
  
  pval_vec <- rep(NA, n_kernels)
  for(i in 1:n_kernels){
    score <- (tr.Sigma.seq[i]/tr.Sigma.sq.seq[i])*ytSy.seq[i]
    df_val = (tr.Sigma.seq[i]^2)/tr.Sigma.sq.seq[i]
    pval_vec[i] <- exp(pchisq(score, df_val, lower.tail=FALSE, log.p=TRUE))
  }
  return(list(pval_vec=pval_vec))
}