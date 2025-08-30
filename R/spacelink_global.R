#################################################
#' @param Y Normalized count matrix. Rows and columns indicate genes and spots, respectively.
#' @param spatial_coords Spatial coordinates matrix. Rows indicate spots.
#' @param X Spatial features expected to include intercept. If null, X = rep(1,N) where N is number of spots.
#' @param lengthscales If NULL, length-scales are calculated by using our two-step approach.
#' @param n_lengthscales Number of length-scales (i.e., number of kernels)
#' @param M It controls the minimum length-scale. Minimum length-scale is set to be the minimum distance multiplied by M.
#' @param n_workers Number of workers for parallelization.
#################################################


spacelink_global <- function(Y, spatial_coords, X = NULL, lengthscales = NULL, n_lengthscales = 5, M = 1, n_workers = 1){
  # Create matrix of phi (= 1/lengthscale)
  if(!is.null(lengthscales)){
    if(is.matrix(lengthscales)){
      phi_mat <- 1/lengthscales
    }else{
      phi_mat <- t(matrix(1/lengthscales, length(lengthscales), nrow(Y)))
    }
    n_lengthscales <- ncol(phi_mat)
  }else{
    out <- select_lengthscales(Y, spatial_coords, n_lengthscales, M = M, is.inverse = TRUE, n_workers)
    phi_mat <- out$phi_mat
    time0 <- out$time
  }
  colnames(phi_mat) <- paste0("phi", 1:n_lengthscales)
  phi_mat <- as.data.frame(phi_mat)

  if(is.null(X)){
    Y <- Y - rowMeans(Y)
  }else{
    Y <- Y - t(X %*% solve(crossprod(X), crossprod(X,t(Y))))
  }

  out <- bplapply(1:nrow(Y), function(gene_idx) {
    runtime <- system.time({
      y <- as.numeric(Y[gene_idx,])
      phi_seq <- phi_mat[gene_idx,]
      nnls_res <- run_nnls(y, spatial_coords, phi_seq = phi_seq)
      nnls_res <- cbind(nnls_res, phi_seq)

      sigma.sq_vec <- nnls_res[,grep("sigma", colnames(nnls_res))]
      tau.sq <- nnls_res[,grep("tau",colnames(nnls_res))]
      nnls_res$ESV <- calculate_ESV(sigma.sq_vec, tau.sq, spatial_coords, phi_seq = phi_seq)

      test_res <- score_test(y, spatial_coords, phi_seq = phi_seq)
    })

    nnls_res$time <- runtime[["elapsed"]]

    return(list(nnls_res=nnls_res, pval_vec=test_res$pval_vec))
  }, BPPARAM = MulticoreParam(workers = n_workers))

  results <- do.call("rbind", lapply(out, function(x)x$nnls_res))
  if(is.null(lengthscales)){results$time <- results$time + time0}

  pval_mat <- do.call("rbind", lapply(out, function(x)x$pval_vec))
  colnames(pval_mat) <- paste0("pval",1:ncol(pval_mat))
  results <- cbind(results, pval_mat)
  suppressWarnings({results$pval <- apply(pval_mat,1,ACAT)})
  results$padj <- p.adjust(results$pval, method = "BH")

  results$ESV_adj <- results$ESV
  results$ESV_adj[results$padj > 0.05] <- 0

  if(!is.null(rownames(Y))){
    rownames(results) <- make.unique(rownames(Y))
  }

  return(results)
}
