#' SPEAR - SuPervised Bayes fActor for Multi-omics
#' imports:
#'@useDynLib SPEAR, .registration=TRUE
#'@importFrom ordinalNet ordinalNet
#'@importFrom MASS polr
#'@importFrom glmnet glmnet
#'@importFrom reshape2 melt
#'@importFrom Rcpp evalCpp
#'@import R6
#'@import parallel
#'@import ggplot2
#'@import cowplot
#'@import dplyr
#'@import stringr
#'@export
spear <- function(){
  # Pull parameters from SPEARobject:
  X = self$data$train$X
  Y = self$data$train$Y
  Z = self$data$train$Z
  Xobs = self$data$train$Xobs
  Yobs = self$data$train$Yobs
  family = self$params$family_encoded
  nclasses = self$params$nclasses
  num_factors = self$params$num_factors
  functional_path = self$params$functional_path
  weights_case = self$params$weights.case
  weights = self$params$weights
  inits_type = self$params$inits_type
  inits_post_mu = self$params$inits_post_mu
  sparsity_upper = self$params$sparsity_upper
  warm_up = self$params$warm_up
  max_iter = self$params$max_iter
  thres_elbo = self$params$thres_elbo
  thres_count = self$params$thres_count
  thres_factor = self$params$thres_factor
  print_out = self$params$print_out
  seed = self$params$seed
  a0 = self$inits$a0 
  b0 = self$inits$b0 
  a1 = self$inits$a1 
  b1 = self$inits$b1
  a2 = self$inits$a2 
  b2 = self$inits$b2 
  L1 = self$inits$L1 
  L2 = self$inits$L2
  robust_eps = self$inits$robust_eps
  
  # Set seed:
  set.seed(seed)
  
  # Set up weights:
  if(all(weights[,2] == 1)){
    type_weights = "xonly"
  }else if(all(weights[,1] == 1)){
    type_weights = "yonly"
  }else{
    type_weights = "both"
  }
  all_ws = weights
  if(type_weights != "yonly"){
    one_penalty_idx = which(weights[,1] == 1)[1]
  }else{
    one_penalty_idx = 1
  }
  
  # Dimensions:
  px = ncol(X); py = ncol(Y); pz = ncol(Z); n = nrow(Y)
  interceptsY = list()
  interceptsX = rep(0, px)
  for(j in 1:py){
    interceptsY[[j]] = rep(0, nclasses[j]-1)
  }
  
  # Initialization:
  post_mu = array(0, dim = c(ncol(Z), num_factors))
  post_sigma2 = array(0.1, dim=c(ncol(Z), num_factors));
  post_pi = array(1, dim=c(ncol(Z), num_factors));
  
  # Check for inits_post_mu (if provided):
  if(!is.null(inits_post_mu)){
    if((ncol(inits_post_mu)!=num_factors) | (nrow(inits_post_mu)!= ncol(Z))){
      stop("wrong initialization dimension for post_mu!")
    }
    post_mu = inits_post_mu;
    # Else, check for other types of initialization:
  }else if(inits_type == "None"){
    post_mu = array(rnorm(pz*num_factors), dim = c(pz, num_factors))
    for(k in 1:num_factors){
      post_mu[,k] =post_mu[,k]/sqrt(sum(post_mu[,k])^2)
    }
  }else if(inits_type == "pca"){
    z_svd = svd(Z)
    for(k in 1:num_factors){
      post_mu[,k] = z_svd$v[,k]
    }
  }else if (inits_type == "sparsepca"){
    z_svd = spca(Z, num_factors,  sparse="varnum", type = "predictor",
                 para = rep(min(ceiling(sqrt(nrow(X))), ncol(X)/2),num_factors))
    for(k in 1:num_factors){
      post_mu[,k] = z_svd$v[,k]
    }
  }
  
  # Initialize other variables:
  post_tmuX =array(0, dim=c(px, num_factors));
  post_tsigma2X = array(1e-4, dim=c(px, num_factors));
  post_tpiX = array(1.0, dim=c(px, num_factors));
  post_tpiX_marginal = array(1.0, dim=c(px, num_factors));
  post_tmuY =array(0, dim=c(py, num_factors));
  post_tsigma2Y = array(1e-4, dim=c(py, num_factors));
  post_tpiY = array(1.0, dim=c(py, num_factors));
  tauY = array(1, dim=c(py,num_factors));
  tauZ = array(1, dim=c(pz,num_factors));
  tauZ[-c(1:px),] = 1
  post_tmuY[,1] = 1;
  log_pi =array(log(.5), dim=c(pz, num_factors));
  log_minus_pi = array(log(.5), dim=c(pz, num_factors));
  nuYmat = array(2, dim = c(n, py))
  nuXmat = array(2, dim = c(n, px))
  meanFactors = array(0, dim=c(n, num_factors));
  post_a0 = matrix(1, ncol = num_factors, nrow = length(functional_path))
  post_a1 = matrix(1, ncol = num_factors, nrow = length(functional_path))
  post_b0 = matrix(1, ncol = num_factors, nrow = length(functional_path))
  post_b1 = matrix(1, ncol = num_factors, nrow = length(functional_path))
  post_a2x = rep(1, ncol(X))
  post_b2x = rep(1, ncol(X))
  post_a2y = rep(1, ncol(Y))
  post_b2y = rep(1, ncol(Y))
  
  # Record the model estimated with weights all 1 for initial start of y
  one_post_mu = post_mu; one_post_pi = post_pi;one_post_sigma2 = post_sigma2;
  one_post_tmuX =array(0, dim=c(px, num_factors));
  one_post_tsigma2X = array(1e-4, dim=c(px, num_factors));
  one_post_tpiX = array(1.0, dim=c(px, num_factors));
  one_post_tpiX_marginal = array(1.0, dim=c(px, num_factors));
  one_post_tmuY =array(0, dim=c(py, num_factors));
  one_post_tsigma2Y = array(1e-4, dim=c(py, num_factors));
  one_post_tpiY = array(1.0, dim=c(py, num_factors));
  one_tauY = array(1, dim=c(py,num_factors));
  one_tauZ = array(1, dim=c(pz,num_factors));
  one_tauZ[-c(1:px),] = 1
  one_post_tmuY[,1] = 1;
  one_log_pi =array(log(.5), dim=c(pz, num_factors));
  one_log_minus_pi = array(log(.5), dim=c(pz, num_factors));
  one_nuYmat = array(2, dim = c(n, py))
  one_nuXmat = array(2, dim = c(n, px))
  one_meanFactors = array(0, dim=c(n, num_factors));
  one_post_a0 = matrix(1, ncol = num_factors, nrow = length(functional_path))
  one_post_a1 = matrix(1, ncol = num_factors, nrow = length(functional_path))
  one_post_b0 = matrix(1, ncol = num_factors, nrow = length(functional_path))
  one_post_b1 = matrix(1, ncol = num_factors, nrow = length(functional_path))
  one_post_a2x = rep(1, ncol(X))
  one_post_b2x = rep(1, ncol(X))
  one_post_a2y = rep(1, ncol(Y))
  one_post_b2y = rep(1, ncol(Y))
  lowers = rep(0, length(all_ws))
  post_betas = array(NA, dim  = c(ncol(X),num_factors , nrow(all_ws)))
  post_bys = array(NA, dim = c(num_factors, ncol(Y), nrow(all_ws)))
  post_bxs = array(NA, dim = c(ncol(X), num_factors,  nrow(all_ws)))
  post_pis = array(NA, dim = c(ncol(X), num_factors, nrow(all_ws)))
  post_selections = array(NA, dim = c(ncol(X), num_factors, nrow(all_ws)))
  post_selections_marginal = array(NA, dim = c(ncol(X), num_factors,  nrow(all_ws)))
  
  # Go through all weights, run SPEAR each time:
  for(idx_w in 1:nrow(all_ws)){
    weights = rep(all_ws[idx_w,1], ncol(X))
    weights_y = rep(all_ws[idx_w,2], ncol(Y))
    lowers[idx_w] = max(1-all_ws[idx_w,1], 0)
    set.seed(seed)
    if(idx_w == 1){
      warm_up1 = warm_up
    }else{
      warm_up1 = 1
    }
    if((type_weights != "yonly") & (idx_w == (nrow(all_ws)+1))){
      post_mu = one_post_mu
      post_sigma2 =  one_post_sigma2
      post_pi =  one_post_pi
      post_tmuX =one_post_tmuX
      post_tsigma2X = one_post_tsigma2X
      post_tpiX = one_post_tpiX
      post_tpiX_marginal = one_post_tpiX_marginal
      post_tmuY =one_post_tmuY
      post_tsigma2Y = one_post_tsigma2Y
      post_tpiY = one_post_tpiY
      tauY = one_tauY
      tauZ = one_tauZ
      tauZ =one_tauZ
      post_tmuY= one_post_tmuY
      log_pi =one_log_pi
      log_minus_pi = one_log_minus_pi
      nuYmat =  one_nuYmat
      nuXmat = one_nuXmat 
      meanFactors = one_meanFactors
      post_a0 =one_post_a0 
      post_a1 =one_post_a1
      post_b0 =one_post_b0
      post_b1 =one_post_b1
      post_a2x = one_post_a2x 
      post_b2x = one_post_b2x
      post_a2y = one_post_a2y 
      post_b2y = one_post_b2y
      meanFactors = one_meanFactors
    }
    
    if(print_out > 0)
      cat(paste0("\n--- ", private$color.text(paste0("Running weight.x = ", all_ws[idx_w,1], " | weight.y = ",all_ws[idx_w,2]), "green"), " ---\n"))
    
    spear_(family  = family, Y = Y, X = X, Yobs = Yobs, Xobs = Xobs, Z = Z,
           nclasses =  nclasses,  functional_path = functional_path,
           weights = weights,  weights0 = weights_y, 
           weights_case = weights_case,
           num_factors = num_factors, warm_up = warm_up1,
           max_iter = max_iter, thres_elbo = thres_elbo,  thres_count = thres_count,
           thres_factor = thres_factor,  a0  = a0, b0 = b0,
           a1 = a1, b1 = b1,a2 = a2,b2 = b2, lower = lowers[idx_w], print_out = print_out,
           interceptsX = interceptsX, interceptsY = interceptsY,
           post_mu = post_mu, post_sigma2 = post_sigma2,post_pi = post_pi, 
           post_tmuX = post_tmuX, post_tsigma2X = post_tsigma2X, post_tpiX = post_tpiX, post_tpiX_marginal = post_tpiX_marginal,
           post_tmuY = post_tmuY, post_tsigma2Y = post_tsigma2Y, post_tpiY = post_tpiY,
           tauY = tauY, tauZ = tauZ,log_pi = log_pi,log_minus_pi = log_minus_pi, 
           nuXmat = nuXmat, nuYmat = nuYmat,
           post_a0 = post_a0, post_b0 = post_b0,
           post_a1 = post_a1, post_b1 = post_b1,
           post_a2x = post_a2x, post_b2x = post_b2x,
           post_a2y = post_a2y, post_b2y = post_b2y,
           meanFactors = meanFactors, 
           seed0 = seed,robust_eps =robust_eps, alpha0 = sparsity_upper, L = L1,L2 = L2)
    
    if(idx_w==one_penalty_idx){
      one_post_mu = post_mu
      one_post_sigma2 = post_sigma2
      one_post_pi = post_pi
      one_post_tmuX =post_tmuX
      one_post_pi = post_pi
      one_post_tsigma2X = post_tsigma2X
      one_post_tpiX = post_tpiX
      one_post_tpiX_marginal = post_tpiX_marginal
      one_post_tmuY = post_tmuY
      one_post_tsigma2Y = post_tsigma2Y
      one_post_tpiY =post_tpiY
      one_tauY = tauY
      one_tauZ = tauZ
      one_log_pi =log_pi
      one_log_minus_pi = log_minus_pi
      one_nuYmat =  nuYmat
      one_nuXmat = nuXmat 
      meanFactors = one_meanFactors
      one_post_a0 = post_a0 
      one_post_a1 = post_a1
      one_post_b0 = post_b0
      one_post_b1 = post_b1
      one_post_a2x = post_a2x 
      one_post_b2x = post_b2x
      one_post_a2y = post_a2y 
      one_post_b2y = post_b2y
      one_meanFactors = meanFactors
    }
    ###return both the factors after re-order and sign-fliping
    post_beta =array(0, dim = dim(post_mu))
    post_bx =  post_tmuX *  post_tpiX
    post_beta = post_mu*post_pi
    meanFactors = Z%*%post_beta
    post_by = post_tmuY
    cors = matrix(0, nrow  = num_factors, ncol = ncol(Y))
    if(family != 0){
      ##ordinal regression
      for(j in (1:ncol(Y))){
        y = Y[,j]
        for(k in 1:num_factors){
          data = data.frame(cbind(meanFactors[,k],y))
          colnames(data) = c("x", "y")
          data[,2] = as.factor(data[,2])
          labels = unique(data[,2])
          cors[k,j] = cor(y, meanFactors[,k], method = "spearman")
          cors[k,j] = cors[k,j] * sqrt(sum(meanFactors[,k]^2))
        }
      }
    }else{
      cors = cov(meanFactors, Y)
    }
    cors_abs = sqrt(apply(cors^2,1,sum))
    ordering = order(cors_abs, decreasing = T)
    for(k in 1:num_factors){
      k0 = ordering[k]
      aligning = sum(cors[k0,])
      post_beta[,k] = (post_mu[,k0] * post_pi[,k0])
      post_bx[,k] = post_tmuX[,k0] *post_tpiX[,k0]
      post_by[,k] = post_tmuY[,k0] *post_tpiY[,k0]
      if(aligning < 0){
        post_beta[,k] = - post_beta[,k]
        post_bx[,k] = -post_bx[,k]
        post_by[,k] = -post_by[,k]
      }
    }
    post_mu = post_mu[,ordering]
    post_pi = post_pi[,ordering]
    post_tpiX = post_tpiX[,ordering]
    post_tpiX_marginal = post_tpiX_marginal[,ordering]
    post_betas[,,idx_w] = post_beta
    post_bys[,,idx_w] = post_by
    post_bxs[,,idx_w] = post_bx
    post_pis[,,idx_w] = post_pi
    post_selections[,,idx_w] = post_tpiX
    post_selections_marginal[,,idx_w]  = post_tpiX_marginal
  }
  post_selections_joint = ifelse(post_selections<=post_selections_marginal, post_selections, post_selections_marginal)
  
  return(list(post_betas = post_betas, 
              post_bxs = post_bxs, 
              post_bys = post_bys,
              post_pis = post_pis, 
              post_selections = post_selections, 
              post_selections_marginal = post_selections_marginal,
              post_selections_joint = post_selections_joint,
              interceptsX = interceptsX, 
              interceptsY = interceptsY))
}



#' Run SPEAR
#' @examples
#' SPEARobj <- make_spear_model(...)
#' 
#' SPEARobj$run.spear()
#' 
#'@export
run.spear <- function(){
  
  tmp <- private$spear()
  
  self$fit <- list(regression.coefs = tmp$post_betas, 
                   projection.coefs.x = tmp$post_bxs, 
                   projection.coefs.y = tmp$post_bys,
                   nonzero.probs = tmp$post_pis, 
                   projection.probs = tmp$post_selections, 
                   marginal.probs = tmp$post_selections_marginal,
                   joint.probs = tmp$post_selections_joint,
                   intercepts.x = tmp$interceptsX, 
                   intercepts.y = tmp$interceptsY,
                   type = "regular")
  
  private$update.dimnames()
  
  cat("\n", private$color.text("NOTE: ", "yellow"), "Setting current.weight.idx to 1\nUse SPEARobject$set.weights(...) to choose a different model.\n", sep = "")
  self$options$current.weight.idx = 1
}

