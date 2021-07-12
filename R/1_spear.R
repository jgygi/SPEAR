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
spear <- function(X = NULL, Y = NULL, Z = NULL, Xobs = NULL, Yobs = NULL, weights_case = NULL, silence = NULL){
  # Pull parameters from SPEARobject:
  if(is.null(X)){X = self$data$train$X}
  if(is.null(Y)){Y = self$data$train$Y}
  if(is.null(Z)){Z = self$data$train$Z}
  if(is.null(Xobs)){Xobs = self$data$train$Xobs}
  if(is.null(Yobs)){Yobs = self$data$train$Yobs}
  if(is.null(weights_case)){weights_case = self$params$weights.case}
  family = self$params$family_encoded
  nclasses = self$params$nclasses
  num_factors = self$params$num_factors
  functional_path = self$params$functional_path
  weights = self$params$weights
  inits_type = self$params$inits_type
  inits_post_mu = self$params$inits_post_mu
  sparsity_upper = self$params$sparsity_upper
  warm_up = self$params$warm_up
  max_iter = self$params$max_iter
  thres_elbo = self$params$thres_elbo
  thres_count = self$params$thres_count
  thres_factor = self$params$thres_factor
  if(is.null(silence)){silence = self$options$quiet}
  if(silence){print_out = 0}else{print_out = self$params$print_out}
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


cv.spear <- function(){
  
  fold_ids = sort(unique(self$params$fold.ids))
  num_factors = self$params$num_factors
  num.cores = self$params$num.cores
  if(num.cores > (self$params$num.folds + 1)){
    num.cores = self$params$num.folds + 1
  }
  
  # add 0 fold_id:
  fold_ids = c(0, fold_ids)
  
  if(is.null(self$params$inits_post_mu)){
    inits_post_mu = matrix(0, nrow = ncol(self$data$train$X), ncol = num_factors)
  }
  if(self$params$inits_type == "None"){
    for(k in 1:num_factors){
      inits_post_mu[,k] = rnorm(ncol(self$data$train$X))
      inits_post_mu[,k] = inits_post_mu[,k]/sqrt(sum(inits_post_mu[,k]^2))
    }
  }else if (self$params$inits_type == "pca"){
    x_svd = svd(self$data$train$Z)
    for(k in 1:num_factors){
      inits_post_mu[,k] = x_svd$v[,k]
    }
  }
  
  run_parallel <- function(fold_id){
    if(fold_id == 0){
      X = self$data$train$X
      Y = self$data$train$Y
      Z = self$data$train$Z
      Xobs = self$data$train$Xobs
      Yobs = self$data$train$Yobs
      res = private$spear(X = X, Y = Y, Z = Z, Xobs = Xobs, Yobs = Yobs)
    }else{
      subsets = which(self$params$fold.ids != fold_id)
      Ycv = self$data$train$Y;
      Xcv = self$data$train$X;
      Zcv = self$data$train$Z;
      Yobs_cv = self$data$train$Yobs;
      Xobs_cv = self$data$train$Xobs;
      if(self$params$only.cross.Y){
        for(j in 1:py){
          Ycv[self$params$fold.ids==fold_id,j] = 0
          Yobs_cv[self$params$fold.ids==fold_id,j] = 0
          weights_case = self$params$weights.case
        }
      }else{
        Xcv = Xcv[subsets,,drop = F]
        Ycv = Ycv[subsets,,drop = F]
        Zcv = Zcv[subsets,,drop = F]
        Yobs_cv = Yobs_cv[subsets,,drop = F]
        Xobs_cv = Xobs_cv[subsets,,drop = F]
        weights_case = self$params$weights.case[subsets]
      }
      fit <- try(private$spear(X = Xcv, Y = Ycv, Z = Zcv, Xobs = Xobs_cv, Yobs = Yobs_cv, weights_case = weights_case, silence = TRUE))
      if(class(fit)=="try-error"){
        stop(paste0("fold",fold_id,":C++failure."))
      }
      res = list(post_betas = fit$post_betas, post_bys = fit$post_bys, 
                 fold_id = fold_id)
    }
    return(res)
  }
  if(self$options$quiet){
    cl <- parallel::makeCluster(num.cores)
  } else {
    cl <- parallel::makeCluster(num.cores, outfile = "")
  }
  parallel::clusterExport(cl, "spear_")
  a <- system.time(
    #results <- parallel::mclapply(fold_ids, run_parallel, mc.cores = numCores)
    results <- parallel::parLapply(cl, fold_ids, fun = run_parallel)
  )
  on.exit(parallel::stopCluster(cl))
  
  if(!self$options$quiet){cat("\n--- All runs finished in ", as.numeric(round(a['elapsed'], 2)), " seconds\n")}
  
  # Make two new return lists:
  factors_coefs = array(0, dim = c(ncol(self$data$train$X), num_factors, max(self$params$fold.ids), nrow(self$params$weights)));
  projection_coefs = array(0, dim = c(num_factors, ncol(self$data$train$Y), max(self$params$fold.ids), nrow(self$params$weights)));
  for(k in 1:(length(results)-1)){
    factors_coefs[,,k,] =results[[k+1]]$post_betas
    projection_coefs[,,k,] = results[[k+1]]$post_bys
  }
  #fit = results[[1]]
  #fit[['post_betas_cv']] = factors_coefs
  #fit[['post_bys_cv']] = projection_coefs
  #return(fit)
  
  return(list(results = results[[1]],
              factors_coefs = factors_coefs,
              projection_coefs = projection_coefs))
}



#' Run SPEAR
#' @examples
#' SPEARobj <- make.SPEARobject(...)
#' 
#' SPEARobj$run.spear()
#' 
#'@export
run.spear <- function(){
  
  tmp <- private$spear()
  
  if(!is.null(self$fit)){
    cat("\n", private$color.text("NOTE:", "yellow"), " this SPEARobject has been trained before. Overwriting previous $fit...\n", sep = "")
  }
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
  
  if(!self$options$quiet){
    cat("\n", private$color.text("NOTE: ", "yellow"), "Setting current.weight.idx to 1\nUse SPEARobject$set.weights(...) to choose a different model.\n", sep = "")
  }
  self$options$current.weight.idx = 1
  if(!self$options$quiet){cat("\n")}
  return(invisible(self))
}



#' Run CV SPEAR
#' @examples
#' SPEARobj <- make.SPEARobject(...)
#' 
#' SPEARobj$run.cv.spear()
#' 
#'@export
run.cv.spear <- function(fold.ids = NULL, num.folds = NULL, only.cross.Y = FALSE, num.cores = NULL){
  # Reset from previous runs...
  self$params$fold.ids <- NULL
  self$params$num.folds <- NULL
  self$params$only.cross.Y = only.cross.Y
  # num.cores:
  if(is.null(num.cores)){self$params$num.cores = parallel::detectCores()}else{self$params$num.cores = num.cores}
  
  # If fold.ids are NULL, generate them...
  if(is.null(fold.ids)){
    # Set seed:
    set.seed(self$params$seed)
    if(is.null(num.folds)){
      self$params$num.folds = 5
    } else {
      self$params$num.folds = num.folds
    }
    fold.ids = sample(1:self$params$num.folds, nrow(self$data$train$X), replace = TRUE)
    if(self$params$family != "gaussian"){
      flag = TRUE
      counter = 0
      while(flag){
        flag = private$check.fold.ids(fold.ids)
        # make new ones if there is an issue
        if(flag){
          fold.ids = sample(1:self$params$num.folds, nrow(self$data$train$X), replace = TRUE)
          counter = counter + 1
          if(counter > 100){
            stop("ERROR: Tried 100 times to generate fold.ids where each class is represented at least twice in the training folds. Try increasing num.folds or removing classes with very few samples.")
          }
        }
      }
    }
  } else {
    # TODO: Check that fold.ids work!
    if(length(fold.ids) != nrow(self$data$train$Y)){
      stop("ERROR: fold.ids provided have a different length (", length(fold.ids), ") than number of samples in $data$train (", nrow(self$data$train$Y), "). Need to match.")
    }
    if(is.null(num.folds)){
      self$params$num.folds = max(fold.ids)
    } else {
      self$params$num.folds = num.folds
    }
    if(any(!fold.ids %in% 1:self$params$num.folds)){
      stop("ERROR: fold.ids provided are outside of the range of possible fold.ids: 1 - ", self$params$num.folds)
    }
    if(length(unique(fold.ids)) != self$params$num.folds){
      stop("ERROR: not enough fold.ids provided for each possible fold (1 - ", self$params$num.folds, "). Need to provide at least one instance for each fold.")
    }
    if(private$check.fold.ids(fold.ids)){
      stop("ERROR: fold.ids provided do not have at least 2 of each possible response class in the training folds. Consider increasing the num.folds argument or removing classes with too few instances.")
    }
  }
  self$params$fold.ids <- fold.ids
  
  # Run cv.spear:
  if(!self$options$quiet){
    cat("Beginning cross-validated SPEAR with K = ", self$params$num.folds, " folds\n", sep = "")
    cat(private$color.text("NOTE: ", "yellow"), " Only printing out progress of one fold for interpretability\n", sep = "")
  }
  tmp <- private$cv.spear()
  
  if(!is.null(self$fit)){
    cat("\n", private$color.text("NOTE: ", "yellow"), " this SPEARobject has already been trained before. Overwriting previous $fit...\n", sep = "")
  }

  self$fit <- list(regression.coefs = tmp$results$post_betas, 
                   projection.coefs.x = tmp$results$post_bxs, 
                   projection.coefs.y = tmp$results$post_bys,
                   nonzero.probs = tmp$results$post_pis, 
                   projection.probs = tmp$results$post_selections, 
                   marginal.probs = tmp$results$post_selections_marginal,
                   joint.probs = tmp$results$post_selections_joint,
                   intercepts.x = tmp$results$interceptsX, 
                   intercepts.y = tmp$results$interceptsY,
                   regression.coefs.cv = tmp$factors_coefs,
                   projection.coefs.y.cv = tmp$projection_coefs,
                   type = "cv")
   
  private$update.dimnames()
  
  
  if(!self$options$quiet){
    cat("\n", private$color.text("NOTE: ", "yellow"), "Setting current.weight.idx to 1\nUse SPEARobject$set.weights(...) to choose a different model.\n", sep = "")
  }
  self$options$current.weight.idx = 1
  if(!self$options$quiet){cat("\n")}
  return(invisible(self))
}



#' Evaluate cv SPEAR object
#' @param nlambda Number of lambdas (defaults to 100)
#' @param calculate.factor.contributions Calculate factor contributions? When `$params$family == "multinomial"` or `"ordinal"` can save time to put `FALSE`. Defaults to `TRUE`.
#' @param max_iter Maximum number of iterations (defaults to 10000)
#' @param multinomial_loss Type of loss for when `$params$family == "multinomial"`. Can be `"deviance"` (default) or `"misclassification"`
#'@export
cv.evaluate <- function(nlambda = 100, calculate.factor.contributions = TRUE, max_iter = 1e4, multinomial_loss = "deviance"){
    if(self$fit$type != "cv"){stop("ERROR: $evaluate.cv must be used after $run.cv.spear. Proper $fit not found.")}
    if(!self$options$quiet){cat("Beginning evaluation of cv.spear...\n")}
    time.start = Sys.time()
    fitted.obj = self$fit
    cv.fact_coefs = fitted.obj$regression.coefs.cv;
    cv.projection_coefs = fitted.obj$projection.coefs.y.cv;
    fact_coefs = fitted.obj$regression.coefs
    projection_coefs = fitted.obj$projection.coefs.y
    X = self$data$train$X
    Y = self$data$train$Y
    Z = self$data$train$Z
    foldid = self$params$fold.ids
    family = self$params$family_encoded
    nclasses = self$params$nclasses
    n = nrow(Y);
    px = ncol(X);
    py = ncol(Y);
    pz = ncol(Z);
    standardize_family = c(2,3)
    nfolds = length(unique(foldid))
    #estimate the across validation error for each of the weight
    num_factors = self$params$num_factors
    num_weights = nrow(self$params$weights)
    intercepts = list()
    scale.factors = array(1, dim  = c(num_factors,nfolds,num_weights))
    for(j in 1:py){
      intercepts[[j]] = matrix(NA, nrow = num_weights,ncol = nclasses[j]-1)
    }
    if(family!=3){
      cvm = matrix(NA, nrow = num_weights, ncol = py)
      cvsd =matrix(NA, nrow = num_weights, ncol = py)
    }else{
      cvm = matrix(NA, nrow = num_weights, ncol = 1)
      cvsd =matrix(NA, nrow = num_weights, ncol = 1)
    }
    
    #rescale the overall coefficients
    Yhat = array(NA, dim = c(n, py, num_weights))
    cmin = 0
    factors.keep = array(NA, dim = c(n, num_factors, py, num_weights))
    Uhat.cv = array(0, dim = c(n,num_factors,num_weights))
    Yhat.cv = array(0, dim = c(n, py, nfolds, num_weights))
    for(l in 1:num_weights){
      for(k in 1:nfolds){
        ucv = array(0, dim = c(n, num_factors))
        beta = cv.fact_coefs[,,k,l]
        ucv = Z[foldid==k,,drop = F]%*%beta
        Uhat.cv[foldid==k,,l] =  ucv
        for(j in 1:py){
          factors.keep[foldid==k,,j,l] = t(apply(ucv, 1, function(z) z*cv.projection_coefs[,j,k,l]))
        }
      }
    }
    for(l in 1:num_weights){
      U0 = array(0, dim = c(n, num_factors))
      U0 = Z%*% fact_coefs[,,l]
      #find the best scaling factors for each weight and response
      if(family != 3){
        r2norm = rep(1, py)
        if(py == 1){
          Yhat[,1,l] =  (U0 %*% projection_coefs[,1,l])
          if(family %in%  standardize_family){
            r2norm =  sqrt(mean(Yhat[,1,l]^2))
            Yhat[,1,l] = Yhat[,1,l]/r2norm
          }
        }else{
          Yhat[,,l] =  (U0 %*% projection_coefs[,,l])
          if(family %in%  standardize_family){
            r2norm =  sqrt(apply(Yhat[,,l]^2,2,function(z) mean(z^2)))
            Yhat[,,l] = apply(Yhat[,,l],2,function(z) z/sqrt(mean(z^2)))
          }
        }
        r2norm_cv = matrix(1, nrow = nfolds, ncol = py)
        for(k in 1:nfolds){
          beta = cv.fact_coefs[,,k,l]
          ucv = Z%*%beta
          b = cv.projection_coefs[,,k,l]
          if(py == 1){
            Yhat.cv[,1,k,l] =  (ucv %*% b)
            if(family %in%  standardize_family){
              r2norm_cv[k,1] = sqrt(mean(Yhat.cv[,1,k,l]^2))
              Yhat.cv[,1,k,l] = Yhat.cv[,1,k,l]/sqrt(mean(Yhat.cv[,1,k,l]^2))
            }
          }else{
            Yhat.cv[,,k,l] =  (ucv %*% b)
            if(family %in%  standardize_family){
              r2norm_cv[k,j] = sqrt(apply(Yhat.cv[,,k,l]^2,2,function(z) mean(z^2)))
              Yhat.cv[,,k,l] =apply(Yhat.cv[,,k,l],2,function(z) z/sqrt(mean(z^2)))
            }
          }
        }
        for(j in 1:py){
          ##note that scaling is only required for Gaussian and logistic!
          y = Y[,j]
          yhat = Yhat[,j,l]
          if(family == 0){
            tmp = lm(y~Yhat[,j,l])
            cmax = tmp$coefficients[2]
          }else if(family==1){
            tmp = glmnet::glmnet(cbind(yhat,rep(0,length(yhat))),y, family = "binomial",lambda = 1e-4)
            cmax = tmp$beta[1,ncol(tmp$beta)]
          }else if(family == 2){
            reverseY = max(y) - y
            reverseY = as.factor(reverseY)
            tmp = ordinalNet(cbind(yhat,rep(0,length(yhat))), reverseY,lambdaVals=1e-4)
            tmp = tmp$coefs
            cmax = tmp[length(tmp)-1]
          }
          chats = seq(0, cmax, length.out = nlambda)
          errs = array(NA, dim = c(n,length(chats)))
          for(k in 1:nfolds){
            yhat_cv = Yhat.cv[,j,k,l]
            if(family == 1){
              tmp0 = glmnet::glmnet(cbind(yhat_cv[foldid!=k],rep(0,sum(foldid!=k))),y[foldid!=k], family = "binomial",lambda = 1e-4)
              a0 = tmp0$a0
              c0 = tmp0$beta[1]
            }else if(family == 2){
              tmp0 = ordinalNet(cbind(yhat_cv[foldid!=k],rep(0,sum(foldid!=k))),reverseY[foldid!=k],lambdaVals=1e-4)
              a0 = tmp0$coefs[1:(nclasses[j]-1)]
              c0 = tmp0$coefs[nclasses[j]]
            }
            for(ll in 1:length(chats)){
              if(family == 0){
                a = mean(y[foldid!=k] - chats[ll] * yhat_cv[foldid!=k])
                errs[foldid==k,ll] = (y[foldid == k] - a - chats[ll] * yhat_cv[foldid==k])^2
              }else if(family == 1){
                if(ll == 1){
                  tmp <-glm(y[foldid!=k] ~ offset(-yhat_cv[foldid!=k]*chats[ll]), family = "binomial")
                  a = tmp$coefficients[1]
                }else{
                  if(abs(c0) >=abs(chats[ll])){
                    tmp <-glm(y[foldid!=k] ~ offset(-yhat_cv[foldid!=k]*chats[ll]), family = "binomial",  start = a)
                    a = tmp$coefficients[1]
                  }else{
                    a = a0
                  }
                }
                Probs = 1/(1+exp(-chats[ll] * yhat_cv[foldid == k]-a))
                errs[foldid==k, ll] = -2*(y[foldid == k] * log(Probs+1e-10) +
                                            (1 - y[foldid == k]) * log(1 - Probs+1e-10));
              }else if(family == 2){
                if(ll == 1){
                  tmp <-polr(reverseY[foldid!=k] ~ offset(-chats[ll]*yhat_cv[foldid !=k]))
                  a = tmp$zeta
                }else{
                  if(abs(c0) >=abs(chats[ll])){
                    tmp <-polr(reverseY[foldid!=k] ~ offset(-chats[ll]*yhat_cv[foldid !=k]),start = a)
                    a = tmp$zeta
                  }else{
                    a = a0
                  }
                }
                a1 = a[(nclasses[j]-1):1]
                ##deviance loss
                Pmat0 = matrix(0, ncol = max(y), nrow = sum(foldid == k))
                Pmat = matrix(0, ncol = max(y)+1, nrow = sum(foldid == k))
                y_extended = matrix(0, ncol = max(y)+1, nrow = sum(foldid == k))
                for(kk in 1:length(a)){
                  Pmat0[,kk] = 1/(1+exp(-chats[ll] * yhat_cv[foldid == k]-a1[kk]))
                }
                for(kk in 1:ncol(Pmat)){
                  y_extended[y[foldid==k]==(kk-1),kk]=1
                  if(kk==1){
                    Pmat[,kk] = 1 - Pmat0[,kk]
                  }else if(kk==ncol(Pmat)){
                    Pmat[,kk] = Pmat0[,kk-1]
                  }else{
                    Pmat[,kk] = Pmat0[,kk-1] - Pmat0[,kk]
                  }
                }
                errs[foldid==k,ll] = -2 * apply(log(Pmat+1e-10) * y_extended,1,sum)
              }
            }
          }
          cv_tmp = apply(errs,2,mean)
          cvsd_tmp = apply(errs,2,sd)/sqrt(n-1)
          cvm[l,j] = min(cv_tmp)
          cvsd[l,j] = cvsd_tmp[which.min(cv_tmp)]

          projection_coefs[,j,l] = projection_coefs[,j,l] *chats[which.min(cv_tmp)]
          if(family %in% standardize_family){
            projection_coefs[,j,l] = projection_coefs[,j,l]/r2norm[j]
          }
          if(family == 0){
            intercepts[[j]][l,] = mean(y - mean(yhat *chats[which.min(cv_tmp)] ))
          }else if(family==1){
            tmp = glm(y ~ offset(chats[which.min(cv_tmp)] * yhat), family = "binomial")
            a = tmp$coefficients[1]
            intercepts[[j]][l,] = a
          }else{
            tmp = polr(reverseY ~ offset(-chats[which.min(cv_tmp)] * yhat))
            a = tmp$zeta
            a = a[length(a):1]
            intercepts[[j]][l,] = a
          }
          ###scale the projection coefficients post_by in the original data object
          for(foldind in 1:nfolds){
            cv.projection_coefs[,j,foldind,l] =cv.projection_coefs[,j,foldind,l]/r2norm_cv[foldind,j]
            cv.projection_coefs[,j,foldind,l]= cv.projection_coefs[,j,foldind,l] *chats[which.min(cv_tmp)]
          }
        }
      }else{
        if(family %in% standardize_family){
          r2norm = sqrt(apply(U0,2,function(z) mean(z^2)))
          U0std  = U0
        }else{
          U0std = U0
          r2norm = rep(1, ncol(U0))
        }
        if(ncol(Y) <=2){
          stop('multinomial class < 3, please use binomial.')
        }
        ycollapsed = rep(0, nrow(Y))
        for(j in 1:ncol(Y)){
          ycollapsed[Y[,j]==1] = (j-1)
        }
        #get penalty lists
        fitted_multinomial = glmnet::glmnet(x =U0std, y = ycollapsed, standardize = F, family = 'multinomial', nlambda = nlambda)
        lambdas =fitted_multinomial$lambda
        # perform cv fit with given penalty lists
        probs_predictions = array(NA, dim  = c(nrow(Y),ncol(Y),length(lambdas)))
        tmp.fit = list()
        for(fold_id in 1:nfolds){
          test_id = which(foldid == fold_id)
          train_id = which(foldid != fold_id)
          U0cv = array(0, dim = c(n, num_factors))
          #' This is where the error occurs previously!
          #' changed from U0cv * scale.factors[,fold_id,l] to
          #'  U0cv %*%diag(scale.factors[,fold_id,l])
          U0cv = Z%*% cv.fact_coefs[,,fold_id, l]
          if(family %in% standardize_family){
            tmp = sqrt(apply(U0cv,2,function(z) mean(z^2)))
            scale.factors[,fold_id,l] = r2norm/tmp
            U0cv.std  = U0cv %*%diag(scale.factors[,fold_id,l])
          }else{
            U0cv.std  = U0cv
          }
          tmp.fit[[fold_id]] = glmnet::glmnet(x =U0cv.std[train_id,], y = ycollapsed[train_id], standardize = F, family = 'multinomial', lambda = lambdas)
          preds = predict(tmp.fit[[fold_id]], U0cv.std)
          total = apply(preds,c(1,3),function(z) matrixStats::logSumExp(z))
          for(kkk in 1:dim(preds)[2]){
            preds[,kkk,] =  preds[,kkk,] -total
          }
          preds = exp(preds)
          probs_predictions[test_id,,1:dim(preds)[3]] = preds[test_id,,]
        }
        if(multinomial_loss=="deviance"){
          cvs = apply(probs_predictions,c(3), function(z) -2*apply(as.matrix(log(z+1e-8) * Y),1,sum))
        }else{
          #classification loss
          cvs = apply(probs_predictions,c(3), function(z){
            tmp = apply(z,1,which.max)
            tmp1 = apply(Y==1, 1, which)
            ifelse(tmp==tmp1, 0, 1)
          } )
        }
        # choose penalty and model refit
        cvs0 = apply(cvs,2,mean)
        idx = which.min(cvs0)
        cvm[l] =cvs0[idx]
        cvsd[l] = sd(cvs[,idx])/sqrt(nrow(cvs))
        for(j in 1:py){
          projection_coefs[,j,l] = fitted_multinomial$beta[[j]][,idx]
          intercepts[[j]][l] = fitted_multinomial$a0[j,idx]
          for(fold_id in 1:nfolds){
            cv.projection_coefs[,j,fold_id,l] = tmp.fit[[fold_id]]$beta[[j]][,idx]/scale.factors[,fold_id,l]
          }
        }
      }
    }
    ##regression coefficients
    reg_coefs = array(0, dim = c(pz,py,num_weights))
    for(l in 1:num_weights){
      for(j in 1:py){
        reg_coefs[,j,l] = fact_coefs[,,l]%*%projection_coefs[,j,l]
      }
    }
    #replace deviance contribution to spearman correlation
    factor_contributions = array(NA,dim = c(num_factors, py, num_weights))
    factor_contributions_pvals = array(NA,dim = c(num_factors, py, num_weights))
    if(calculate.factor.contributions){
      for(l in 1:num_weights){
        if(!self$options$quiet)
          cat(paste0("--- Calculating contribution for ", private$color.text(paste0("weight.x = ", self$params$weights[l,1], " | weight.y = ",self$params$weights[l,2]), "green"), " - (", l, "/", num_weights, ")", " ---\n"))

        for(k in 1:num_factors){
          for(j in 1:py){
            yhat = Uhat.cv[,k,l]
            y = Y[,j]
            for(fold_id in 1:nfolds){
              b = cv.projection_coefs[k,j,fold_id,l]
              yhat[foldid==fold_id] = yhat[foldid==fold_id] *b
            }
            yhat = yhat + rnorm(n = length(yhat), mean = 0, sd = sqrt(var(y))*1/n)
            tmp_pearson = cor(y, yhat)
            suppressWarnings( tmp_spearman <- cor.test(y, yhat, method = 'spearman') )
            if( tmp_pearson<0){
              factor_contributions[k,j,l] = 0
              factor_contributions_pvals[k,j,l] = 1
            }else{
              factor_contributions[k,j,l] = (tmp_spearman$estimate)^2
              factor_contributions_pvals[k,j,l] = tmp_spearman$p.value
            }
          }
        }
      }
    }
    
    if(!self$options$quiet){cat("\n--- Finished evaluation in ", round(as.numeric(Sys.time() - time.start), 4), " seconds\n")}
    
    self$fit$cv.eval = list(projection_coefs_scaled = projection_coefs,
                cv.projection_coefs_scaled = cv.projection_coefs,
                reg_coefs = reg_coefs,
                intercepts = intercepts,
                cvm = cvm, cvsd = cvsd,
                factor_contributions = factor_contributions,
                factor_contributions_pvals = factor_contributions_pvals)
    
    self$fit$projection.coefs.y.scaled = projection_coefs
    self$fit$projection.coefs.y.cv.scaled = cv.projection_coefs
    self$fit$intercepts.scaled = intercepts
    self$fit$factor.contributions = factor_contributions
    self$fit$factor.contributions.pvals = factor_contributions_pvals
    self$fit$loss <- list(cvm = cvm, cvsd = cvsd)
    
    if(!self$options$quiet){cat("\n")}
    return(invisible(self))
}



