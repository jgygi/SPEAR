#' SuPervised Bayes fActor for Multi-omics
#'@useDynLib SPEAR, .registration=TRUE
#'@importFrom ordinalNet ordinalNet
#'@importFrom MASS polr
#'@importFrom glmnet glmnet
#'@importFrom reshape2 melt
#'@importFrom GGally ggpairs
#'@importFrom Rcpp evalCpp
#'@import parallel
#'@import ggplot2
#'@import cowplot
#'@import dplyr
#'@import stringr
#'@param Y Response matrix (can be multidimensional for gaussian data).
#'@param X Assay matrix.
#'@param Z Complete feature matrix (usually the features are the imputed version of X, other features are attached to the end).
#'@param family 0=gaussian (multiple response); 1 = binary(multiple response); 2 = ordinal (multiple response); 3 = multinomial (single response)
#'@param ws A vector of weights that you want to try out, default is ws = c(0).
#'@param num.factors Number of factors estimated.
#'@param functional_path Grouping structure.
#'@param pattern_samples Sample indexes for each missing pattern, default NULL (one pattern). It should be a partition of all 1-n indexes.
#'@param pattern_features Feature indexes for each missing pattern, default NULL (one pattern). 
#'@param inits_type Initialization type, can be None, pca, sparsepac.
#'@param warm_up: warm up iterations for the inference.
#'@param max_iter: max number of iterations.
#'@param thres_elbo: if EBLO increase by less than thres_elbo, clock +1; otherwise, the clock is reset to 0.
#'@param thres_count: stop if clock has reached thres_count.
#'@param print_out: print out the progress for every print_out iterations.
#'@param a0: hyper parametr, no need to tune usually.
#'@param b0: hyper parametr, no need to tune usually.
#'@param a1: hyper parametr, no need to tune usually.
#'@param b1: hyper parametr, no need to tune usually.
#'@param a2: hyper parametr, no need to tune usually.
#'@param b2: hyper parametr, no need to tune usually.
#'@param seed: random seed number.
#'@param robust_eps: robust_eps
#'@param sparsity_upper: Sparsity parameter for feature selection
#'@param L: parameter
#'@export
spear <- function(X, Xobs, Y, Yobs, Z, family, nclasses, ws, num_factors, 
                  functional_path, case.weights = NULL, ws_y = NULL,
pattern_samples = NULL, pattern_features = NULL,
                  inits_type = "pca", warm_up = 100, max_iter = 1000,
                  thres_elbo = 0.01, thres_count = 5, thres_factor = 1e-8, print_out = 10,
                  a0 = 1e-2, b0 = 1e-2, a1 = sqrt(nrow(X)), b1 = sqrt(nrow(X)),
                  a2= sqrt(nrow(X)), b2 = sqrt(nrow(X)), 
                  inits_post_mu = NULL,seed = 1, robust_eps = 1.0/(nrow(X)), 
                  sparsity_upper = 0.5, L = nrow(X)/log(ncol(X))){
  if(is.null(dim(Y))){
    Y = matrix(Y, ncol = 1)
    Yobs = matrix(Yobs, ncol = 1)
  }
  px = ncol(X); py = ncol(Y); pz = ncol(Z); n = nrow(Y)
  interceptsY = list()
  interceptsX = rep(0, px)
  for(j in 1:py){
      interceptsY[[j]] = rep(0, nclasses[j]-1)
  }
  if(is.null(pattern_samples) | is.null(pattern_features)){
    pattern_samples = list()
    pattern_features = list()
    pattern_samples[[1]] = c(1:n)
    pattern_features[[1]] = c(1:px)
  }else if(length(pattern_samples)!=length(pattern_features)){
    stop("feature patterns and sample patterns do not match!")
  }else{
    tmp1 = pattern_samples[[1]]
    tmp2 = pattern_samples[[1]]
    for(k in 1:length(pattern_samples)){
      tmp1 = intersect(tmp1, pattern_samples[[k]])
      tmp2 = sort(union(tmp2, pattern_samples[[k]]))
      if((length(tmp1) > 0 & k > 1) | length(tmp2)!=n){
        stop("pattern_samples is not a partition of all samples!")
      }
    }
  }
  num_patterns = length(pattern_samples)
  post_mu = array(0, dim = c(ncol(Z), num_factors, length(pattern_samples)))
  post_sigma2 = array(0.1, dim=c(ncol(Z), num_factors, length(pattern_samples)));
  post_pi = array(1, dim=c(ncol(Z), num_factors, length(pattern_samples)));
  if(!is.null(inits_post_mu)){
    if((ncol(inits_post_mu)!=num_factors) | (nrow(inits_post_mu)!= ncol(Z))){
      stop("wrong initialization dimension for post_mu!")
    }
    for(k in 1:num_patterns){
      post_mu[,,k] = inits_post_mu;
    }
  }else if(inits_type == "None"){
    post_mu = array(rnorm(pz*num_factors*num_patterns), dim = c(pz, num_factors, num_factors))
    for(k in 1:num_factors){
      for(j in 1:num_patterns){
        post_mu[,k,j] =post_mu[,k,j]/sqrt(sum(post_mu[,k,j])^2)
      }
    }
  }else if(inits_type == "pca"){
    z_svd = svd(Z)
    for(k in 1:num_factors){
      for(j in 1:num_patterns){
        post_mu[,k,j] = z_svd$v[,k]
      }
    }
  }else if (inits_type == "sparsepca"){
    z_svd = spca(Z, num_factors,  sparse="varnum", type = "predictor",
                 para = rep(min(ceiling(sqrt(nrow(X))), ncol(X)/2),num_factors))
    for(k in 1:num_factors){
      for(j in 1:num_patterns){
        post_mu[,k,j] = z_svd$v[,k]
      }
    }
  }
  post_tmuX =array(0, dim=c(px, num_factors));
  post_tsigma2X = array(1e-4, dim=c(px, num_factors));
  post_tpiX = array(1.0, dim=c(px, num_factors));
  post_tpiX_marginal = array(1.0, dim=c(px, num_factors));
  post_tmuY =array(0, dim=c(py, num_factors));
  post_tsigma2Y = array(1e-4, dim=c(py, num_factors));
  post_tpiY = array(1.0, dim=c(py, num_factors));
  tauY = array(1, dim=c(py,num_factors));
  tauZ = array(0.01, dim=c(pz,num_factors));
  tauZ[-c(1:px),] = 1
  post_tmuY[,1] = 1;
  log_pi =array(log(.5), dim=c(pz, num_factors));
  log_minus_pi = array(log(.5), dim=c(pz, num_factors));
  nuYmat = array(2, dim = c(n, px))
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
  post_betas = array(NA, dim  = c(ncol(X),num_factors , num_patterns, length(ws)))
  post_bys = array(NA, dim = c(num_factors, ncol(Y), length(ws)))
  post_bxs = array(NA, dim = c(ncol(X), num_factors, length(ws)))
  post_pis = array(NA, dim = c(ncol(X), num_factors, num_patterns, length(ws)))
  post_selections = array(NA, dim = c(ncol(X), num_factors, length(ws)))
  post_selections_marginal = array(NA, dim = c(ncol(X), num_factors, length(ws)))
  lowers = rep(0, length(ws))
  for(idx_w in 1:length(ws)){
    weights = rep(ws[idx_w], ncol(X))
    lowers[idx_w] = max(1-ws[idx_w], 0)
    set.seed(seed)
    if(idx_w == 1){
      warm_up1 = warm_up
    }else{
      warm_up1 = 1
    }
    set.seed(seed)
    spear_(family  = family, Y = Y, X = X, Yobs = Yobs, Xobs = Xobs, Z = Z,
           nclasses =  nclasses,  functional_path = functional_path,
           pattern_samples = pattern_samples, pattern_features = pattern_features,
           weights = weights,  num_factors = num_factors, warm_up = warm_up1,
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
           meanFactors = meanFactors, seed0 = seed,robust_eps =robust_eps, alpha0 = sparsity_upper, L = L)
    ###return both the factors after re-order and sign-fliping
    post_beta =array(0, dim = dim(post_mu))
    post_bx =  post_tmuX *  post_tpiX
    for(j in 1:num_patterns){
      ii = pattern_samples[[j]]
      jj = pattern_features[[j]]
      for(k in 1:num_factors){
        post_beta[,k,j] = (post_mu[,k,j] * post_pi[,k,j])
        meanFactors[ii,k] = Z[ii,] %*%post_beta[,k,j]
      }
    }
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
      for(j in 1:num_patterns){
        post_beta[,k,j] = (post_mu[,k0,j] * post_pi[,k0,j])
      }
      post_bx[,k] = post_tmuX[,k0] *post_tpiX[,k0]
      post_by[,k] = post_tmuY[,k0] *post_tpiY[,k0]
      if(aligning < 0){
        for(j in 1:num_patterns){
          post_beta[,k,j] = - post_beta[,k,j]
        }
        post_bx[,k] = -post_bx[,k]
        post_by[,k] = -post_by[,k]
      }
    }
    for(j in 1:num_patterns){
      post_mu[,,j] = post_mu[,ordering,j]
      post_pi[,,j] = post_pi[,ordering,j]
    } 
    post_tpiX = post_tpiX[,ordering]
    post_tpiX_marginal = post_tpiX_marginal[,ordering]
    post_betas[,,,idx_w] = post_beta
    post_bys[,,idx_w] = post_by
    post_bxs[,,idx_w] = post_bx
    post_pis[,,,idx_w] = post_pi
    post_selections[,,idx_w] = post_tpiX
    post_selections_marginal[,,idx_w]  = post_tpiX_marginal
    cat(paste0("*** weights: ",ws[idx_w], "------------------------"))
  }
  post_selections_joint = ifelse(post_selections<=post_selections_marginal, post_selections, post_selections_marginal)
  # hist(post_selections_marginal[,1,idx_w], breaks = 100)
  # hist(post_selections[,1,idx_w], breaks = 100)
  # hist(post_selections_joint[,1,idx_w], breaks = 100)
  return(list(post_betas = post_betas, post_bys = post_bys, post_bxs =post_bxs,
              post_pis = post_pis, post_selections = post_selections, 
              post_selections_marginal = post_selections_marginal,
              post_selections_joint = post_selections_joint,
              interceptsX = interceptsX, interceptsY = interceptsY))
  
}


#' Cross-fit of SPEAR
#'@param Y Response matrix (can be multidimensional for gaussian data).
#'@param X Feature matrix.
#'@param ws A vector of weights that you want to try out, default is ws = c(0).
#'@param num.factors Number of factors estimated.
#'@param functional_path Grouping structure.
#'@param foldid CV foldid
#'@param inits_type Initialization type, can be None, pca, sparsepac.
#'@param warm_up: warm up iterations for the inference.
#'@param max_iter: max number of iterations.
#'@param thres_elbo: if EBLO increase by less than thres_elbo, clock +1; otherwise, the clock is reset to 0.
#'@param thres_count: stop if clock has reached thres_count.
#'@param print_out: print out the progress for every print_out iterations.
#'@param a0: hyper parametr, no need to tune usually.
#'@param b0: hyper parametr, no need to tune usually.
#'@param a1: hyper parametr, no need to tune usually.
#'@param b1: hyper parametr, no need to tune usually.
#'@param a2: hyper parametr, no need to tune usually.
#'@param b2: hyper parametr, no need to tune usually.
#'@param seed: random seed number.
#'@param sparsity_upper: Sparsity parameter
#'@param robust_eps: robust_eps
#'@param run.debug: debug?
#'@param L: parameter
#'@export
cv.spear <- function(X, Xobs, Y, Yobs, Z, family, nclasses, ws, num_factors, 
                     functional_path, foldid = foldid, case.weights = NULL, ws_y = NULL, 
                     pattern_samples = NULL, pattern_features = NULL,
                     inits_type = "pca", warm_up = 100, max_iter = 1000,
                     thres_elbo = 0.01, thres_count = 5, thres_factor = 1e-8, print_out = 10,
                     a0 = 1e-2, b0 = 1e-2, a1 = sqrt(nrow(X)), b1 = sqrt(nrow(X)),
                     a2= sqrt(nrow(X)), b2 = sqrt(nrow(X)), robust_eps =1/nrow(X),
                     sparsity_upper = 0.1, L = nrow(X)/log(nrow(X)),inits_post_mu = NULL,seed = 1, crossYonly = F, numCores = NULL, run.debug = FALSE){
  fold_ids = sort(unique(foldid))
  fold_ids = c(0, fold_ids)
  px = ncol(X); py = ncol(Y); pz = ncol(Z); n = nrow(Y)
  num_patterns = length(pattern_samples)
  if(is.null(pattern_samples) | is.null(pattern_features)){
    pattern_samples = list()
    pattern_features = list()
    pattern_samples[[1]] = c(1:n)
    pattern_features[[1]] = c(1:px)
  }else if(length(pattern_samples)!=length(pattern_features)){
    stop("feature patterns and sample patterns do not match!")
  }else{
    tmp1 = pattern_samples[[1]]
    tmp2 = pattern_samples[[1]]
    for(k in 1:length(pattern_samples)){
      tmp1 = intersect(tmp1, pattern_samples[[k]])
      tmp2 = sort(union(tmp2, pattern_samples[[k]]))
      if((length(tmp1) > 0 & k > 1) | length(tmp2)!=n){
        stop("pattern_samples is not a partition of all samples!")
      }
    }
  }
  if(is.null(inits_post_mu)){
    inits_post_mu = matrix(0, nrow = ncol(X), ncol = num_factors)
  }
  if(inits_type == "None"){
    for(k in 1:num_factors){
      inits_post_mu[,k] = rnorm(ncol(X))
      inits_post_mu[,k] = inits_post_mu[,k]/sqrt(sum(inits_post_mu[,k]^2))
    }
  }else if (inits_type == "pca"){
    x_svd = svd(Z)
    for(k in 1:num_factors){
      inits_post_mu[,k] = x_svd$v[,k]
    }
  }
  run_parallel <- function(fold_id){
    if(fold_id == 0){
      res = spear(family  = family, Y = Y, X = X, Yobs = Yobs, Xobs = Xobs, Z = Z,
                   nclasses =  nclasses,  functional_path = functional_path,
                   pattern_samples = pattern_samples, pattern_features = pattern_features,
                   ws = ws,  num_factors = num_factors, warm_up = warm_up,
                   max_iter = max_iter, thres_elbo = thres_elbo,  thres_count = thres_count,
                   thres_factor = thres_factor,  print_out = print_out, a0  = a0, b0 = b0,
                   a1 = a1, b1 = b1,a2 = a2,b2 = b2, inits_post_mu = inits_post_mu, seed = seed,
                  robust_eps=robust_eps, sparsity_upper = sparsity_upper, L = L)
      
    }else{
      subsets = which(foldid != fold_id)
      Ycv = Y;
      Xcv = X;
      Zcv = Z;
      Yobs_cv = Yobs;
      Xobs_cv = Xobs;
      pattern_samples_cv = pattern_samples;
      if(crossYonly){
        for(j in 1:py){
          Ycv[foldid==fold_id,j] = 0
          Yobs_cv[foldid==fold_id,j] = 0
        }
      }else{
        Xcv = Xcv[subsets,]
        Ycv = Ycv[subsets,]
        Zcv = Zcv[subsets,]
        Yobs_cv = Yobs_cv[subsets,]
        Xobs_cv = Xobs_cv[subsets,]
        ids = 1:length(subsets)
        ids_names = subsets
        names(ids) = subsets
        for(k in 1:num_patterns){
          ll = pattern_samples_cv[[k]]
          ll1 = intersect(ll, subsets)
          pattern_samples_cv[[k]] = as.integer(ids[as.character(ll1)])
        }
      }
      fit <- try(spear(family  = family, Y = Ycv, X = Xcv, Yobs = Yobs_cv, Xobs = Xobs_cv, Z = Zcv,
                  nclasses =  nclasses,  functional_path = functional_path,
                  pattern_samples = pattern_samples_cv, pattern_features = pattern_features,
                  ws = ws,  num_factors = num_factors, warm_up = warm_up,
                  max_iter = max_iter, thres_elbo = thres_elbo,  thres_count = thres_count,
                  thres_factor = thres_factor,  print_out = print_out, a0  = a0, b0 = b0,
                  a1 = a1, b1 = b1,a2 = a2,b2 = b2, inits_post_mu = inits_post_mu, seed = seed,robust_eps=robust_eps,
                  sparsity_upper = sparsity_upper, L = L))
      if(class(fit)=="try-error"){
        stop(paste0("fold",fold_id,":C++failure."))
      }
      res = list(post_betas = fit$post_betas, post_bys = fit$post_bys, 
                 fold_id = fold_id)
    }
    return(res)
  }
  if(is.null(numCores)){
    numCores <- detectCores()
  }
  a <- system.time(
   results <- mclapply(fold_ids, run_parallel, mc.cores = numCores)
   #results <- sapply(fold_ids, run_parallel)
  )
  # a <- system.time(
  #   results <- run_parallel(0)
  # )
  print(a)
  if(run.debug){
    print(results)
  }
  factors_coefs = array(0, dim = c(ncol(X), num_factors, num_patterns,max(foldid), length(ws)));
  projection_coefs = array(0, dim = c(num_factors, ncol(Y), max(foldid), length(ws)));
  for(k in 1:(length(results)-1)){
    factors_coefs[,,,k,] =results[[k+1]]$post_betas
    projection_coefs[,,k,] = results[[k+1]]$post_bys
  }
  return(list(results = results[[1]],
              factors_coefs = factors_coefs,
              projection_coefs = projection_coefs, foldid = foldid))
  
}


