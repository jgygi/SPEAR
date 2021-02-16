dataGen_gassian <- function(N = 500, Ntest = 2000, P = 500, D = 4, seed = 123, num_factors = 5, c = 1, 
                            pi = 0.2, eta = 1, num_specific =D-2, Ymodel = "factor",
                            pi_reg = 0.05){
  set.seed(seed)
  Theta0 = list(); Gamma0 = list()
  X = list(); Xte = list()
  for(d in 1:D){
    Theta0[[d]] = matrix(rnorm(P*num_factors, sd = 1), ncol = P) * c
    Gamma0[[d]] = matrix(rbinom(P*num_factors, size = 1, prob  = pi), ncol = P)
  }
  for(k in 1:num_factors){
    ii = sample(1:D, num_specific)
    for(d in 1:D){
      if(!(d%in%ii)){
        Gamma0[[d]][k,] = 0
        Theta0[[d]] = Theta0[[d]] * Gamma0[[d]]
      }
    }
  }
  U0 = matrix(rnorm(N*num_factors), nrow = N)
  U0te = matrix(rnorm(Ntest*num_factors), nrow = Ntest)
  X = list(); Xte = list()
  scale_mean.X = c(); scale_sd.X = c()
  for(d in 1:D){
    X[[d]] = scale(U0 %*% Theta0[[d]]+matrix(rnorm(N*P), ncol = P))
    Xte[[d]] = scale(U0te %*% Theta0[[d]]+matrix(rnorm(Ntest*P), ncol = P))
    tmp1 = apply(X[[d]], 2, mean);
    tmp2 =  apply(X[[d]], 2, sd);
    scale_mean.X =c(scale_mean.X, tmp1)
    scale_sd.X = c(scale_sd.X, tmp2)
    X[[d]] = t(apply(X[[d]], 1, function(z) (z - tmp1)/tmp2))
    Xte[[d]] = t(apply(Xte[[d]], 1, function(z) (z - tmp1)/tmp2))
  }
  if(Ymodel == "factor"){
    Y = (U0[,1] + U0[,2])*sqrt(eta/2) + rnorm(N)
    Yte = (U0te[,1] + U0te[,2])* sqrt(eta/2)+rnorm(Ntest)
    Ytruth = (U0[,1] + U0[,2])*sqrt(eta/2)
    Ytruth_te =  (U0te[,1] + U0te[,2])* sqrt(eta/2)
  }else{
    Xcombine = X[[1]]
    Xcombine.te = Xte[[1]]
    for(d in 2:D){
      Xcombine = cbind(Xcombine, X[[d]])
      Xcombine.te = cbind(Xcombine.te, Xte[[d]])
    }
    beta =rnorm(ncol(Xcombine)) * rbinom(ncol(Xcombine), size = 1, prob = pi_reg) 
    beta = beta/sqrt(sum(beta^2)) * sqrt(eta)
    Ytruth =  Xcombine%*%beta
    Ytruth_te =  Xcombine.te%*% beta
    Y = Ytruth + rnorm(N)
    Yte = Ytruth_te+rnorm(Ntest)
  }
  scale_mean.y = mean(Y); scale_sd.y = sd(Y);
  y = (Y - scale_mean.y)/scale_sd.y; yte = (Yte - scale_mean.y)/scale_sd.y;
  ytruth = (Ytruth - scale_mean.y)/scale_sd.y; ytruth.te = (Ytruth_te - scale_mean.y)/scale_sd.y;
  data.tr = preparation(Y = y, X = X, family = 0, path.type = "assay")
  data.te = preparation(Y = yte, X = Xte, family = 0,path.type = "assay")
  data.tr$truth = ytruth
  data.tr$xlist = X
  data.tr$U = U0
  data.te$truth = ytruth.te
  data.te$U = U0te
  data.te$xlist = Xte
  return(list(data.tr = data.tr, data.te = data.te))
}


dataGen_ordinal <- function(N = 500, Ntest = 2000, P = 500, levels = 7,
                             D = 4, seed = 123, num_factors = 5, c = 1, 
                             pi = 0.2, eta = 1, num_specific =D-2, Ymodel = "factor",
                             pi_reg = 0.05){
  set.seed(seed)
  Theta0 = list(); Gamma0 = list()
  X = list(); Xte = list()
  for(d in 1:D){
    Theta0[[d]] = matrix(rnorm(P*num_factors, sd = 1), ncol = P) * c
    Gamma0[[d]] = matrix(rbinom(P*num_factors, size = 1, prob  = pi), ncol = P)
  }
  for(k in 1:num_factors){
    ii = sample(1:D, num_specific)
    for(d in 1:D){
      if(!(d%in%ii)){
        Gamma0[[d]][k,] = 0
        Theta0[[d]] = Theta0[[d]] * Gamma0[[d]]
      }
    }
  }
  U0 = matrix(rnorm(N*num_factors), nrow = N)
  U0te = matrix(rnorm(Ntest*num_factors), nrow = Ntest)
  X = list(); Xte = list()
  scale_mean.X = c(); scale_sd.X = c()
  for(d in 1:D){
    X[[d]] = scale(U0 %*% Theta0[[d]]+matrix(rnorm(N*P), ncol = P))
    Xte[[d]] = scale(U0te %*% Theta0[[d]]+matrix(rnorm(Ntest*P), ncol = P))
    tmp1 = apply(X[[d]], 2, mean);
    tmp2 =  apply(X[[d]], 2, sd);
    scale_mean.X =c(scale_mean.X, tmp1)
    scale_sd.X = c(scale_sd.X, tmp2)
    X[[d]] = t(apply(X[[d]], 1, function(z) (z - tmp1)/tmp2))
    Xte[[d]] = t(apply(Xte[[d]], 1, function(z) (z - tmp1)/tmp2))
  }
  if(Ymodel == "factor"){
    Ytruth = (U0[,1] + U0[,2])*sqrt(eta/2)
    Ytruth_te =  (U0te[,1] + U0te[,2])* sqrt(eta/2)
  }else{
    Xcombine = X[[1]]
    Xcombine.te = Xte[[1]]
    for(d in 2:D){
      Xcombine = cbind(Xcombine, X[[d]])
      Xcombine.te = cbind(Xcombine.te, Xte[[d]])
    }
    beta =rnorm(ncol(Xcombine)) * rbinom(ncol(Xcombine), size = 1, prob = pi_reg) 
    beta = beta/sqrt(sum(beta^2)) * sqrt(eta)
    Ytruth =  Xcombine%*%beta
    Ytruth_te =  Xcombine.te%*% beta
  }
  ###
  tmp = c(Ytruth, Ytruth_te)
  tmp1 = c((2*levels-1):1)/(2*levels)
  tmp1 = c(min(tmp1), sort(sample(tmp1[2:(length(tmp1)-1)], levels-3)),max(tmp1))
  tmp1 = sort(tmp1, decreasing = T)
  intercepts = quantile(tmp, tmp1)
  Pmat0 = matrix(0, ncol = levels-1, nrow = N)
  Pmat0.te = matrix(0, ncol = levels-1, nrow = Ntest)
  for(l in 1:(levels-1)){
    Pmat0[,l] = 1/(1+exp(-Ytruth-intercepts[l]))
    Pmat0.te[,l] = 1/(1+exp(-Ytruth_te-intercepts[l]))
  }
  Pmat = matrix(0, ncol = levels, nrow = N)
  Pmat.te = matrix(0, ncol = levels, nrow = Ntest)
  for(l in 1:levels){
    if(l == 1){
      Pmat[,l] = 1 - Pmat0[,l]
      Pmat.te[,l] = 1-Pmat0.te[,l]
    }else if(l == levels){
      Pmat[,l] = Pmat0[,l-1]
      Pmat.te[,l] = Pmat0.te[,l-1]
    }else{
      Pmat[,l] = Pmat0[,l-1] -  Pmat0[,l]
      Pmat.te[,l] = Pmat0.te[,l-1] -Pmat0.te[,l]
    }
  }
  Y = apply(Pmat, 1,function(z) which(rmultinom(n=1, size = 1, prob = z)==1))-1
  Yte =  apply(Pmat.te, 1,function(z) which(rmultinom(n=1, size = 1, prob = z)==1))-1
  data.tr = preparation(Y = Y, X = X, family = 2, path.type = "assay")
  data.te = preparation(Y = Yte, X = Xte, family = 2,path.type = "assay")
  data.tr$truth = Ytruth
  data.tr$xlist = X
  data.tr$U = U0
  data.te$truth = Ytruth_te
  data.te$U = U0te
  data.te$xlist = Xte
  return(list(data.tr = data.tr, data.te = data.te))
}

prob_calculation <- function(yhat, intercept){
  Pmat0 = matrix(0, ncol = length(intercept), nrow = length(yhat))
  Pmat = matrix(0, ncol = length(intercept)+1, nrow = length(yhat))
  for(k in 1:ncol(Pmat0)){
    Pmat0[,k] = 1.0/(1+exp(-yhat - intercept[k]))
  }
  for(k in 1:ncol(Pmat)){
    if(k == 1){
      Pmat[,k] = 1-Pmat0[,k]
    }else if(k == ncol(Pmat)){
      Pmat[,k] = Pmat0[,k-1]
    }else{
      Pmat[,k] = Pmat0[,k-1] - Pmat0[,k]
    }
  }
  return(Pmat)
}


loss_calculation <- function(Phat, y, type = "deviance"){
  S = array(0, dim  = dim(Phat))
  for(k in 1:ncol(Phat)){
    S[y == (k-1),k] = 1
  }
  if(type == "deviance"){
    loss = -2 *sum(S * log(Phat) + (1-S) * log(1-Phat))
  }else{
    loss = mean((apply(Phat,1,which.max) - 1) != y)
  }
  loss
}


#' Create the data from lists, and create the functional pathway list.
#'@param Y Response array.
#'@param List
#'@param refit Whether to refit the model with elasticNet. If FALSE, the coefficients are from the spear model.
#'@param rule Rule to use for picking cross validation model: "min" or "1se".
#'@param standardize Whether to standardize the data.
#'@param alpha alpha = 1 corresponds to lasso and alpha = 0 corresponds to ridge.
#'@export
preparation <- function(Y,  X, family, pattern_samples = NULL, pattern_assays = NULL,
                        path.type = "assay", 
                        other.path = NULL){
  ###check input types
  if(is.null(dim(Y))){
    Y = matrix(Y, ncol = 1)
  }
  py = ncol(Y)
  nclasses = rep(2, ncol(Y))
  for(j in 1:py){
    if(family != 0){
      labels = sort(unique(Y[,j]))
      labels.correct = 0:(length(labels)-1)
      if(sum(labels!= labels.correct) != 0){
        stop("class labels are not consecutive integers starting from 0")
      }
      nclasses[j] = length(labels)
    }else{
    }
  }
  if(!(path.type %in% c("assay", "none", "other"))){
    stop("Unsupported grouping structure.")
  }else if(path.type == "other" & is.null(other.path)){
    stop("other.grouping is missing.")
  }
  ##prepare the data
  X_ = X[[1]]
  n = nrow(X[[1]])
  p = rep(0, length(X))
  for(d in 1:length(X)){
    p[d] = ncol(X[[d]])
    if(d > 1){
      X_ = cbind(X_, X[[d]])
    }
  }
  functional_path = list()
  if(path.type == "assay"){
    for(i in 1:length(p)){
      if(i == 1){
        functional_path[[i]]  = 1:p[i]
      }else{
        functional_path[[i]] = (sum(p[1:(i-1)])+1):sum(p[1:(i)])
      }
      functional_path[[i]]  = functional_path[[i]] -1
    }
  }else if(group.type == "other"){
    functional_path= other.grouping;
  }else{
    functional_path[[1]] = c(1:sum(p));
  }
  px = ncol(X_)
  if(is.null(pattern_samples) | is.null(pattern_assays)){
    pattern_samples = list()
    pattern_assays = list()
    pattern_samples[[1]] = c(1:n)
    pattern_assays[[1]] = c(1:length(p))
  }else if(length(pattern_samples)!=length(pattern_assays)){
    stop("feature patterns and sample patterns do not match!")
  }else{
    tmp1 = pattern_samples[[1]]
    tmp2 = pattern_samples[[1]]
    for(k in 1:length(pattern_samples)){
      tmp1 = intersect(tmp1, pattern_samples[[k]])
      tmp2 = sort(union(tmp2, pattern_samples[[k]]))
      if(length(tmp1) > 0 | length(tmp2)!=n){
        stop("pattern_samples is not a partition of all samples!")
      }
    }
  }
  pattern_features = list()
  psum = cumsum(p)
  psum = c(0, psum)
  for(k in 1:length(pattern_assays)){
    pattern_features[[k]] = c(NA)
    for(l in pattern_assays[[k]]){
      pattern_features[[k]] = c(pattern_features[[k]], (psum[l]+1):psum[l+1])
    }
    pattern_features[[k]] = pattern_features[[k]][-1]
  }
  return(list(Y = Y, X = X_, functional_path = functional_path,
              pattern_samples = pattern_samples, pattern_features = pattern_features,
              nclasses = nclasses
              ))
}


#cross-validation: select models with high accuracy for predicting Y.
#' We consider two types the model:
#' (1) yhat = (U * b) * chat (refit = F)
#' (2) yhat = U * bhat (refit = T)
#'@param fitted.obj Fitted object from cv.spear
#'@param X Feature array.
#'@param Y Response array.
#'@param rule Rule to use for picking cross validation model: "min" or "1se".
#'@param standardize Whether to standardize the data.
#'@param alpha alpha = 1 corresponds to lasso and alpha = 0 corresponds to ridge.
#'@export
cv.evaluation <- function(fitted.obj, X, Y, Z, family, nclasses, 
                          pattern_samples, pattern_features, nlambda = 100){
  n = nrow(Y);
  px = ncol(X);
  py = ncol(Y);
  pz = ncol(Z);
  foldid = fitted.obj$foldid;
  cv.fact_coefs = fitted.obj$factors_coefs;
  cv.projection_coefs = fitted.obj$projection_coefs;
  nfolds = length(unique(foldid))
  #estimate the across validation error for each of the weight
  fact_coefs = fitted.obj$results$post_betas
  projection_coefs = fitted.obj$results$post_bys
  num_patterns = dim(fitted.obj$factors_coefs)[3]
  num_factors = dim(fitted.obj$factors_coefs)[2]
  num_weights = dim(fitted.obj$factors_coefs)[5]
  intercepts = list()
  for(j in 1:py){
    intercepts[[j]] = matrix(NA, nrow = num_weights,ncol = nclasses[j]-1)
  }
  cvm = matrix(NA, nrow = num_weights, ncol = py)
  cvsd =matrix(NA, nrow = num_weights, ncol = py)
  #rescale the overall coefficients
  Yhat = array(NA, dim = c(n, py, num_weights))
  cmin = 0
  Yhat.keep = array(NA, dim = c(n, py, num_weights))
  factors.keep = array(NA, dim = c(n, num_factors, py, num_weights))
  for(l in 1:num_weights){
    for(k in 1:nfolds){
      ucv = array(0, dim = c(n, num_factors))
      for(kk in 1:num_patterns){
        ii = pattern_samples[[kk]]
        jj = pattern_features[[kk]]
        beta = cv.fact_coefs[jj,,kk,k,l]
        ucv[ii,] = Z[ii,jj]%*%beta
      }
      for(j in 1:py){
        factors.keep[foldid==k,,j,l] = t(apply(ucv[foldid==k,], 1, function(z) z*cv.projection_coefs[,j,k,l]))
      }
    }
  }

  for(l in 1:num_weights){
    U0 = array(0, dim = c(n, num_factors))
    for(k in 1:num_patterns){
      ii = pattern_samples[[k]]
      jj = pattern_features[[k]]
      U0[ii,] = Z[ii,jj]%*% fact_coefs[jj,,k,l]
    }
    r2norm = rep(0, py)
    if(py == 1){
      Yhat[,1,l] =  (U0 %*% projection_coefs[,1,l])
      r2norm =  sqrt(mean(Yhat[,1,l]^2))
      if(family == 2){
        Yhat[,1,l] = Yhat[,1,l]/r2norm
      }
    }else{
      Yhat[,,l] =  (U0 %*% projection_coefs[,,l])
      r2norm =  sqrt(apply(Yhat[,,l]^2,2,function(z) mean(z^2)))
      if(family == 2){
        Yhat[,,l] = apply(Yhat[,,l],2,function(z) z/sqrt(mean(z^2)))
      }
    }
    Yhat.cv = array(0, dim = c(n, py, nfolds, num_weights))
    for(k in 1:nfolds){
      ucv = array(0, dim = c(n, num_factors))
      for(kk in 1:num_patterns){
        ii = pattern_samples[[kk]]
        jj = pattern_features[[kk]]
        beta = cv.fact_coefs[jj,,kk,k,l]
        ucv[ii,] = Z[ii,jj]%*%beta
      }
      b = cv.projection_coefs[,,k,l]
      if(py == 1){
        Yhat.cv[,1,k,l] =  (ucv %*% b)
        if(family == 2){
          Yhat.cv[,1,k,l] = Yhat.cv[,1,k,l]/sqrt(mean(Yhat.cv[,1,k,l]^2))
        }
      }else{
        Yhat.cv[,,k,l] =  (ucv %*% b)
        if(family == 2){
        Yhat.cv[,,k,l] =apply(Yhat.cv[,,k,l],2,function(z) z/sqrt(mean(z^2)))
        }
      }
    }
    for(j in 1:py){
      y = Y[,j]
      yhat = Yhat[,j,l]
      if(family == 0){
        tmp = lm(y~Yhat[,j,l])
        cmax = tmp$coefficients[2]
      }else if(nclasses[j] == 2){
        tmp = glm(y ~ yhat, family = "binomial")
        cmax = tmp$coefficients[2]
      }else{
        reverseY = max(y) - y
        reverseY = as.factor(reverseY)
        tmp = polr(reverseY ~ yhat)
        cmax = -tmp$coefficients
        zeta.init = tmp$zeta
      }
      chats = seq(0, cmax, length.out = nlambda)
      errs = array(NA, dim = c(n,length(chats)))
      for(ll in 1:length(chats)){
        for(k in 1:nfolds){
          yhat_cv = Yhat.cv[,j,k,l]
          if(family == 0){
            a = mean(y[foldid!=k] - chats[ll] * yhat_cv[foldid!=k])
            errs[foldid==k,ll] = (y[foldid == k] - a - chats[ll] * yhat_cv[foldid==k])^2
          }else if(nclasses[j] == 2){
            tmp = glm(y[foldid!=k] ~ offset(chats[ll] * yhat_cv[foldid!=k]), family = "binomial")
            a= tmp$coefficients[1]
            Probs = 1/(1+exp(-chats[ll] * yhat_cv[foldid == k]-a))
            errs[foldid==k, ll] = y[foldid == k] * log(Probs+1e-10) +
              (1 - y[foldid == k]) * log(1 - Probs+1e-10);
          }else{
            tmp <-try(polr(reverseY[foldid!=k] ~ offset(-chats[ll] * yhat_cv[foldid !=k]),
                           start = c(zeta.init)))
            if(class(tmp) == "try-error"){
              a = zeta.init
            }else{
              a = tmp$zeta
            }
            a = a[length(a):1]
            ##deviance loss
            Pmat0 = matrix(0, ncol = max(y), nrow = sum(foldid == k))
            Pmat = matrix(0, ncol = max(y)+1, nrow = sum(foldid == k))
            y_extended = matrix(0, ncol = max(y)+1, nrow = sum(foldid == k))
            for(kk in 1:length(a)){
              Pmat0[,kk] = 1/(1+exp(-chats[ll] * yhat_cv[foldid == k]-a[kk]))
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
        for(foldind in 1:nfolds){
          Yhat.keep[foldid == foldind,j,l] = Yhat.cv[foldid == foldind,j,foldind,l]
        }

      factors.keep[,,j,l] = factors.keep[,,j,l] *chats[which.min(cv_tmp)]
      projection_coefs[,j,l] = projection_coefs[,j,l] *chats[which.min(cv_tmp)]
      if(family == 2){
        projection_coefs[,j,l] = projection_coefs[,j,l]/r2norm[j]
      }
      if(family == 0){
        intercepts[[j]][l,] = mean(y - mean(yhat *chats[which.min(cv_tmp)] ))
      }else if(nclasses[j]==2){
        tmp = glm(y ~ offset(chats[which.min(cv_tmp)] * yhat))
        a = tmp$coefficients[1]
        intercepts[[j]][l,] = a
      }else{
        tmp = polr(reverseY ~ offset(-chats[which.min(cv_tmp)] * yhat))
        a = tmp$zeta
        a = a[length(a):1]
        intercepts[[j]][l,] = a
      }
    }
  }
  ###increament contributions
  factor_contributions = array(NA,dim = c(num_factors, py, num_weights))
  for(l in 1:num_weights){
    by = projection_coefs[,,l]
    if(is.null(dim(by))){
      by = matrix(by, ncol = 1)
    }
    U0 = array(0, dim = c(n, num_factors))
    for(k in 1:num_patterns){
      ii = pattern_samples[[k]]
      jj = pattern_features[[k]]
      U0[ii,] = Z[ii,jj]%*% fact_coefs[jj,,k,l]
    }
    if(family == 0){
      for(j in 1:py){
        factor_contributions[,j,l] = apply(factors.keep[,,j,l],2,function(z) var((z-Y[,j])))
        factor_contributions[,j,l] =  1-factor_contributions[,j,l]/var(Y[,j])
      }
      
    }else if(nclasses[j] == 2){
      for(j in 1:py){
        y = Y[,j]
        yhat = Yhat[,j,l]
        null_model = glm(y~offset(rep(0, length(y))), family = "binomial")
        Delta0 = null_model$deviance
        for(k in 1:num_factors){
          yhat1 = factors.keep[,,j,l]
          tmp = glm(y~offset(yhat1), family = "binomial")
          Delta1 = tmp$deviance
          factor_contributions[k,j,l]  =  (Delta0 - Delta1)/Delta0
        }
      }
    }else{
      for(j in 1:py){
        y = Y[,j]
        yhat = Yhat[,j,l]
        reverseY = max(y) - y
        reverseY = as.factor(reverseY)
        null_model = polr(reverseY~offset(rep(0, length(reverseY))))
        Delta0 = null_model$deviance
        for(k in 1:num_factors){
          yhat1 = -factors.keep[,k,j,l]
          tmp = polr(reverseY~offset(yhat1))
          Delta1 = tmp$deviance
          factor_contributions[k,j,l]  = (Delta0 - Delta1)/Delta0
        }
      }
    }
  }
  ##regression coefficients 
  reg_coefs = array(0, dim = c(pz,num_patterns,py,num_weights))
  for(l in 1:num_weights){
    for(j in 1:py){
      for(k in 1:num_patterns){
        reg_coefs[,k,j,l] = fact_coefs[,,k,l]%*%projection_coefs[,j,l]
      }
    }
  }
  return(list(projection_coefs = projection_coefs,
              reg_coefs = reg_coefs,
              intercepts = intercepts,
              cvm = cvm, cvsd = cvsd, 
              factor_contributions = factor_contributions,
              Yhat.keep = Yhat.keep
  ))
}
