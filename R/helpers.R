dataGen_gaussian <- function(N = 500, Ntest = 2000, P = 500, D = 4, seed = 123, num_factors = 5, c = 1, 
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
  if(Ymodel == "factor"){
      data.tr$true_bx =Theta0
      data.tr$true_by =rep(0, num_factors)
      data.tr$true_by[1:2] = 1
      data.tr$true_beta = NULL
  }else{
      data.tr$true_bx =NULL
      data.tr$true_by =NULL
      data.tr$true_beta =beta
  }
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
#'@param X description
#'@param refit Whether to refit the model with elasticNet. If FALSE, the coefficients are from the spear model.
#'@param rule Rule to use for picking cross validation model: "min" or "1se".
#'@param standardize Whether to standardize the data.
#'@param alpha alpha = 1 corresponds to lasso and alpha = 0 corresponds to ridge.
#'@export
preparation <- function(Y,  X, family, pattern_samples = NULL, pattern_assays = NULL,
                        path.type = "assay", bx = NULL,
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
  bx_aggregated = NULL
  if(!is.null(bx)){
      bx_aggregated = bx[[1]]
      for(d in 1:length(X)){
        if(d > 1){
            bx_aggregated = cbind(bx_aggregated, bx[[d]])
        }
      }
      bx_aggregated = t(bx_aggregated)
  }
  functional_path = list()
  if(path.type == "assay"){
    for(i in 1:length(p)){
      if(i == 1){
        functional_path[[i]]  = 1:p[i]
      }else{
        functional_path[[i]] = (sum(p[1:(i-1)])+1):sum(p[1:(i)])
      }
      functional_path[[i]]  = functional_path[[i]]
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
              nclasses = nclasses, bx_aggregated = bx_aggregated))
}



#' Create the data from lists, and create the functional pathway list.
#'@param X description
#'@param Y  description
#'@param Z  description
#'@param Xobs  description
#'@param Yobs  description
#'@param foldid  description
#'@param weights  description
#'@param family  description
#'@param inits.type  description
#'@param num.factors  description
#'@param seed  description
#'@param scale.x  description
#'@param scale.y  description
#'@param num.folds  description
#'@param warmup.iterations  description
#'@param max.iterations  description
#'@param elbo.threshold  description
#'@param elbo.threshold.count  description
#'@param cv.nlambda  description
#'@param print.out  description
#'@param save.model  description
#'@param save.path  description
#'@param save.name  description
#'@param sparsity_upper  sparsity_upper, defaults to .5
#'@param L0 parameter
#'@param factor_contribution Calculate factor contributions? Defaults to TRUE
#'@param run.debug debug?
#'@export
run_cv_spear <- function(X, Y, Z = NULL, Xobs = NULL, Yobs = NULL, foldid = NULL, weights = NULL, family = "gaussian", inits.type = "pca",
                         num.factors = NULL, seed = NULL, scale.x = FALSE, scale.y = FALSE, num.folds = 5, 
                         warmup.iterations = NULL, max.iterations = NULL, elbo.threshold = NULL, elbo.threshold.count = NULL, cv.nlambda = 100, print.out = 100,
                         save.model = TRUE, save.path = NULL, save.name = NULL, run.debug = FALSE, robust_eps = NULL, sparsity_upper = .1, L0 = 1, factor_contribution = TRUE){
  
  success.color <- "light green"
  update.color <- "green"
  info.color <- "light green"
  tilde.color <- "yellow"
  omic.color <- "light cyan"
  response.color <- "light red"
  
  cat("
          __________________________________________
          |                                        |
          |           |''---____                   |
          |     ______|  SPEAR  '''''-----______   |
          |     ''''''|  v.1.0  _____-----''''''   |
          |           |__---''''                   |
          |________________________________________|\n",
 "\nSPEAR version 1.0.   Please direct all questions to Jeremy Gygi \n(", SPEAR.color_text("jeremy.gygi@yale.edu", "green"), ") or Leying Guan (", SPEAR.color_text("leying.guan@yale.edu", "green"), ")\n\n")
  
  
  # Prepare SPEAR object
  cat("******************\n", SPEAR.color_text("Loading Datasets", info.color), "\n******************\n")
  
  if(family == "gaussian"){
    family.encoded <- 0
  } else if(family == "binomial"){
    family.encoded <- 1
  } else if(family == "ordinal"){
    family.encoded <- 2
  } else if(family == "multinomial"){
    family.encoded <- 3
  }else {
    cat("***", paste0(" 'family' parameter not recognized (", family, "). Assuming 'gaussian'. Acceptable values are 'gaussian', 'binomial', and 'multinomial'.\n"))
    family.encoded <- 0
    family <- "gaussian"
  }
  
  cat(SPEAR.color_text("\nPreparing omics data...\n", update.color))
  
  # Make sure X is a list:
  if(scale.x){
    X.scaled <- lapply(X, scale)
  } else {
    X.scaled <- X
  }
  if(is.null(names(X.scaled))){
    cat("***", paste0(" Names for datasets in X not provided. Renaming to ", paste(paste0("X", 1:length(X.scaled)), collapse = ", "), "\n"))
    names(X.scaled) <- paste0("X", 1:length(X.scaled))
  }
  for(d in 1:length(X.scaled)){
    if(is.null(colnames(X.scaled[[d]]))){
      cat("***", paste0(" Feature names in ", names(X.scaled)[d], " not provided. Renaming to ", paste0(names(X.scaled)[d], "_feat", 1), " ... ", paste0(names(X.scaled)[d], "_feat", ncol(X.scaled[[d]])), "\n"))
      colnames(X.scaled[[d]]) <- paste0(names(X.scaled)[d], "_feat", 1:ncol(X.scaled[[d]]))
    }
    if(is.null(rownames(X.scaled[[d]]))){
      cat("***", paste0(" Sample names in ", names(X.scaled)[d], " not provided. Renaming to sample1... sample", nrow(X.scaled[[d]]), "\n"))
      rownames(X.scaled[[d]]) <- paste0("sample", 1:nrow(X.scaled[[d]]))
    }
  }
  cat(paste0("Detected ", length(X.scaled), " datasets:\n"))
  for(i in 1:length(X.scaled)){
    cat(SPEAR.color_text(names(X.scaled)[i], omic.color), "\tSubjects: ", nrow(X.scaled[[i]]), "\tFeatures: ", ncol(X.scaled[[i]]), "\n")
  }
  
  if(scale.y & family.encoded == 0){
    Y.scaled <- scale(Y)
  } else {
    Y.scaled <- Y
  }
  if(is.null(colnames(Y.scaled))){
    cat("***", paste0(" Column names for response Y not provided. Renaming to ", paste(paste0("Y", 1:ncol(Y.scaled)), collapse = ", "), "\n"))
    colnames(Y.scaled) <- paste0("Y", 1:ncol(Y.scaled))
  }
  if(is.null(rownames(Y.scaled))){
    cat("***", paste0(" Row names for response Y not provided. Renaming to sample1... sample", nrow(Y.scaled), "\n"))
    rownames(Y.scaled) <- paste0("Y", 1:ncol(Y.scaled))
  }
  cat(paste0("Detected ", ncol(Y.scaled), " response ", ifelse(ncol(Y.scaled) == 1, "variable", "variables"), ":\n"))
  for(i in 1:ncol(Y.scaled)){
    cat(SPEAR.color_text(colnames(Y.scaled)[i], response.color), "\tSubjects: ", sum(!is.na(Y.scaled[,i])), "\tType: ", family, "\n")
  }
  
  # Quickly ensure that each sample is lined up:
  for(d in 1:length(X.scaled)){
    if(any(rownames(X.scaled[[d]]) != rownames(Y.scaled))){
      stop("ERROR: Rownames for ", names(X.scaled)[d], " are not equal to the others (at least to the response Y). Ensure they are all the same and try again.")
    }
  }
  
  
  cat("\n")
  
  # Run Preparation Function:
  data <- preparation(Y = Y.scaled, X = X.scaled, family = family.encoded)
  data$xlist <- X.scaled
  
  # Parameters:
  cat(SPEAR.color_text("Preparing SPEAR parameters...\n", update.color))
  
  # Seed:
  if(is.null(seed)){
    cat("***", paste0(" seed not provided. Consider using a seed (i.e. 123) for reproducibility.\n"))
  } else {
    cat(SPEAR.color_text("~~~", tilde.color), paste0(" seed set to ", seed, ".\n"))
    set.seed(seed)
  }
  # Set up matrices that indicate which samples are missing (since none are missing, just fill with 1's)
  if(is.null(Xobs)){
    cat("***", paste0(" Xobs not provided. Assuming full observations in X.\n"))
    Xobs <- array(1, dim  = dim(data$X))
  }
  if(is.null(Yobs)){
    cat("***", paste0(" Yobs not provided. Assuming full observations in Y.\n"))
    Yobs <- array(1, dim  = dim(data$Y))
  }
  # Assigning fold-id's to each of the N subjects in the training set:
  if(is.null(foldid)){
    cat("***", paste0(" foldid not provided. Assigning folds randomly from 1 to ", num.folds, ".\n"))
    foldid = sample(1:num.folds, nrow(data$X), replace = T)
  } else {
    if(length(foldid) == nrow(data$X)){
      cat(SPEAR.color_text("~~~", tilde.color), paste0(" foldid provided | passed\n"))
    } else {
      stop(paste0("SPEAR preparation ERROR: length(foldid) = ", length(foldid), " != Number of Subjects (", nrow(data$X), "). Cannot use provided foldids.\n"))
    }
  }
  if(is.null(weights)){
    cat("***", paste0(" weights not provided. Using c(2, 1.5, 1.0, 0.5, 0.0)\n"))
    weights <- c(2, 1.5, 1.0, 0.5, 0.0)
  } else {
    cat(SPEAR.color_text("~~~", tilde.color), paste0(" weights - using c(", paste(round(weights, 2), collapse = ", "), ")\n"))
  }
  if(is.null(num.factors)){
    cat("***", paste0(" num.factors not provided. Defaulting to K = 5 factors\n"))
    num.factors <- 5
  } else {
    cat(SPEAR.color_text("~~~", tilde.color), paste0(" num.factors - using K = ", num.factors, " factors\n"))
  }
  
  if(is.null(Z)){
    Z <- data$X
  }
  
  other.params <- c("warmup", "max.iterations", "elbo.threshold", "elbo.threshold.count")
  if(is.null(warmup.iterations)){
    warm_up <- 50
  } else {
    warm_up <- warmup.iterations
  }
  if(is.null(max.iterations)){
    max_iter <- 1000
  } else {
    max_iter <- max.iterations
  }
  if(is.null(elbo.threshold)){
    thres_elbo <- 0.1
  } else {
    thres_elbo <- elbo.threshold
  }
  if(is.null(elbo.threshold.count)){
    thres_count <- 5
  } else {
    thres_count <- elbo.threshold.count
  }
  if(is.null(robust_eps)){
    robust_eps = 1.0/sqrt(nrow(data$X))
  }
  
  # Debug statements to print:
  if(run.debug){
    for(k in 1:num.folds){
      print(table(data$Y[foldid==k]))
    }
    print(foldid)
  }
  
  
  
  cat("\n*****************\n", SPEAR.color_text(" Running SPEAR", info.color), "\n*****************\n")
  
  
  cat(SPEAR.color_text("Starting parallel workers...\n", success.color))
  if(.Platform$OS.type == "windows"){
    cat(SPEAR.color_text("*NOTE:* Windows machine detected. SPEAR uses the mclapply function for parallelization, which is not supported on Windows.
        Consider using SPEAR on a unix operating system for boosted performance. Setting numCores = 1\n", "red"))
    numCores <- 1
  } else {
    numCores <- parallel::detectCores()
  }
  
  cat("NOTE: Only printing out results from one fold of the CV for simplicity...\n")
  
  spear_fit <- cv.spear(X = as.matrix(data$X),
                        Y = as.matrix(data$Y),
                        Xobs = Xobs, 
                        Yobs = Yobs,
                        Z = Z, 
                        pattern_samples = data$pattern_samples,
                        pattern_features = data$pattern_features,
                        family = family.encoded, 
                        nclasses = data$nclasses,
                        ws_x = weights,
                        foldid = foldid, 
                        num_factors = num.factors,  
                        functional_path = data$functional_path,
                        print_out = print.out, 
                        inits_type = inits.type, 
                        warm_up= warm_up, 
                        max_iter = max_iter, 
                        seed = seed,
                        sparsity_upper = sparsity_upper, 
                        L = nrow(data$X)/L0,
                        thres_elbo = thres_elbo, 
                        thres_count = thres_count,
                        crossYonly = F,
                        numCores = numCores,
                        run.debug = run.debug,
                        a0 = 1e-2, b0 = 1e-2,
                        a1 = sqrt(nrow(data$X)), b1 =sqrt(nrow(data$X)),
                        a2= sqrt(nrow(data$X)), b2 = sqrt(nrow(data$X)))
  
  # Run cv.eval:
  
  cat(SPEAR.color_text("\nFinished running SPEAR.\n", success.color))
  cat("\n*****************\n", SPEAR.color_text("Evaluating CV:", info.color), "\n*****************\n")
  
  cv.eval <- cv.evaluation(fitted.obj = spear_fit, 
                           X = data$X, 
                           Y = data$Y, 
                           Z = Z, 
                           family = family.encoded, 
                           nclasses = data$nclasses, 
                           pattern_samples = data$pattern_samples, 
                           pattern_features = data$pattern_features,
                           factor_contribution = factor_contribution,
                           nlambda = cv.nlambda)
  
  # Return a SPEARobject:
  params <- list()
  params$weights <- weights
  params$foldid <- foldid
  params$seed <- seed
  params$num_factors = num.factors
  params$max_iter = max_iter
  params$thres_elbo = thres_elbo
  params$thres_count = thres_count
  params$family <- family
  params$sparsity_upper <- sparsity_upper
  # Add colors:
  colorlist <- list()
  num.omics <- length(data$xlist)
  num.resp <- ncol(data$Y)
  cols <- c("#9E0142", "#F46D43", "#ABDDA4", "#3288BD", "#5E4FA2")
  omic.cols <- grDevices::colorRampPalette(cols)(num.omics)
  names(omic.cols) <- names(data$xlist)
  cols <- c("#F8AFA8", "#D0C7E7", "#DBD7AF", "#F1AD75", "#E47A7F", "#85D4E3")
  resp.cols <- grDevices::colorRampPalette(cols)(num.resp)
  names(resp.cols) <- colnames(data$Y)
  params$colors <- list(X = omic.cols, Y = resp.cols)
  
  SPEARobj <- list(fit = spear_fit,
                   cv.eval = cv.eval,
                   params = params,
                   data = data)
  
  cat(SPEAR.color_text("\nSPEAR has been run successfully.\n\n", success.color))
  
  # Save SPEAR results:
  if(save.model){
    # Check save.path:
    if(is.null(save.path)){
      dir.create(paste0("SPEAR_model_", format(Sys.time(), "%m%d%Y_%H_%M_%S")))
      cat("***", paste0(" No directory specificied for saving SPEAR results. Generated temporary folder at:\n", getwd(), "/SPEAR_model_", format(Sys.time(), "%m%d%Y_%H_%M_%S"), "/\n"))
      save.path <- paste0(getwd(), "/SPEAR_model_", format(Sys.time(), "%m%d%Y_%H_%M_%S"), "/")
    } else if(save.path == ""){
      save.path <- paste0(getwd(), "/")
    } else if(!endsWith(save.path, "/")){
      save.path <- paste0(save.path, "/")
    }
    if(is.null(save.name)){
      cat("***", paste0(" No filename specified for saving SPEAR results. Saving as ", paste0(save.path, "SPEARobject_", format(Sys.time(), "%m%d%Y_%H_%M_%S"), ".rds"), "\n"))
      save.name <- paste0(save.path, "SPEARobject_", format(Sys.time(), "%m%d%Y_%H_%M_%S"), ".rds")
    } else if(!endsWith(save.name, ".rds")){
      save.name <- paste0(save.name, ".rds")
    }
    saveRDS(SPEARobj, file = paste0(save.path, save.name))
    cat(paste0("Saved SPEARobject to RDS file at ", save.name, "\n\n"))
  }
  
  return(SPEARobj)
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
                          pattern_samples, pattern_features, nlambda = 100,
                          factor_contribution = F, weights = NULL, max_iter = 1e4,
                          multinomical_loss = "deviance"){
  n = nrow(Y);
  px = ncol(X);
  py = ncol(Y);
  pz = ncol(Z);
  standardize_family = c(2,3)
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
  
  chats.return = list()
  #rescale the overall coefficients
  Yhat = array(NA, dim = c(n, py, num_weights))
  cmin = 0
  Yhat.keep = array(NA, dim = c(n, py, num_weights))
  factors.keep = array(NA, dim = c(n, num_factors, py, num_weights))
  Uhat.cv = array(0, dim = c(n,num_factors,num_weights))
  Yhat.cv = array(0, dim = c(n, py, nfolds, num_weights))
  for(l in 1:num_weights){
    for(k in 1:nfolds){
      ucv = array(0, dim = c(n, num_factors))
      for(kk in 1:num_patterns){
        ii = pattern_samples[[kk]]
        jj = pattern_features[[kk]]
        beta = cv.fact_coefs[jj,,kk,k,l]
        ucv[ii,] = Z[ii,jj]%*%beta
      }
      Uhat.cv[foldid==k,,l] = ucv[foldid==k,]
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
      chats.return.temp <- list()
      for(j in 1:py){
        ##note that scaling is only required for Gaussian and logistic!
        y = Y[,j]
        yhat = Yhat[,j,l]
        if(family == 0){
          tmp = lm(y~Yhat[,j,l])
          cmax = tmp$coefficients[2]
        }else if(family==1){
          tmp = glmnet(cbind(yhat,rep(0,length(yhat))),y, family = "binomial",lambda = 1e-4)
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
            tmp0 = glmnet(cbind(yhat_cv[foldid!=k],rep(0,sum(foldid!=k))),y[foldid!=k], family = "binomial",lambda = 1e-4)
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
        for(foldind in 1:nfolds){
          Yhat.keep[foldid == foldind,j,l] = Yhat.cv[foldid == foldind,j,foldind,l]
        }
        projection_coefs[,j,l] = projection_coefs[,j,l] *chats[which.min(cv_tmp)]
        if(family %in% standardize_family){
          projection_coefs[,j,l] = projection_coefs[,j,l]/r2norm[j]
        }
        chats.return.temp[[j]] <- list(chats = chats, minval = which.min(cv_tmp))
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
      fitted_multinomial = glmnet(x =U0std, y = ycollapsed, standardize = F, family = 'multinomial', nlambda = nlambda)
      lambdas =fitted_multinomial$lambda 
      # perform cv fit with given penalty lists
      probs_predictions = array(NA, dim  = c(nrow(Y),ncol(Y),length(lambdas)))
      tmp.fit = list()
      for(fold_id in 1:nfolds){
        test_id = which(foldid == fold_id)
        train_id = which(foldid != fold_id)
        U0cv = array(0, dim = c(n, num_factors))
        for(k in 1:num_patterns){
          ii = pattern_samples[[k]]
          jj = pattern_features[[k]]
          #' This is where the error occurs previously!
          #' changed from U0cv * scale.factors[,fold_id,l] to 
          #'  U0cv %*%diag(scale.factors[,fold_id,l])
          U0cv[ii,] = Z[ii,jj]%*% cv.fact_coefs[jj,,k,fold_id, l]
        }
        if(family %in% standardize_family){
          tmp = sqrt(apply(U0cv,2,function(z) mean(z^2)))
          scale.factors[,fold_id,l] = r2norm/tmp
          U0cv.std  = U0cv %*%diag(scale.factors[,fold_id,l])
        }else{
          U0cv.std  = U0cv
        }
        tmp.fit[[fold_id]] = glmnet(x =U0cv.std[train_id,], y = ycollapsed[train_id], standardize = F, family = 'multinomial', lambda = lambdas)
        preds = predict(tmp.fit[[fold_id]], U0cv.std)
        total = apply(preds,c(1,3),function(z) matrixStats::logSumExp(z))
        for(kkk in 1:dim(preds)[2]){
          preds[,kkk,] =  preds[,kkk,] -total
        }
        preds = exp(preds)
        probs_predictions[test_id,,1:dim(preds)[3]] = preds[test_id,,]
      }
      if(multinomical_loss=="deviance"){
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
  reg_coefs = array(0, dim = c(pz,num_patterns,py,num_weights))
  for(l in 1:num_weights){
    for(j in 1:py){
      for(k in 1:num_patterns){
        reg_coefs[,k,j,l] = fact_coefs[,,k,l]%*%projection_coefs[,j,l]
      }
    }
  }
  #replace deviance contriution to spearman correlation
  factor_contributions = array(NA,dim = c(num_factors, py, num_weights))
  factor_contributions_pvals = array(NA,dim = c(num_factors, py, num_weights))
  if(factor_contribution){
    for(l in 1:num_weights){
      if(!is.null(weights)){
        print(paste0("~~~ Calculating contribution for w = ",weights[l], " | ", l, "/", num_weights))
      } else {
        print(paste0("~~~ Calculating contribution for weights: ", l, "/", num_weights))
      }
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
  
  
  return(list(projection_coefs_scaled = projection_coefs,
              cv.projection_coefs_scaled = cv.projection_coefs,
              reg_coefs = reg_coefs,
              intercepts = intercepts, chats = chats.return,
              cvm = cvm, cvsd = cvsd, 
              factor_contributions = factor_contributions,
              factor_contributions_pvals = factor_contributions_pvals,
              Yhat.keep = Yhat.keep))
}



