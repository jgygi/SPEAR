### SPEARkit Functions ###

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# RUN SPEAR COMBO
#
run_cv_spear <- function(X, Y, Z = NULL, Xobs = NULL, Yobs = NULL, foldid = NULL, weights = NULL, family = 0, inits.type = "pca",
                         num.factors = NULL, seed = NULL, scale.x = TRUE, scale.y = TRUE, num.folds = 5, 
                         warmup.iterations = NULL, max.iterations = NULL, elbo.threshold = NULL, elbo.threshold.count = NULL, cv.nlambda = 100, print.out = 100,
                         save.model = TRUE, save.path = NULL, save.name = NULL){
  cat("
          __________________________________________
          |                                        |
          |           |''---____                   |
          |     ______|  SPEAR  '''''-----______   |
          |     ''''''|  v.1.0  _____-----''''''   |
          |           |__---''''                   |
          |________________________________________|
      
 SPEAR version 1.0. Please direct all questions to Jeremy Gygi (jeremy.gygi@yale.edu) or Leying Guan (leying.guan@yale.edu).\n\n")


  # Prepare SPEAR object
  cat("*****************\n Preparing SPEAR\n*****************\n")
  
  cat("\nPreparing omics data...\n")
  
  # Make sure X is a list:
  if(scale.x){
    X.scaled <- lapply(X, scale)
  } else {
    X.scaled <- X
  }
  if(is.null(names(X.scaled))){
    cat(paste0("*** Names for datasets in X not provided. Renaming to ", paste(paste0("X", 1:length(X.scaled)), collapse = ", "), "\n"))
    names(X.scaled) <- paste0("X", 1:length(X.scaled))
  }
  for(d in 1:length(X.scaled)){
    if(is.null(colnames(X.scaled[[d]]))){
      cat(paste0("*** Feature names in ", names(X.scaled)[d], " not provided. Renaming to ", paste0(names(X.scaled)[d], "_feat", 1), " ... ", paste0(names(X.scaled)[d], "_feat", ncol(X.scaled[[d]])), "\n"))
      colnames(X.scaled[[d]]) <- paste0(names(X.scaled)[d], "_feat", 1:ncol(X.scaled[[d]]))
    }
  }
  cat(paste0("Detected ", length(X.scaled), " datasets:\n"))
  for(i in 1:length(X.scaled)){
    cat(paste0(names(X.scaled)[i], "\tSubjects: ", nrow(X.scaled[[i]]), "\tFeatures: ", ncol(X.scaled[[i]]), "\n"))
  }
  
  if(scale.y){
    Y.scaled <- scale(Y)
  } else {
    Y.scaled <- Y
  }
  if(is.null(colnames(Y.scaled))){
    cat(paste0("*** Names for response Y not provided. Renaming to ", paste(paste0("Y", 1:ncol(Y.scaled)), collapse = ", "), "\n"))
    colnames(Y.scaled) <- paste0("Y", 1:ncol(Y.scaled))
  }
  cat(paste0("Detected ", ncol(Y.scaled), " response ", ifelse(ncol(Y.scaled) == 1, "variable", "variables"), ":\n"))
  for(i in 1:ncol(Y.scaled)){
    cat(paste0(colnames(Y.scaled)[i], "\tSubjects: ", sum(!is.na(Y.scaled[,i])), "\tType: ", ifelse(family == 0, "Gaussian", "Ordinal"), "\n"))
  }
  
  
  
  
  cat("\n")

  # Run Preparation Function:
  data <- SPEARcomplete::preparation(Y = Y.scaled, X = X.scaled, family = family)
  data$xlist <- X.scaled
  
  # Parameters:
  cat(paste0("Preparing SPEAR parameters...\n"))
  
  # Seed:
  if(is.null(seed)){
    cat(paste0("*** seed not provided. Consider using a seed (i.e. 123) for reproducibility.\n"))
  } else {
    cat(paste0("~~~ seed set to ", seed, ".\n"))
    set.seed(seed)
  }
  # Set up matrices that indicate which samples are missing (since none are missing, just fill with 1's)
  if(is.null(Xobs)){
    cat(paste0("*** Xobs not provided. Assuming full observations in X.\n"))
    Xobs <- array(1, dim  = dim(data$X))
  }
  if(is.null(Yobs)){
    cat(paste0("*** Yobs not provided. Assuming full observations in Y.\n"))
    Yobs <- array(1, dim  = dim(data$Y))
  }
  # Assigning fold-id's to each of the N subjects in the training set:
  if(is.null(foldid)){
    cat(paste0("*** foldid not provided. Assigning folds randomly from 1 to ", num.folds, ".\n"))
    foldid = sample(1:num.folds, nrow(data$X), replace = T)
  } else {
    if(length(foldid) == nrow(data$X)){
      cat(paste0("~~~ foldid provided | passed\n"))
    } else {
      cat(paste0("~~~ foldid provided | ERROR: length(foldid) = ", length(foldid), " != Number of Subjects (", nrow(data$X), "). Cannot use provided foldids.\n"))
      stop(paste0("SPEAR preparation ERROR: length(foldid) = ", length(foldid), " != Number of Subjects (", nrow(data$X), "). Cannot use provided foldids.\n"))
    }
  }
  
  if(is.null(weights)){
    cat(paste0("*** weights not provided. Using c(2, 1.8, 1.6, 1.4, 1.2, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0)\n"))
    weights <- c(2, 1.8, 1.6, 1.4, 1.2, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0)
  } else {
    cat(paste0("~~~ weights - using c(", paste(round(weights, 2), collapse = ", "), ")\n"))
  }
  if(is.null(num.factors)){
    cat(paste0("*** num.factors not provided. Defaulting to K = 5 factors\n"))
    num.factors <- 5
  } else {
    cat(paste0("~~~ num.factors - using K = ", num.factors, " factors\n"))
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
  
  
  
  cat("\n*****************\n  Running SPEAR\n*****************\n\n")
  
  if(.Platform$OS.type == "windows"){
    cat(paste0("*NOTE:* Windows machine detected. SPEAR uses the mclapply function for parallelization, which is not supported on Windows.
        Consider using SPEAR on a unix operating system for boosted performance.\n"))
    numCores <- 1
  } else {
    numCores <- parallel::detectCores()
  }
  
  # Run cv.spear:
  spear_fit <- cv.spear(X = data$X, 
                        Y = data$Y,
                        Xobs = Xobs, 
                        Yobs = Yobs,
                        Z = Z, 
                        pattern_samples = data$pattern_samples,
                        pattern_features = data$pattern_features,
                        family = family, 
                        nclasses = data$nclasses,
                        ws = weights,
                        foldid = foldid, 
                        num_factors = num.factors,  
                        functional_path = data$functional_path,
                        print_out = print.out, 
                        inits_type = inits.type, 
                        warm_up= warm_up, 
                        max_iter = max_iter, 
                        seed = seed, 
                        thres_elbo = thres_elbo, 
                        thres_count = thres_count, 
                        crossYonly = F,
                        numCores = numCores)
  
  # Run cv.eval:
  
  cat("\nFinished running SPEAR.\n")
  cat("\n*****************\n  Evaluating CV:\n*****************\n")
  
  cv.eval <- cv.evaluation(fitted.obj = spear_fit, 
                                        X = data$X, 
                                        Y = data$Y, 
                                        Z = Z, 
                                        family = family, 
                                        nclasses = data$nclasses, 
                                        pattern_samples = data$pattern_samples, 
                                        pattern_features = data$pattern_features, 
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
  
  SPEARobj <- list(fit = spear_fit,
                   cv.eval = cv.eval,
                   params = params,
                   data = data)
  
  cat("\nSPEAR has been run successfully.\n\n")
  
  # Save SPEAR results:
  if(save.model){
    # Check save.path:
    if(is.null(save.path)){
      dir.create(paste0("SPEAR_model_", format(Sys.time(), "%m%d%Y_%H_%M_%S")))
      cat(paste0("*** No directory specificied for saving SPEAR results. Generated temporary folder at:\n", getwd(), "/SPEAR_model_", format(Sys.time(), "%m%d%Y_%H_%M_%S"), "/\n"))
      save.path <- paste0(getwd(), "/SPEAR_model_", format(Sys.time(), "%m%d%Y_%H_%M_%S"), "/")
    } else if(save.path == ""){
      save.path <- paste0(getwd(), "/")
    } else if(!endsWith(save.path, "/")){
      save.path <- paste0(save.path, "/")
    }
    if(is.null(save.name)){
      cat(paste0("*** No filename specified for saving SPEAR results. Saving as ", paste0(save.path, "SPEARobject_", format(Sys.time(), "%m%d%Y_%H_%M_%S"), ".rds"), "\n"))
      save.name <- paste0(save.path, "SPEARobject_", format(Sys.time(), "%m%d%Y_%H_%M_%S"), ".rds")
    } else if(!endsWith(save.name, ".rds")){
      save.name <- paste0(save.name, ".rds")
    }
    saveRDS(SPEARobj, file = paste0(save.path, save.name))
    cat(paste0("Saved SPEARobject to RDS file at ", paste0(save.path, save.name)))
  }
  
  return(SPEARobj)
}


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# GET NEW SPEAR PREDICTIONS
#
SPEAR.predict_new_samples <- function(SPEARobj, X, scale.x = TRUE){
  w.idx <- which.min(apply(SPEARobj$cv.eval$cvm,1,sum))
  # Quickly check that the dimensions in X match:
  if(length(X) != length(SPEARobj$data$xlist)){
    stop(paste0("ERROR: Object 'X' passed in has ", length(X), " datasets, whereas the training data for this SPEAR object contained ", length(SPEARobj$data$xlist), " datasets. Incompatible dimensions."))
  }
  for(d in 1:length(X)){
    if(ncol(X[[d]]) != ncol(SPEARobj$data$xlist[[d]])){
      stop(paste0("ERROR: Incompatible dimensions. New dataset ", d, " has ", ncol(X[[d]]), " features, whereas training dataset ", d, " has ", ncol(SPEARobj$data$xlist[[d]]), " features. Cannot predict new responses."))
    }
  }
  # Collapse X:
  xlist.te <- do.call("cbind", X)
  # Scale X by feature:
  if(scale.x){
    xlist.te <- scale(xlist.te)
  }
  # Right now, assumes 1 group...
  g <- 1
  preds <- xlist.te %*% SPEARobj$cv.eval$reg_coefs[,g,,w.idx] + SPEARobj$cv.eval$intercepts[[1]][w.idx,]
  colnames(preds) <- colnames(SPEARobj$data$Y)
  rownames(preds) <- rownames(xlist.te)
  res <- list(predictions = preds, w = SPEARobj$params$weights[w.idx])
  return(res)
}
#SPEAR.predict_new_samples(SPEARobj, SPEARobj$data$xlist)
#stop("test")


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# GET BEST SPEAR WEIGHTS:
# include.best.overall - add a final element to the named vector for the best overall weight
SPEAR.get_best_weights <- function(SPEARobj, include.best.overall = FALSE){
  indices = apply(SPEARobj$cv.eval$cvm, 2, which.min)
  ws <- SPEARobj$params$weights[indices]
  names(ws) <- colnames(SPEARobj$data$Y)
  if(include.best.overall){
    idx = which.min(apply(SPEARobj$cv.eval$cvm,1,sum))
    ws <- c(ws, SPEARobj$params$weights[idx])
    names(ws)[length(ws)] <- "best.overall.weight"
  }
  return(ws)
}

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# GET FACTOR CONTRIBUTIONS:
# Threshold: minimum value for a factor contribution to be "relevant"
SPEAR.get_factor_contributions <- function(SPEARobj, threshold = .01){
  factor.contributions <- list()
  ws <- apply(SPEARobj$cv.eval$cvm, 2, which.min)
  relevant.factors = matrix(0,ncol = length(ws), nrow = dim(SPEARobj$cv.eval$factor_contributions)[1])
  for(i in 1:length(ws)){
    factor.contributions[[i]] <- SPEARobj$cv.eval$factor_contributions[, i, ws[i]]
    relevant.factors[factor.contributions[[i]] >= threshold, i] <- 1
  }
  factor.contributions <- matrix(unlist(factor.contributions), nrow = dim(SPEARobj$cv.eval$factor_contributions)[1])
  colnames(factor.contributions) <- colnames(SPEARobj$data$Y)
  rownames(factor.contributions) <- paste0("Factor", 1:nrow(factor.contributions))

  return(list(factor.contributions = factor.contributions, relevant.factors = relevant.factors, threshold = threshold))
}


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# PLOT FACTOR CONTRIBUTIONS:
# Threshold: minimum value for a factor contribution to be "relevant"
# show.labels - should the contribution values be shown on the graph?
SPEAR.plot_factor_contributions <- function(SPEARobj, threshold = .01, show.labels = TRUE, show.irrelevant = FALSE){
  factor.contributions <- list()
  ws <- apply(SPEARobj$cv.eval$cvm, 2, which.min)
  relevant.factors = matrix(0,ncol = length(ws), nrow = dim(SPEARobj$cv.eval$factor_contributions)[1])
  for(i in 1:length(ws)){
    factor.contributions[[i]] <- SPEARobj$cv.eval$factor_contributions[, i, ws[i]]
    relevant.factors[factor.contributions[[i]] >= threshold, i] <- 1
  }
  factor.contributions <- matrix(unlist(factor.contributions), nrow = dim(SPEARobj$cv.eval$factor_contributions)[1])
  colnames(factor.contributions) <- colnames(SPEARobj$data$Y)
  rownames(factor.contributions) <- paste0("Factor", 1:nrow(factor.contributions))
  
  if(show.irrelevant){
    relevant.factors <- matrix(1, ncol = length(ws), nrow = dim(SPEARobj$cv.eval$factor_contributions)[1])
  }
  df <- reshape2::melt(factor.contributions * relevant.factors)
  g <- ggplot(data = df) +
    geom_tile(aes(y = Var2, x = Var1, fill = value, alpha = abs(value)), color = "black") +
    geom_text(aes(y = Var2, x = Var1, label = round(ifelse(value==0, NA, value), 2)), color = ifelse(show.labels, "black", NA)) +
    scale_fill_gradient(low = "white", high = "#74149E") +
    scale_alpha_continuous(guide = F) +
    ggtitle("Factor Contribution") +
    ylab("") +
    xlab("") +
    coord_equal() +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 12))
  return(g)
}


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# GET CV PREDICTION ERRORS:
# Verbose: return a matrix with both best values of w AND cv.errors? Or just errors?
SPEAR.get_cv_prediction_error <- function(SPEARobj, show.all.w = FALSE, verbose = FALSE){
  if(show.all.w){
    cv.errors <- SPEARobj$cv.eval$cvm
    cv.errors <- cbind(cv.errors, SPEARobj$params$weights)
    colnames(cv.errors) <- c(colnames(SPEARobj$data$Y), "w")
    return(cv.errors)
  }
  else if(verbose){
    cv.errors <- apply(SPEARobj$cv.eval$cvm,2,min)
    cv.error.mat <- matrix(0, nrow = 2, ncol = length(cv.errors))
    cv.error.mat[1,] <- SPEARobj$params$weights[apply(SPEARobj$cv.eval$cvm, 2, which.min)]
    cv.error.mat[2,] <- cv.errors
    rownames(cv.error.mat) <- c("best.w", "cv.error")
    colnames(cv.error.mat) <- colnames(SPEARobj$data$Y)
    return(cv.error.mat)
  }
  else{
    cv.errors <- apply(SPEARobj$cv.eval$cvm,2,min)
    names(cv.errors) <- colnames(SPEARobj$data$Y)
    return(cv.errors)
  }
}


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# PLOT FACTOR LOADINGS:
# Threshold: minimum value for a factor contribution to be "relevant"
# plot.per.omic: should each omic be a different plot?
# show.feature.names: show feature names (on Y axis)? Recommended as FALSE for large datasets
# return.list - return a list of all plots instead of a cowplot
SPEAR.plot_factor_loadings <- function(SPEARobj, threshold = .01, plot.per.omic = FALSE, return.list = FALSE, show.feature.names = FALSE, plot.irrelevant.factors = FALSE){
  plot.list <- list()
  ws <- apply(SPEARobj$cv.eval$cvm, 2, which.min)
  for(k in 1:ncol(SPEARobj$data$Y)){
    relevant.factors <- SPEARobj$cv.eval$factor_contributions[, k, ws[k]] >= threshold
    names(relevant.factors) <- paste0("Factor", 1:length(relevant.factors))
    w.loadings = SPEARobj$fit$results$post_selections[,,ws[k]]
    if(!plot.irrelevant.factors){
      w.loadings[,!relevant.factors] <- NA
    }
    colnames(w.loadings) <- paste0("F", 1:ncol(w.loadings))
    rownames(w.loadings) <- colnames(SPEARobj$data$X)
    # Generate plots
    if(!plot.per.omic){
      # Plot all together
      df <- reshape2::melt(w.loadings)
      colnames(df)[3] <- "Probability"
      g.temp <- ggplot(df) +
        geom_tile(aes(x = Var2, y = Var1, fill = Probability)) +
        scale_fill_distiller(palette = "Reds", direction = 1, na.value="#F8F8F8") +
        ggtitle(paste0("Factor Loadings | ", colnames(SPEARobj$data$Y)[k])) +
        ylab("Features") +
        xlab("") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 12))
      
      if(!show.feature.names){
        g.temp <- g.temp + theme(axis.text.y=element_blank(),
                                 axis.ticks.y=element_blank())
      }
        
      # Save
      plot.list[[length(plot.list) + 1]] <- g.temp
    }
    else{
      # Split w.loadings into separate omics
      num.omics <- length(SPEARobj$data$xlist)
      s <- 1
      for(o in 1:num.omics){
        w.loadings.current <- w.loadings[s:(s + ncol(SPEARobj$data$xlist[[o]]) - 1),]
        s <- s + ncol(SPEARobj$data$xlist[[o]])
        # Make a separate plot for each:
        df <- reshape2::melt(w.loadings.current)
        colnames(df)[3] <- "Probability"
        g.temp <- ggplot(df) +
          geom_tile(aes(x = Var2, y = Var1, fill = Probability)) +
          scale_fill_distiller(palette = "Reds", direction = 1, na.value="#F8F8F8") +
          ggtitle(paste0("Factor Loadings | ", colnames(SPEARobj$data$Y)[k])) +
          ylab("Features") +
          xlab("") +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5, size = 12))
        
        if(!show.feature.names){
          g.temp <- g.temp + theme(axis.text.y=element_blank(),
                                   axis.ticks.y=element_blank())
        }
        
        # Save
        plot.list[[length(plot.list) + 1]] <- g.temp
      }
    }
  }
  # Return plot or list:
  if(return.list){
    return(plot.list)
  }
  else{
    if(plot.per.omic){
      p <- cowplot::plot_grid(plotlist = plot.list, ncol = ncol(SPEARobj$data$Y), byrow = FALSE)
    }
    else{
      p <- cowplot::plot_grid(plotlist = plot.list, ncol = length(plot.list))
    }
    return(p)
  }
}


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# GET FEATURES:
# Threshold: minimum value for a factor contribution to be "relevant"
# cutoff - posterior selection probability cutoff (defaults to .5)
SPEAR.get_factor_features <- function(SPEARobj, cutoff = .5, threshold = .01){
  feature.list <- list()
  ws <- apply(SPEARobj$cv.eval$cvm, 2, which.min)
  for(k in 1:ncol(SPEARobj$data$Y)){
    features.temp <- list()
    relevant.factors <- SPEARobj$cv.eval$factor_contributions[, k, ws[k]] >= threshold
    names(relevant.factors) <- paste0("Factor", 1:length(relevant.factors))
    w.loadings = SPEARobj$fit$results$post_selections[,,ws[k]]
    w.loadings[,!relevant.factors] <- NA
    colnames(w.loadings) <- paste0("F", 1:ncol(w.loadings))
    rownames(w.loadings) <- colnames(SPEARobj$data$X)
    # by factor
    if(sum(relevant.factors) > 0){
      factor.features <- list()
      for(f in 1:sum(relevant.factors)){
        temp <- list()
        factor <- which(relevant.factors)[f]
        num.omics <- length(SPEARobj$data$xlist)
        s <- 1
        for(o in 1:num.omics){
          omic.features <- list()
          w.loadings.current <- w.loadings[s:(s + ncol(SPEARobj$data$xlist[[o]]) - 1),factor]
          s <- s + ncol(SPEARobj$data$xlist[[o]])
          w.loadings.current <- w.loadings.current[which(w.loadings.current > cutoff)]
          w.loadings.current <- sort(w.loadings.current, decreasing = TRUE)
          if(length(w.loadings.current) > 0){
            temp[[names(SPEARobj$data$xlist)[o]]]$features <- names(w.loadings.current)
            temp[[names(SPEARobj$data$xlist)[o]]]$probabilities <- w.loadings.current
          }
        }
        factor.features[[names(factor)]] <- temp
      }
      feature.list[[colnames(SPEARobj$data$Y)[k]]] <- factor.features
    }
    else{
      feature.list[[colnames(SPEARobj$data$Y)[k]]] <- "No contributing factors found."
    }
  }
  return(feature.list)
}




# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# GET CV PREDICTIONS
# SPEARobj - a SPEAR object
# w - 'best' (choose best w for each response through cvm)
#     'overall' (choose best overall w for ALL responses through cvm)
#     'all' (return an extra dimension with all weights being shown)
SPEAR.get_cv_predictions <- function(SPEARobj, w = "best"){
  if(w == 'all'){
    w.idxs <- 1:length(SPEARobj$params$weights)
    results <- list()
    for(i in 1:length(w.idxs)){
      pred.mat <- SPEARobj$cv.eval$Yhat.keep[,,w.idxs[i]]
      if(is.null(dim(pred.mat))){
        pred.mat <- matrix(pred.mat, ncol = ncol(SPEARobj$data$Y))
      }
      colnames(pred.mat) <- colnames(SPEARobj$data$Y)
      rownames(pred.mat) <- rownames(SPEARobj$data$X)
      ws <- rep(SPEARobj$params$weights[w.idxs[i]], ncol(SPEARobj$data$Y))
      names(ws) <- colnames(SPEARobj$data$Y)
      temp <- list()
      temp$predictions <- pred.mat
      temp$weights <- ws
      results[[i]] <- temp
    }
  } else if(w == 'overall'){
      w.idx = which.min(apply(SPEARobj$cv.eval$cvm,1,sum))
      pred.mat <- SPEARobj$cv.eval$Yhat.keep[,,w.idx]
      colnames(pred.mat) <- colnames(SPEARobj$data$Y)
      rownames(pred.mat) <- rownames(SPEARobj$data$X)
      ws <- rep(SPEARobj$params$weights[w.idx], ncol(SPEARobj$data$Y))
      names(ws) <- colnames(SPEARobj$data$Y)
      temp <- list()
      temp$predictions <- pred.mat
      temp$weights <- ws
      results <- temp
  } else if(w == 'best'){
    w.idxs <- apply(SPEARobj$cv.eval$cvm,2,which.min)
    pred.mat <- matrix(0, nrow = nrow(SPEARobj$data$Y), ncol = ncol(SPEARobj$data$Y))
    for(i in 1:ncol(SPEARobj$data$Y)){
      pred.mat[,i] <- SPEARobj$cv.eval$Yhat.keep[,i,w.idxs[i]]
    }
    colnames(pred.mat) <- colnames(SPEARobj$data$Y)
    rownames(pred.mat) <- rownames(SPEARobj$data$X)
    ws <- SPEARobj$params$weights[w.idxs]
    names(ws) <- colnames(SPEARobj$data$Y)
    temp <- list()
    temp$predictions <- pred.mat
    temp$weights <- ws
    results <- temp
  }
  return(results)
}



# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# PLOT CV PREDICTIONS:
# SPEARobj - a SPEAR object
# ncol - number of columns to display in the final plot. Defaults to auto
# nrow - number of rows to display in the final plot. Defaults to 1
SPEAR.plot_cv_predictions <- function(SPEARobj, ncol = NULL, nrow = 1){
  plotlist <- list()
  w.idxs = apply(SPEARobj$cv.eval$cvm,2,which.min)
  for(k in 1:ncol(SPEARobj$data$Y)){
    x <- SPEARobj$cv.eval$Yhat.keep[,k,w.idxs[k]]
    y <- SPEARobj$data$Y[,k]
    fit <- lm(y ~ x)
    r2 <- summary(fit)$r.squared
    g <- ggplot() +
      geom_point(aes(x = SPEARobj$cv.eval$Yhat.keep[,k,w.idxs[k]], y = SPEARobj$data$Y[,k])) + 
      geom_smooth(aes(x = SPEARobj$cv.eval$Yhat.keep[,k,w.idxs[k]], y = SPEARobj$data$Y[,k]), color = "#74149E", method = "lm") +
      xlab("SPEAR Pred.") +
      ylab("Actual") +
      ggtitle(paste0(colnames(SPEARobj$data$Y)[k], " | w=", round(SPEARobj$params$weights[w.idxs[k]], 2), " | r2=", round(r2, 3))) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
    plotlist[[k]] <- g
  }
  pg <- cowplot::plot_grid(plotlist = plotlist, nrow = nrow, ncol = ncol)
  return(pg)
}



# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# GET FACTOR SCORES
# SPEARobj - a SPEAR object
# w - 'best' (choose best w for each response through cvm)
#     'overall' (choose best overall w for ALL responses through cvm)
#     'all' (return an extra dimension with all weights being shown)
SPEAR.get_factor_scores <- function(SPEARobj, w = "best"){
  if(w == 'all'){
    w.idxs <- 1:length(SPEARobj$params$weights)
    results <- list()
    for(i in 1:length(w.idxs)){
      Uhat <- SPEARobj$data$X %*% SPEARobj$fit$results$post_betas[,,,w.idxs[i]]
      colnames(Uhat) <- paste0("Factor", 1:SPEARobj$params$num_factors)
      rownames(Uhat) <- rownames(SPEARobj$data$X)
      temp <- list()
      temp$factor.scores <- Uhat
      temp$weight <- SPEARobj$params$weights[w.idxs[i]]
      results[[i]] <- temp
    }
  } 
  else if(w == 'overall'){
    w.idx = which.min(apply(SPEARobj$cv.eval$cvm,1,sum))
    results <- list()
    Uhat <- SPEARobj$data$X %*% SPEARobj$fit$results$post_betas[,,,w.idx]
    colnames(Uhat) <- paste0("Factor", 1:SPEARobj$params$num_factors)
    rownames(Uhat) <- rownames(SPEARobj$data$X)
    temp <- list()
    temp$factor.scores <- Uhat
    temp$weight <- SPEARobj$params$weights[w.idx]
    results <- temp
  } 
  else if(w == 'best'){
    w.idxs <- apply(SPEARobj$cv.eval$cvm,2,which.min)
    results <- list()
    for(i in 1:length(w.idxs)){
      Uhat <- SPEARobj$data$X %*% SPEARobj$fit$results$post_betas[,,,w.idxs[i]]
      colnames(Uhat) <- paste0("Factor", 1:SPEARobj$params$num_factors)
      rownames(Uhat) <- rownames(SPEARobj$data$X)
      temp <- list()
      temp$factor.scores <- Uhat
      temp$weight <- SPEARobj$params$weights[w.idxs[i]]
      results[[i]] <- temp
    }
    names(results) <- colnames(SPEARobj$data$Y)
  }
  return(results)
}




# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# PLOT FACTOR GRID
# SPEARobj - a SPEAR object
# w - 'best' (choose best w for each response through cvm)
#     'overall' (choose best overall w for ALL responses through cvm)
# groups - metadata, must be a named list of metadata where names(groups) = subject names
# factors - a list of factors to include. If NULL, will plot all factors
# return.as.list - return a list of grob objects (each one named with factors) (defaults to FALSE)
# include.legend - should a legend be put on the bottom of the plot (defaults to TRUE)
SPEAR.plot_factor_grid <- function(SPEARobj, w = "overall", groups = NULL, factors = NULL, return.as.list = FALSE, include.legend = TRUE){
  if(is.null(factors)){
    factors <- 1:SPEARobj$params$num_factors
  }
  if(w == 'overall'){
    w.idx = which.min(apply(SPEARobj$cv.eval$cvm,1,sum))
    results <- list()
    Uhat <- as.data.frame(SPEARobj$data$X %*% SPEARobj$fit$results$post_betas[,,,w.idx])
    Uhat <- Uhat[factors]
    colnames(Uhat) <- paste0("Factor", factors)
    rownames(Uhat) <- rownames(SPEARobj$data$X)
    Uhat$Group <- "NULL"
    if(!is.null(groups)){
      # Check for mapping of metadata to subjects:
      if(is.null(names(groups)) | any(!names(groups) %in% rownames(Uhat))){
        return("ERROR. 'groups' needs to be a named list to ensure the correct mapping to subjects. Please check that your subject names match.")
      }
      else{
        Uhat$Group <- sapply(rownames(Uhat), function(subject){return(groups[which(names(groups) == subject)])})
      }
    }
    # Generate factor vs. factor plots:
    plotlist <- list()
    for(j in 1:length(factors)){
      for(i in 1:length(factors)){
        if(i == j){
          g <- ggplot(data = Uhat) +
            theme_bw() +
            theme(plot.title = element_blank(), axis.title = element_blank())
          
          if(any(is.numeric(Uhat$Group))){
            g <- g + geom_density(aes_string(x = paste0("Factor",factors[i])), fill = "black", alpha = .15) +
              scale_color_distiller(palette = "RdBu", guide = FALSE)
          } else {
            g <- g + geom_density(aes_string(x = paste0("Factor",factors[i]), group = "Group", color = "Group", fill = "Group"), alpha = .15) +
              scale_color_brewer(palette = "Set1", guide = FALSE) +
              scale_fill_brewer(palette = "Set1", guide = FALSE)
          }
          
          plotlist[[length(plotlist) + 1]] <- g
        }
        else {
          g <- ggplot(data = Uhat) +
            geom_point(aes_string(x = paste0("Factor",factors[i]), y = paste0("Factor",factors[j]), group = "Group", color = "Group"), size = 1) + 
            #geom_smooth(aes(x = Uhat[,i], y = Uhat[,j]), color = "#74149E", method = "lm", formula = Uhat[,j] ~ Uhat[,i]) +
            theme_bw() +
            theme(plot.title = element_blank(), axis.title = element_blank())
          
          if(any(is.numeric(Uhat$Group))){
            g <- g + scale_color_distiller(palette = "RdBu", guide = FALSE)
          } else {
            g <- g + scale_color_brewer(palette = "Set1", guide = FALSE)
          }
          
          plotlist[[length(plotlist) + 1]] <- g
          names(plotlist)[length(plotlist)] <- paste0("X.Factor",factors[i],"_Y.Factor",factors[j])
        }
      }
      # Y Labels:
      g <- ggplot() +
        geom_text(aes(x = 0, y = 0), label = paste0("Factor ", factors[j]), angle = 270) + 
        ggtitle("") +
        theme_map() +
        theme(plot.title = element_blank(), axis.title.y = element_blank())
      
      plotlist[[length(plotlist) + 1]] <- g
      names(plotlist)[length(plotlist)] <- paste0("Y.Label.Factor", factors[j])
    }
    # X Labels:
    for(i in 1:length(factors)){
      g <- ggplot() +
        geom_text(aes(x = 0, y =99), label = paste0("Factor ", factors[i]), angle = 0) + 
        ggtitle("") +
        theme_map() +
        theme(plot.title = element_blank(), axis.title.y = element_blank())
      
      plotlist[[length(plotlist) + 1]] <- g
      names(plotlist)[length(plotlist)] <- paste0("X.Label.Factor", factors[i])
    }
    # Legend:
    if(include.legend){
      g <- ggplot(data = Uhat) +
        geom_point(aes_string(x = paste0("Factor",factors[i]), y = paste0("Factor",factors[j]), group = "Group", color = "Group")) +
        scale_fill_brewer(palette = "Set1", guide = FALSE) +
        theme_bw() +
        theme(legend.key.size = unit(1, 'cm'), legend.position = "top")
      if(any(is.numeric(Uhat$Group))){
        g <- g + scale_color_distiller(palette = "RdBu") +
          guides(color = guide_colourbar(barheight = 0.5))
      } else {
        g <- g + scale_color_brewer(palette = "Set1")
      }
      legend <- get_legend(g)
    }
  } 
  else if(w == 'best'){ w.idx
    w.idxs <- apply(SPEARobj$cv.eval$cvm,2,which.min)
    results <- list()
    for(i in 1:length(w.idxs)){
      Uhat <- SPEARobj$data$X %*% SPEARobj$fit$results$post_betas[,,,w.idxs[k]]
      colnames(Uhat) <- paste0("Factor", 1:SPEARobj$params$num_factors)
      rownames(Uhat) <- rownames(SPEARobj$data$X)
      temp <- list()
      temp$factor.scores <- Uhat
      temp$weight <- SPEARobj$params$weights[w.idxs[i]]
      results[[i]] <- temp
    }
    names(results) <- colnames(SPEARobj$data$Y)
  }
  if(return.as.list){
    return(plotlist)
  }
  else{
    title <- ggdraw() + 
      draw_label(
        paste0("SPEAR Factor Scores for w=", round(SPEARobj$params$weights[w.idx], 2)),
        fontface = 'bold',
        x = 0,
        hjust = 0
      ) +
      theme(
        # add margin on the left of the drawing canvas,
        # so title is aligned with left edge of first plot
        
        plot.margin = margin(0, 0, 0, 15)
      )

    p <- cowplot::plot_grid(plotlist = plotlist, ncol = length(factors) + 1, rel_heights = c(rep(1, length(factors)), .3), rel_widths = c(rep(1, length(factors)), .3))
    if(include.legend){
      p <- cowplot::plot_grid(title, p, legend, ncol = 1, rel_heights = c(1, 10, 1))
    } else {
      p <- cowplot::plot_grid(title, p, ncol = 1, rel_heights = c(3, 20))
    }
    return(p)
  }
}



# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# PLOT FACTOR SCORES
# SPEARobj - a SPEAR object
# w - 'best' (choose best w for each response through cvm)
#     'overall' (choose best overall w for ALL responses through cvm)
# groups - metadata, must be a named list of metadata where names(groups) = subject names
# factors - a list of factors to include. If NULL, will plot all factors
# return.as.list - return a list of grob objects (each one named with factors) (defaults to FALSE)
# include.legend - should a legend be put on the bottom of the plot (defaults to TRUE)
SPEAR.plot_factor_scores <- function(SPEARobj, w = "overall", groups = NULL, factors = NULL, return.as.list = FALSE, include.legend = TRUE){
  if(is.null(factors)){
    factors <- 1:SPEARobj$params$num_factors
  }
  if(w == 'overall'){
    w.idx = which.min(apply(SPEARobj$cv.eval$cvm,1,sum))
    results <- list()
    Uhat <- as.data.frame(SPEARobj$data$X %*% SPEARobj$fit$results$post_betas[,,,w.idx])
    Uhat <- Uhat[factors]
    colnames(Uhat) <- paste0("Factor", factors)
    rownames(Uhat) <- rownames(SPEARobj$data$X)
    Uhat$Subject <- rownames(Uhat)
    Uhat$Group <- ""
    if(!is.null(groups)){
      # Check for mapping of metadata to subjects:
      if(is.null(names(groups)) | any(!names(groups) %in% rownames(Uhat))){
        return("ERROR. 'groups' needs to be a named list to ensure the correct mapping to subjects. Please check that your subject names match.")
      }
      else{
        Uhat$Group <- sapply(rownames(Uhat), function(subject){return(groups[which(names(groups) == subject)])})
      }
    }
    
    if(class(Uhat$Group) == "factor" | class(Uhat$Group) == "character"){
      group.type <- "categorical"
    } else {
      group.type <- "numerical"
    }
    
    # Melt dataframe:
    Uhat.melt <- melt(Uhat, id.vars = c("Subject", "Group"))
    
    # Generate factor scores plot:
    if(group.type == "categorical"){
      # categorical, plot with boxplots
      g <- ggplot(filter(Uhat.melt, variable %in% paste0("Factor", factors))) +
        geom_jitter(aes(x = Group, y = value, color = Group), pch = 20) +
        geom_boxplot(aes(x = Group, y = value), alpha = .5, outlier.alpha = 0) +
        scale_color_brewer(palette = "RdBu", direction = -1) +
        xlab(NULL) +
        ylab("Factor Score") +
        theme_bw() +
        facet_wrap(vars(variable), nrow = 1)
    } else {
      # Numeric, plot with regression
      g <- ggplot(filter(Uhat.melt, variable %in% paste0("Factor", factors))) +
        geom_point(aes(x = Group, y = value, fill = Group), pch = 21, lwd = 1.5) +
        geom_smooth(aes(x = Group, y = value), color = "black", lwd = .6, alpha = .5, method = "lm") +
        scale_fill_distiller(palette = "RdBu") +
        xlab(NULL) +
        ylab("Factor Score") +
        theme_bw() +
        facet_wrap(vars(variable), nrow = 1)
    }

    return(g)
  }
}
# Test numeric:
#g <- as.vector(SPEARobj$data$Y)
#names(g) <- rownames(SPEARobj$data$Y)
# Test categorical:
#g <- factor(as.character(floor(as.vector(SPEARobj$data$Y))), levels = c("-3", "-2", "-1", "0", "1", "2"))
#names(g) <- rownames(SPEARobj$data$Y)
#SPEAR.plot_factor_scores(SPEARobj, groups = g)




# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# GET PATHWAY INFO FROM GENES:
# example

#feature.list <- SPEAR.get_factor_features(SPEARobj)
#gene.list <- feature.list[['FHA']]$Factor2$X3$features
## Modify the features to be regular Ensembl identifiers:
#gene.list <- unname(sapply(gene.list, function(ensemblid){
#  return(stringr::str_split(ensemblid, pattern = '\\.')[[1]][1])
#}))
#gostres <- gprofiler2::gost(query = gene.list, organism = "hsapiens")
#res_table <- as_tibble(gostres$result) %>% select(source, term_id, term_name, term_size, intersection_size, p_value) %>% arrange(p_value)
#res_table

#s <- run_cv_spear(data.total$data.tr$xlist, data.total$data.tr$Y, seed = 123)


