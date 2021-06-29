#' Update the dimension names for all SPEARobject matrices. Used internally.
#'@export
update.dimnames = function(){
  dimnames(self$fit$regression.coefs) = list(colnames(self$data$train$X), paste0("Factor", 1:self$params$num_factors), paste0("w.idx.", 1:nrow(self$params$weights)))
  dimnames(self$fit$projection.coefs.x) = list(colnames(self$data$train$X), paste0("Factor", 1:self$params$num_factors), paste0("w.idx.", 1:nrow(self$params$weights))) 
  dimnames(self$fit$projection.coefs.y) = list(paste0("Factor", 1:self$params$num_factors), colnames(self$data$train$Y), paste0("w.idx.", 1:nrow(self$params$weights))) 
  dimnames(self$fit$nonzero.probs) = list(colnames(self$data$train$X), paste0("Factor", 1:self$params$num_factors), paste0("w.idx.", 1:nrow(self$params$weights))) 
  dimnames(self$fit$projection.probs) = list(colnames(self$data$train$X), paste0("Factor", 1:self$params$num_factors), paste0("w.idx.", 1:nrow(self$params$weights))) 
  dimnames(self$fit$marginal.probs) = list(colnames(self$data$train$X), paste0("Factor", 1:self$params$num_factors), paste0("w.idx.", 1:nrow(self$params$weights))) 
  dimnames(self$fit$joint.probs) = list(colnames(self$data$train$X), paste0("Factor", 1:self$params$num_factors), paste0("w.idx.", 1:nrow(self$params$weights)))
  if(self$fit$type == "cv"){
    dimnames(self$fit$regression.coefs.cv) = list(colnames(self$data$train$X), paste0("Factor", 1:self$params$num_factors), paste0("Fold", 1:self$params$num.folds), paste0("w.idx.", 1:nrow(self$params$weights)))
    dimnames(self$fit$projection.coefs.y.cv) = list(paste0("Factor", 1:self$params$num_factors), colnames(self$data$train$Y), paste0("Fold", 1:self$params$num.folds), paste0("w.idx.", 1:nrow(self$params$weights)))
  }
}



#' Add a dataset to a SPEARobject. Needs to have the same number of datasets/columns as `SPEARobject$data$train$Xlist`. Response not required (will be recorded as NULL). Can be accessed via `SPEARobject$data$__name__`.
#' @param X List of explanatory dataset matrices. Needs to have same column names and length (number of matrices) as SPEARobject$data$Xlist.
#' @param Y Response matrix. Needs to have the same number of rows/length as the number of rows in parameter `X`. Defaults to `NULL` (not required).
#' @param Z Full matrix with all data. Cannot have missing values; must be imputed. If left blank (`NULL`) will be equivalent to `do.call("cbind", X)` after imputation.
#' @param name Name for the stored dataset in the SPEARobject. Access via `SPEARobject$data$__name__`.
#' @examples
#' SPEARobj <- make.SPEARobject(...)
#' 
#' X.test = ... # Define X.test as a list of matrices that have the same columns as SPEARobj$data$train$Xlist
#' SPEARobj$add.data(X = X.test) # will default name to "dataset1"
#' 
#' X.test2 = ...
#' Y.test2 = ...
#' SPEARobj$add.data(X = X.test2, Y = Y.test2, name = "test2")
#'@export
add.data = function(X, Y = NULL, Z = NULL, name = NULL){
  
  # Set name:
  if(is.null(name)){
    d.tmp = 1
    name = paste0("dataset", d.tmp)
    while(name %in% names(self$data)){
      d.tmp = d.tmp + 1
      name = paste0("dataset", d.tmp)
    }
  } else if(name %in% names(self$data)){
    cat("WARNING: name provided already exists in names(SPEARobject$data). Overwriting...\n")
  }
  
  if(name != "train"){
    # X...
    if(any(sapply(X, function(X.d){return(class(X.d)[1])}) != "matrix")){
      X <- lapply(X, as.matrix)
    }
    # Check that length is the same
    if(length(X) != length(self$data$train$Xlist)){
      stop("ERROR: Number of datasets provided (", length(X), ") is not equal to the number of datasets used to train the SPEARobject (", length(self$data$train$Xlist), "). Needs to be equal.")
    }
    # Check that column dimensions are the same
    if(any(sapply(1:length(X), function(d){return(ncol(X[[d]]) != ncol(self$data$train$Xlist[[d]]))}))){
      stop("ERROR: One or more datasets provided in X has a different number of columns than those in SPEARobject$data$train$Xlist. Needs to match.")
    }
  }
  # Process:
  if(is.null(names(X))){
    names(X) <- paste("dataset", 1:length(X))
  }
  if(all(sapply(X, function(X.d){return(is.null(rownames(X.d)))}))){
    for(d in 1:length(X)){
      rownames(X[[d]]) <- paste0("sample_", 1:nrow(X[[d]]))
    }
  }
  for(d in 1:length(X)){
    if(!is.null(rownames(X[[d]]))){
      rownames.temp <- rownames(X[[d]])
      break
    }
  }
  if(any(sapply(X, function(X.d){return(is.null(rownames(X.d)))}))){
    for(d in 1:length(X)){
      rownames(X[[d]]) <- rownames.temp
    }
  } else if(any(sapply(X, function(X.d){return(rownames(X.d) != rownames.temp)}))){
    for(d in 1:length(X)){
      rownames(X[[d]]) <- rownames.temp
    }
  }
  for(d in 1:length(X)){
    if(is.null(colnames(X[[d]]))){
      colnames(X[[d]]) <- paste0(names(X)[d], "_feat", 1:ncol(X[[d]]))
    }
  }
  if(any(table(sapply(1:length(X), function(d){return(colnames(d))})) > 1)){
    for(d in 1:length(X)){
      colnames(X[[d]]) <- paste0(names(X)[d], "_", colnames(X[[d]]))
    }
  }
  X_ <- do.call("cbind", X)
  # Check that columns all match columns
  if(name != "train"){
    if(any(colnames(X_) != colnames(self$data$train$X))){
      stop("ERROR: One or more column names in X provided were not found in the data used to train this SPEARobject (SPEARobject$data$train$X). All features must match to add a dataset.")
    }
  }

  # Xobs
  Xobs <- array(1, dim  = dim(X_))
  colnames(Xobs) <- colnames(X_)
  rownames(Xobs) <- rownames(X_)
  if(any(is.na(X_))){
    Xobs[which(is.na(X_))] <- 0
  } else if(any(is.nan(X_))){
    Xobs[which(is.nan(X_))] <- 0
  }
  
  # Y...
  if(!is.null(Y)){
    # Check if Y has no dimensions:
    if(is.null(dim(Y))){
      Y = matrix(Y, ncol = 1)
      colnames(Y) = "Y"
      rownames(Y) = rownames(X)
    } else if(class(Y)[1] != "matrix"){
      Y = as.matrix(Y)
    } else{
      Y = Y
    }
    if(is.null(colnames(Y))){colnames(Y) = paste0("Y", ncol(Y))}
    if(is.null(rownames(Y))){rownames(Y) = rownames(X)}
    if(any(rownames(Y) != rownames(X))){stop("ERROR: rownames for supplied Y do not match rownames for supplied X. Need to match exactly.")}
    Yobs <- array(1, dim  = dim(Y))
    colnames(Yobs) <- colnames(Y)
    rownames(Yobs) <- rownames(Y)
    if(any(is.na(Y))){
      Yobs[which(is.na(Y))] <- 0
    } else if(any(is.nan(Y))){
      Yobs[which(is.nan(Y))] <- 0
    }
  } else {
    Yobs = NULL
  }
  
  if(is.null(Z)){
    if(sum(Xobs) == 0){
      Z = X_
    } else {
      # Need to impute...
      Z = private$impute.z(X_)
    }
  } else {
    # Check if Z has same names, dimensions as X...
    if(nrow(Z) != nrow(X_) | ncol(Z) != ncol(X_)){
      stop("ERROR: Z provided has different dimensions than the provided X when combined (cbind). Needs to be the same.")
    }
    if(any(colnames(Z) != colnames(X_))){
      stop("ERROR: Z provided has different column names than the provided X. If duplicate column names in X were found, they may have been renamed.")
    }
    if(any(is.na(Z)) | any(is.nan(Z))){
      Z = private$impute.z(X_)
    } 
  }
  
  # name and store:
  self$data[[name]] <- list(
    Xlist = X,
    X = X_,
    Xobs = Xobs,
    Y = Y,
    Yobs = Yobs,
    Z = Z
  )
  if(!self$options$quiet & name != "train"){cat("Saved dataset to $data$", name, "\n\n", sep = "")}
  return(invisible(self))
}



#' Impute the rest of a missing matrix (Z). Called internally.
#' @param Z
#' @param method Which imputation method?
#'@export
impute.z = function(Z, method = NULL){
  # TODO...
  #
  #
  return(Z)
}



#' Confirm fold.ids will work in cv.spear 
#' @param fold.ids
#'@export
check.fold.ids = function(fold.ids){
  return(any(sapply(unique(fold.ids), function(k){
    for(j in 1:ncol(self$data$train$Y)){
      subset = self$data$train$Y[which(fold.ids != k),j]
      if(any(table(subset) < 2)){
        return(TRUE)
      }
    }
    return(FALSE)
  })))
}



#' Choose a combination of weights for a SPEARobject. Print `SPEARobject$params$weights` to see all combinations trained.
#' @param w.idx Use a manual weight index (row from `SPEARobject$params$weights`)? Defaults to `NULL`
#' @param method Which method to 
#'@export
set.weights = function(w.idx = NULL, method = NULL){
  # TODO... add other ways
  #
  #
  if(!is.null(w.idx)){
    if(!w.idx %in% 1:nrow(self$params$weights)){
      stop("ERROR: w.idx supplied not valid (must be between 1 and ", nrow(self$params$weights), "). Type SPEARobject$params$weights to see all combinations (rows).")
    } else {
      if(!self$options$quiet){
        cat("Setting current weight index to ", w.idx, "\n", 
            private$color.text("w.x: ", "green"), self$params$weights[w.idx,1], "\n",
            private$color.text("w.y: ", "green"), self$params$weights[w.idx,2], "\n",
            sep = "")
      }
      self$options$current.weight.idx = w.idx
    }
    
  } else if(!is.null(method)){
    cat("To do...")
    
    
    
    
    
    
  }
  if(!self$options$quiet){cat("\n")}
  return(invisible(self))
}



#' Get factor scores from a SPEARobject for a particular dataset `data`.
#' @param data Which dataset to use? Can be any dataset listed under `$data$____`. Defaults to `"train"`.
#' @param cv If `data = "train"`, get factor scores generated from `$run.cv.spear`? If `$run.spear` was used or if `data != "train"` this parameter is ignored. Defaults to `FALSE`. NOTE: There is an element of randomness if the factor scores are not correlated with the response, so it is recomended to view the factor scores with `cv = FALSE`
#'@export
get.factor.scores = function(data = "train", cv = FALSE){
  private$check.fit()
  if(data == "train" & self$fit$type == "cv" & cv){
    factor.scores.cv <- matrix(NA, nrow = nrow(self$data[[data]]$X), ncol = self$params$num_factors)
    for(i in 1:length(self$params$fold.ids)){
      factor.scores.cv[i,] <- self$data[[data]]$X[i,] %*% self$fit$regression.coefs.cv[,,self$params$fold.ids[i],self$options$current.weight.idx]
    }
    rownames(factor.scores.cv) <- rownames(self$data[[data]]$X)
    colnames(factor.scores.cv) <- paste0("Factor", 1:ncol(factor.scores.cv))
    return(factor.scores.cv)
  }
  factor.scores <- self$data[[data]]$X %*% self$fit$regression.coefs[,,self$options$current.weight.idx]
  rownames(factor.scores) <- rownames(self$data[[data]]$X)
  colnames(factor.scores) <- paste0("Factor", 1:ncol(factor.scores))
  return(factor.scores)
}



#' Get predictions from a SPEARobject and a dataset `data`.
#' @param data Which dataset to use? Can be any dataset listed under `$data$____`. Defaults to `"train"`.
#' @param cv If `data = "train"`, get factor scores generated from `$run.cv.spear`? If `$run.spear` was used or if `data != "train"` this parameter is ignored. Defaults to `TRUE`.
#'@export
get.predictions = function(data = "train", cv = TRUE){
  # Three ways to get the same signals:
  #yhat_te = SPEARobj$data$train$X %*% SPEARobj$fit$cv.eval$reg_coefs[,,6]
  #yhat_te = SPEARobj$data$train$X %*% SPEARobj$fit$regression.coefs[,,6] %*% SPEARobj$fit$cv.eval$projection_coefs_scaled[,,6]
  #yhat_te = SPEARobj$set.weights(6)$get.factor.scores(cv = FALSE) %*% SPEARobj$fit$projection.coefs.y.scaled[,,6]
  # signal from CV
  #yhat_te = matrix(NA, nrow=379, ncol = 4)
  #fs <- SPEARobj$get.factor.scores(cv = TRUE)
  #for(i in 1:nrow(yhat_te)){
  #  yhat_te[i,] <- fs[i,] %*% SPEARobj$fit$projection.coefs.y.cv.scaled[,,SPEARobj$params$fold.ids[i],6]
  #}
  # Get preds from signal:
  #intercept = sapply(SPEARobj$fit$cv.eval$intercepts, function(z) z[6])
  #prob_hat_te = t(apply(yhat_te, 1, function(z) exp(z+intercept)/sum(exp(z+intercept))))
  #n.preds = apply(prob_hat_te, 1, which.max)
  
  private$check.fit()
  if(data == "train" & self$fit$type == "cv" & cv){
    factor.scores <- self$get.factor.scores(data = data, cv = TRUE)
    preds = matrix(NA, nrow = nrow(self$data[[data]]$X), ncol = ncol(self$data$train$Y))
    for(i in 1:nrow(preds)){
      preds[i,] <- factor.scores[i,] %*% self$fit$projection.coefs.y.cv.scaled[,,self$params$fold.ids[i],self$options$current.weight.idx]
    }
    rownames(preds) <- rownames(self$data[[data]]$X)
    colnames(preds) <- colnames(self$data$train$Y)
  } else {
    factor.scores <- self$get.factor.scores(data = data, cv = FALSE)
    preds = matrix(NA, nrow = nrow(self$data[[data]]$X), ncol = ncol(self$data$train$Y))
    for(i in 1:nrow(preds)){
      preds[i,] <- factor.scores[i,] %*% self$fit$projection.coefs.y.scaled[,,self$options$current.weight.idx]
    }
    rownames(preds) <- rownames(self$data[[data]]$X)
    colnames(preds) <- colnames(self$data$train$Y)
  }
  if(self$params$family == "gaussian"){
    return(preds)
  } else if(self$params$family == "multinomial"){
    signal <- preds
    intercept = sapply(self$fit$intercepts.scaled, function(z) z[self$options$current.weight.idx])
    Pmat = t(apply(signal, 1, function(z) exp(z+intercept)/sum(exp(z+intercept))))
    predictions <- matrix(apply(Pmat, 1, which.max) - 1, ncol=1)
    predictions.encoded = matrix(0, ncol = ncol(Pmat), nrow = nrow(Pmat))
    for(i in 1:nrow(predictions.encoded)){
      predictions.encoded[i,predictions[i]] <- 1
    }
    rownames(predictions) <- rownames(preds)
    colnames(predictions) <- "pred.class"
    results.combined <- list(predictions = predictions,
                             predictions.encoded = predictions.encoded,
                             probabilities = Pmat,
                             signal = signal)
    return(results.combined)
  }
  
  # If family == "binomial" or "ordinal"
  results = list()
  for(j in 1:ncol(preds)){
    resp.name <- colnames(preds)[j]
    signal <- matrix(preds[,j], ncol = 1)
    rownames(signal) <- rownames(preds)
    colnames(signal) <- c(resp.name)
    intercept = self$fit$intercepts.y[[j]]
    Pmat0 = matrix(0, ncol = length(intercept), nrow = length(signal))
    Pmat = matrix(0, ncol = length(intercept) + 1, nrow = length(signal))
    for(k in 1:ncol(Pmat0)){
      Pmat0[,k] = 1.0/(1+exp(-signal - intercept[k]))
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
    class.predictions <- matrix(apply(Pmat, 1, which.max) - 1, ncol=1)
    rownames(class.predictions) <- rownames(preds)
    colnames(class.predictions) <- c(resp.name)
    rownames(Pmat) <- rownames(preds)
    colnames(Pmat) <- paste0("Class", 0:(ncol(Pmat)-1))
    results[[resp.name]] <- list(class.predictions = class.predictions,
                                       class.probabilities = Pmat,
                                       signal = signal)
  }
  return(results)
}



#' Get features from a SPEARobject.
#' @param Factor
#'@export
get.features = function(rank = "probability", factors = NULL, omics = NULL, coefficient.cutoff = 0, probability.cutoff = .95, correlation = "none"){
  private$check.fit()
  w.idx <- self$options$current.weight.idx
  if(is.null(factors)){
    factors <- 1:self$params$num_factors
  }
  factors.string <- paste0("Factor", factors)
  if(is.null(omics)){
    omics <- names(self$data$train$Xlist)
  }
  features <- list()
  w.joint.probabilities = self$fit$joint.probs[,,w.idx]
  w.marginal.probabilities = self$fit$marginal.probs[,,w.idx]
  w.projection.probabilities = self$fit$projection.probs[,,w.idx]
  w.regression.coefficients = self$fit$regression.coefs[,,w.idx]
  w.projection.coefficients = self$fit$projection.coefs.x[,,w.idx]
  colnames(w.joint.probabilities) <- paste0("Factor", 1:ncol(w.joint.probabilities))
  rownames(w.joint.probabilities) <- colnames(self$data$train$X)
  colnames(w.marginal.probabilities) <- paste0("Factor", 1:ncol(w.marginal.probabilities))
  rownames(w.marginal.probabilities) <- colnames(self$data$train$X)
  colnames(w.projection.probabilities) <- paste0("Factor", 1:ncol(w.projection.probabilities))
  rownames(w.projection.probabilities) <- colnames(self$data$train$X)
  colnames(w.regression.coefficients) <- paste0("Factor", 1:ncol(w.joint.probabilities))
  rownames(w.regression.coefficients) <- colnames(self$data$train$X)
  colnames(w.projection.coefficients) <- paste0("Factor", 1:ncol(w.joint.probabilities))
  rownames(w.projection.coefficients) <- colnames(self$data$train$X)
  # by factor
  num.omics <- length(self$data$train$Xlist)
  factor.features <- list()
  for(f in 1:self$params$num_factors){
    temp <- list()
    factor <- paste0("Factor", f)
    s <- 1
    for(o in 1:num.omics){
      omic.features <- list()
      w.joint.probabilities.current <- w.joint.probabilities[s:(s + ncol(self$data$train$Xlist[[o]]) - 1),factor]
      w.marginal.probabilities.current <- w.marginal.probabilities[s:(s + ncol(self$data$train$Xlist[[o]]) - 1),factor]
      w.projection.probabilities.current <- w.projection.probabilities[s:(s + ncol(self$data$train$Xlist[[o]]) - 1),factor]
      w.regression.coefficients.current <- w.regression.coefficients[s:(s + ncol(self$data$train$Xlist[[o]]) - 1),factor]
      w.projection.coefficients.current <- w.projection.coefficients[s:(s + ncol(self$data$train$Xlist[[o]]) - 1),factor]
      s <- s + ncol(self$data$train$Xlist[[o]])
      # Round loading probabilities to 7 decimals (to rank 1's):
      #w.joint.probabilities.current <- round(w.joint.probabilities.current, 7)
      #w.marginal.probabilities.current <- round(w.marginal.probabilities.current, 7)
      #w.projection.probabilities.current <- round(w.projection.probabilities.current, 7)
      feat.order <- order(w.joint.probabilities.current, abs(w.projection.coefficients.current), decreasing = TRUE)
      w.regression.coefficients.current <- w.regression.coefficients.current[feat.order]
      w.projection.coefficients.current <- w.projection.coefficients.current[feat.order]
      w.joint.probabilities.current <- w.joint.probabilities.current[feat.order]
      w.marginal.probabilities.current <- w.marginal.probabilities.current[feat.order]
      w.projection.probabilities.current <- w.projection.probabilities.current[feat.order]
      temp[[names(self$data$train$Xlist)[o]]]$features <- names(w.joint.probabilities.current)
      temp[[names(self$data$train$Xlist)[o]]]$joint.probabilities <- w.joint.probabilities.current
      temp[[names(self$data$train$Xlist)[o]]]$marginal.probabilities <- w.marginal.probabilities.current
      temp[[names(self$data$train$Xlist)[o]]]$projection.probabilities <- w.projection.probabilities.current
      temp[[names(self$data$train$Xlist)[o]]]$projection.coefficients <- w.projection.coefficients.current
      temp[[names(self$data$train$Xlist)[o]]]$regression.coefficients <- w.regression.coefficients.current
    }
    features[[factor]] <- temp
  }
  fs <- self$get.factor.scores()
  fs.cv <- self$get.factor.scores(cv = TRUE)
  results.list <- list()
  for(f in factors.string){
    for(o in omics){
      feature.vec <- features[[f]][[o]]
      passed.coeff.cutoff <- abs(feature.vec$projection.coefficients) >= coefficient.cutoff
      passed.prob.cutoff <- feature.vec$joint.probabilities >= probability.cutoff
      passed.total.cutoff <- passed.coeff.cutoff & passed.prob.cutoff
      if(sum(passed.total.cutoff) > 0){
        res <- as.data.frame(feature.vec)[passed.coeff.cutoff & passed.prob.cutoff,]
        res$omic <- o
        res$factor <- f
        res$cor.in.sample <- NA
        res$cor.cv <- NA
        res$pval.in.sample <- NA
        res$pval.cv <- NA
        # Correlations:
        if(correlation != "none"){
          for(r in 1:nrow(res)){
            vals <- self$data$train$Xlist[[o]][,which(colnames(self$data$train$Xlist[[o]]) == res$features[r])]
            res$cor.in.sample[r] <- cor(vals, fs[,f], method = correlation)
            res$cor.cv[r] <- cor(vals, fs.cv[,f], method = correlation)
            res$pval.in.sample[r] <- suppressWarnings(cor.test(vals, fs[,f], method = correlation)$p.value)
            res$pval.cv[r] <- suppressWarnings(cor.test(vals, fs.cv[,f], method = correlation)$p.value)
          }
          colnames(res) <- c("Feature", "joint.probability", "marginal.probability", "projection.probability", "projection.coefficient", "regression.coefficient", "Dataset", "Factor", "in.sample.correlation", "cv.correlation", "in.sample.pvalue", "cv.pvalue")
          rownames(res) <- NULL
          res <- dplyr::select(res, 
                               Factor,
                               Feature, 
                               Dataset, 
                               regression.coefficient,
                               projection.coefficient, 
                               joint.probability,
                               marginal.probability,
                               projection.probability,
                               in.sample.correlation, 
                               in.sample.pvalue, 
                               cv.correlation, 
                               cv.pvalue)
        } else {
          colnames(res) <- c("Feature", "joint.probability", "marginal.probability", "projection.probability", "projection.coefficient", "regression.coefficient", "Dataset", "Factor")
          rownames(res) <- NULL
          res <- dplyr::select(res, 
                               Factor,
                               Feature, 
                               Dataset, 
                               regression.coefficient,
                               projection.coefficient, 
                               joint.probability,
                               marginal.probability,
                               projection.probability)
        }
        results.list[[length(results.list) + 1]] <- res
      }
    }
  }
  if(length(results.list) == 0){
    warning("*** Warning: No features found within cutoffs for chosen Factors/Omics. Returning data.frame with 0 rows. Try broadening the criteria.")
    total.results <- data.frame( Factor = character(0),
                                 Feature = character(0), 
                                 Dataset = character(0), 
                                 regression.coefficient = numeric(0),
                                 projection.coefficient = numeric(0), 
                                 joint.probability = numeric(0),
                                 marginal.probability = numeric(0),
                                 projection.probability = numeric(0))
  } else {
    total.results <- do.call("rbind", results.list)
    total.results <- dplyr::arrange(total.results, Factor, -joint.probability, -abs(projection.coefficient))
  }
  
  return(total.results)
}


#' Get factor contributions for a SPEARobject.
#' @param do.X Calculate contributions for X (train)? Defaults to `TRUE`
#' @param do.Y Calculate contributions for Y? Defaults to `TRUE`
#' @param do.Y.pvals Include p.values for y? Defaults to `FALSE`
#'@export
get.contributions = function(do.X = TRUE, do.Y = TRUE, do.Y.pvals = FALSE){
  output <- list()
  w.idx <- self$options$current.weight.idx
  if(do.X){
    # X:
    W.to_X <- self$fit$projection.coefs.x[,,w.idx]
    Uhat <- self$get.factor.scores(data = "train", cv = FALSE)
    W.to_X.list <- list()
    ind <- 1
    for(d in 1:length(self$data$train$Xlist)){
      W.to_X.list[[d]] <- W.to_X[ind:(ind - 1 + ncol(self$data$train$Xlist[[d]])),]
      ind <- ind + ncol(self$data$train$Xlist[[d]])
    }
    factor_contributions = array(NA,dim = c(self$params$num_factors, length(self$data$train$Xlist)))
    factor_contributions_pvals = array(NA,dim = c(self$params$num_factors, length(self$data$train$Xlist)))
    for(k in 1:self$params$num_factors){
      for(o in 1:length(self$data$train$Xlist)){
        factor.contributions <- c()
        X.tilde <- Uhat[,k] %*% t(W.to_X.list[[o]][,k])
        # Get var.explained per feature:
        for(j in 1:ncol(X.tilde)){
          X.tilde.j <- X.tilde[,j]
          X.j <- self$data$train$Xlist[[o]][,j]
          X.tilde.j = X.tilde.j + rnorm(n = length(X.tilde.j), mean = 0, sd = sqrt(var(X.j))*1/nrow(self$data$train$X))
          suppressWarnings( tmp_pearson <- cor.test(X.j, X.tilde.j, method = 'pearson') ) 
          if(is.na(tmp_pearson$estimate)){
            factor.contributions <- c(factor.contributions, 0)
          }else if( tmp_pearson$estimate<0){
            factor.contributions <- c(factor.contributions, 0)
          }else{
            factor.contributions <- c(factor.contributions, tmp_pearson$estimate^2)
          }
        }
        # Return mean as factor_contributions:
        factor_contributions[k,o] <- mean(factor.contributions)
      }
    }
    colnames(factor_contributions) <- names(self$data$train$Xlist)
    rownames(factor_contributions) <- paste0("Factor", 1:nrow(factor_contributions))
    output$X <- factor_contributions
  }
  # Y:
  if(do.Y){
    var_explained_Y <- matrix(NA, nrow = dim(self$fit$cv.eval$factor_contributions)[1], ncol = dim(self$fit$cv.eval$factor_contributions)[2])
    var_explained_Y.pvals <- matrix(NA, nrow = dim(self$fit$cv.eval$factor_contributions)[1], ncol = dim(self$fit$cv.eval$factor_contributions)[2])
    for(j in 1:dim(self$fit$cv.eval$factor_contributions)[2]){
      if(length(self$fit$cv.eval$factor_contributions[,j,w.idx]) != 0){
        var_explained_Y[,j] <- self$fit$cv.eval$factor_contributions[,j,w.idx]
      }
      if(length(self$fit$cv.eval$factor_contributions_pvals[,j,w.idx]) != 0){
        var_explained_Y.pvals[,j] <- self$fit$cv.eval$factor_contributions_pvals[,j,w.idx]
      }
    }
    colnames(var_explained_Y) <- colnames(self$data$train$Y)
    rownames(var_explained_Y) <- paste0("Factor", 1:nrow(var_explained_Y))
    colnames(var_explained_Y.pvals) <- colnames(self$data$train$Y)
    rownames(var_explained_Y.pvals) <- paste0("Factor", 1:nrow(var_explained_Y.pvals))
    output$Y <- var_explained_Y
    if(do.Y.pvals){
      output$Y.pvals <-var_explained_Y.pvals
    }
  }
  return(output)
}



#' Check if a fit object is generated and if a weight has been set...
#'@export
check.fit = function(){
  # Check if a fit object exists:
  if(is.null(self$fit)){
    stop("ERROR: $fit is NULL. Run $run.cv.spear(...) or $run.spear() to generate the $fit object.")
  }
  if(is.null(self$options$current.weight.idx)){
    stop("ERROR: $options$current.weight.idx is NULL. Set this parameter manually [1 through nrow($params$weights)] or use $set.weights()")
  }
}



#' Print out a variety of summary information about a SPEARobject
#' @param type Which type of summary to print? Can be "data". Defaults to NULL.
#' @param remove.formatting Remove text color/bold font face. Defaults to FALSE.
#' @param quiet Do not print anything. Defaults to FALSE.
#' @examples
#' SPEARobj <- make.SPEARobject(...)
#' 
#' SPEARobj$run.spear()
#' 
#'@export
print.out = function(type = NULL, remove.formatting = FALSE, quiet = FALSE){
  if(quiet){
    return(NULL)
  }
  success.color <- "light green"
  update.color <- "green"
  error.color.bold <- "light red"
  error.color <- "red"
  dataset.color <- "light cyan"
  response.color <-"yellow"
  
  if(is.null(type)){
    cat("")
  }
  # ---------------------------------
  if(type == "data"){
    cat(paste0("Detected ", length(self$data$train$Xlist), " datasets:\n"))
    for(i in 1:length(self$data$train$Xlist)){
      cat(private$color.text(names(self$data$train$Xlist)[i], dataset.color), "\tSubjects: ", nrow(self$data$train$Xlist[[i]]), "\tFeatures: ", ncol(self$data$train$Xlist[[i]]), "\n")
    }
    cat(paste0("Detected ", ncol(self$data$train$Y), " response ", ifelse(ncol(self$data$train$Y) == 1, "variable", "variables"), ":\n"))
    for(i in 1:ncol(self$data$train$Y)){
      cat(private$color.text(colnames(self$data$train$Y)[i], response.color), "\tSubjects: ", sum(!is.na(self$data$train$Y[,i])), "\tType: ", self$params$family, "\n")
    }
  }
  # ---------------------------------
  if(type == "help"){
    cat("For assistance, browse the SPEAR vignettes\nType browseVignettes('SPEAR')\n\n")
  }
  
}



#' Color text for output in the terminal
#'@param text A string to be colored
#'@param fg Foreground color
#'@param bg Background color
#'@return Text that has been formatted.
#' @examples
#' color.text("Error", fg = "red", bg = NULL)
#'@export
color.text <- function(text, fg = "black", bg = NULL) {
  if(self$options$remove.formatting){
    return(text)
  }
  term <- Sys.getenv()["TERM"]
  colour_terms <- c("xterm-color","xterm-256color", "screen", "screen-256color")
  .fg_colours <- c(
    "black" = "0;30",
    "blue" = "0;34",
    "green" = "0;32",
    "cyan" = "0;36",
    "red" = "0;31",
    "purple" = "0;35",
    "brown" = "0;33",
    "light gray" = "0;37",
    "dark gray" = "1;30",
    "light blue" = "1;34",
    "light green" = "1;32",
    "light cyan" = "1;36",
    "light red" = "1;31",
    "light purple" = "1;35",
    "yellow" = "1;33",
    "white" = "1;37"
  )
  .bg_colours <- c(
    "black" = "40",
    "red" = "41",
    "green" = "42",
    "brown" = "43",
    "blue" = "44",
    "purple" = "45",
    "cyan" = "46",
    "light gray" = "47"
  )
  if(nchar(Sys.getenv('R_TESTS')) != 0 || !any(term %in% colour_terms, na.rm = TRUE)) {
    return(text)
  }
  col_escape <- function(col) {
    paste0("\033[", col, "m")
  }
  col <- .fg_colours[tolower(fg)]
  if (!is.null(bg)) {
    col <- paste0(col, .bg_colours[tolower(bg)], sep = ";")
  }
  init <- col_escape(col)
  reset <- col_escape("0")
  return(paste0(init, text, reset))
}