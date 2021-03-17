### SPEARkit Functions ###

#' Get best SPEAR weights per response Y
#'@param SPEARobj SPEAR object (returned from \code{run_cv_spear})
#'@param w.method How to choose best weight? Options include \code{"min"} (lowest mean CV error, default) and \code{"sd"} (choose a higher weight within 1 standard deviation of the CV errors).
#'@param return.overall Should the best overall weight be returned? Defaults to \code{FALSE} (where each response can have a different best weight).
#'@param return.as.index Return the index (from SPEARobj$params$weights) rather than the weight itself? Defaults to \code{FALSE}.
#'@return A named vector of weights (or weight indices from \code{SPEARobj$params$weights} if \code{return.as.index} == \code{TRUE}).
#' @examples
#' SPEAR.get_best_weights(SPEARobj, w.method = "min", return.overall = TRUE)
#' SPEAR.get_best_weights(SPEARobj, w.method = "sd")
#'@export
SPEAR.get_best_weights <- function(SPEARobj, w.method = "sd", return.overall = FALSE, return.as.index = FALSE){
  # If only one column, "return.overall" means "best"
  if(ncol(SPEARobj$data$Y) == 1){
    return.overall = FALSE
  }
  if(w.method == "min"){
    if(!return.overall){
      w.idxs = apply(SPEARobj$cv.eval$cvm, 2, which.min)
      ws <- SPEARobj$params$weights[w.idxs]
    }
    else if(return.overall){
      w.idx = which.min(apply(SPEARobj$cv.eval$cvm,1,sum))
      w.idxs <- rep(w.idx, ncol(SPEARobj$data$Y))
      ws <- SPEARobj$params$weights[w.idxs]
    }
  } else if(w.method == "sd"){
    if(!return.overall){
      min.ws.idxs = apply(SPEARobj$cv.eval$cvm, 2, which.min)
      min.ws <- sapply(1:length(min.ws.idxs), function(idx){return(SPEARobj$cv.eval$cvm[min.ws.idxs[idx],idx])})
      min.sd <- sapply(1:length(min.ws.idxs), function(idx){return(SPEARobj$cv.eval$cvsd[min.ws.idxs[idx],idx])})
      min.plus.sd.ws <- min.ws + min.sd
      w.idxs <- sapply(1:length(min.ws.idxs), function(idx){
        w.idx <- min.ws.idxs[idx]
        best.idx <- w.idx
        best.w <- SPEARobj$params$weights[w.idx]
        for(t in 1:nrow(SPEARobj$cv.eval$cvm)){
          if(SPEARobj$cv.eval$cvm[t,idx] < min.plus.sd.ws[idx] & SPEARobj$params$weights[t] > best.w){
            best.idx <- t
            best.w <- SPEARobj$params$weights[t]
          }
        }
        return(best.idx)
      })
      ws <- SPEARobj$params$weights[w.idxs]
    }
    else if(return.overall){
      min.ws.idx <- which.min(rowSums(SPEARobj$cv.eval))
      min.ws.sum <- rowSums(SPEARobj$cv.eval$cvm[min.ws.idx,])
      min.ws.sd <- rowSums(SPEARobj$cv.eval$cvsd[min.ws.idx,])
      min.plus.sd.ws <- min.ws.sum + min.ws.sd
      best.idx <- w.idx
      best.w <- SPEARobj$params$weights[w.idx]
      for(t in 1:nrow(SPEARobj$cv.eval$cvm)){
        if(rowSums(SPEARobj$cv.eval$cvm[t,]) < min.plus.sd.ws & SPEARobj$params$weights[t] > best.w){
          best.idx <- t
          best.w <- SPEARobj$params$weights[t]
        }
      }
      w.idxs <- rep(best.idx, ncol(SPEARobj$data$Y))
      ws <- SPEARobj$params$weights[w.idxs]
    }
  } else {
    stop("ERROR: w.method defined incorrectly. Acceptable values are 'min' and 'sd'.")
  }
  # Return as index?
  if(return.as.index){
    output <- w.idxs
  } else {
    output <- ws
  }
  names(output) <- colnames(SPEARobj$data$Y)
  return(output)
}

#' Predict sample responses using a trained SPEAR model.
#'@param SPEARobj SPEAR object (returned from \code{run_cv_spear})
#'@param X List of omics matrices to be used for prediction. Defaults to \code{NULL} (use training data stored within SPEAR)
#'@param w Weight for SPEAR model. Defaults to \code{"best"}, choosing the best weight per response. Can also be \code{"overall"} (choosing the weight with the best overall mean cross-validated error), or one of the weights used to train SPEAR (\code{SPEARobj$params$weights})
#'@param w.method How to choose best weight? Options include \code{"min"} (lowest mean CV error, default) and \code{"sd"} (choose a higher weight within 1 standard deviation of the CV errors).
#'@param scale.x Should \code{X} be scaled (only used if new \code{X} is supplied). Defaults to \code{FALSE}
#'@return A list of predictions from the SPEAR model. One list for each response predicted, as well as \code{predictions}, \code{probabilites}, and \code{w} used.
#' @examples
#' SPEAR.get_predictions(SPEARobj, X = list.of.test.matrices, scale.x = TRUE)
#' SPEAR.get_predictions(SPEARobj, w = 1)
#'@export
SPEAR.get_predictions <- function(SPEARobj, X = NULL, w = "best", w.method = "sd", scale.x = FALSE){
  # Weight ----------------
  if(w == "best"){
    w.idxs <- SPEAR.get_best_weights(SPEARobj = SPEARobj, w.method = w.method, return.overall = FALSE, return.as.index = TRUE)
  } else if(w == "overall"){
    w.idxs <- SPEAR.get_best_weights(SPEARobj = SPEARobj, w.method = w.method, return.overall = TRUE, return.as.index = TRUE)
  } else {
    if(!any(SPEARobj$params$weights == w)){
      cat("*** Warning: w = ", w, " not found among possible weights (", paste(SPEARobj$params$weights, collapse = ", "), "). Will use closest SPEAR weight...\n")
      w = SPEARobj$params$weights[which.min(abs(SPEARobj$params$weights-w))]
      w.idxs <- which(SPEARobj$params$weights == w)
      cat(paste0("*** Using SPEAR with w = ", round(SPEARobj$params$weights[w.idxs], 3), "\n"))
      w.idxs <- rep(w.idxs, ncol(SPEARobj$data$Y))
    } else {
      w.idxs <- which(SPEARobj$params$weights == w)
      cat(paste0("*** Using SPEAR with w = ", round(SPEARobj$params$weights[w.idxs], 3), "\n"))
      w.idxs <- rep(w.idxs, ncol(SPEARobj$data$Y))
    }
  } #----------------------
  
  # Check if X is NULL (use training data):
  if(is.null(X)){
    X <- SPEARobj$data$X
  } else {
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
    X <- do.call("cbind", X)
    # Scale X by feature:
    if(scale.x){
      X <- scale(X)
    }
  }
  
  # Get Predictions:
  if(SPEARobj$params$family == "gaussian"){
    output <- list()
    for(i in 1:length(w.idxs)){
      response <- colnames(SPEARobj$data$Y)[i]
      w.idx <- w.idxs[i]
      ### Group loop? [,g,,]
      preds <- X %*% SPEARobj$cv.eval$reg_coefs[,,i,w.idx] + SPEARobj$cv.eval$intercepts[[1]][w.idx,]
      names(preds) <- rownames(X)
      ### Add to group loop?
      output[[response]] <- list(predictions = preds, w = SPEARobj$params$weights[w.idx])
    }
    
  } else {
    output <- list()
    for(i in 1:length(w.idxs)){
      response <- colnames(SPEARobj$data$Y)[i]
      w.idx <- w.idxs[i]
      ### Group loop? [,g,,]
      yhat = X %*% SPEARobj$cv.eval$reg_coefs[,,i,w.idx]
      intercept = SPEARobj$cv.eval$intercepts[[1]][w.idx,]
      Pmat0 = matrix(0, ncol = length(intercept), nrow = length(yhat))
      Pmat = matrix(0, ncol = length(intercept) + 1, nrow = length(yhat))
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
      # Get predictions (highest probability)
      predictions <- apply(Pmat, 1, which.max) - 1
      rownames(Pmat) <- rownames(SPEARobj$data$Y)
      colnames(Pmat) <- paste0("Class", 0:(ncol(Pmat)-1))
      names(predictions) <- rownames(X)
      ### Add to group loop?
      output[[response]] <- list(probabilities = Pmat, predictions = predictions, w = SPEARobj$params$weights[w.idxs])
    }
  }
  
  return(output)
}

#' Get loss for each value of \code{w} from training a SPEAR model.
#'@param SPEARobj SPEAR object (returned from \code{run_cv_spear})
#'@param include.sd Return standard deviation of loss as well? Defaults to \code{FALSE}.
#'@return A matrix with the first column being each \code{w} trained on and a column for every response. If \code{include.sd == TRUE}, a list of matrices, named \code{mean} and \code{sd}.
#'@examples
#' SPEAR.get_cv_loss(SPEARobj)
#' SPEAR.get_cv_loss(SPEARobj, include.sd = TRUE)
#'@export
SPEAR.get_cv_loss <- function(SPEARobj, include.sd = FALSE){
  cv.errors <- SPEARobj$cv.eval$cvm
  cv.errors <- cbind(SPEARobj$params$weights, cv.errors)
  colnames(cv.errors) <- c("w", colnames(SPEARobj$data$Y))
  if(include.sd){
    cv.sd <- SPEARobj$cv.eval$cvsd
    cv.sd <- cbind(SPEARobj$params$weights, cv.sd)
    colnames(cv.sd) <- c("w", colnames(SPEARobj$data$Y))
    return(list(mean = cv.errors, sd = cv.sd))
  } else {
    return(cv.errors)
  }
}

#' Get factor scores from a trained SPEAR object on training/test data
#'@param SPEARobj SPEAR object (returned from \code{run_cv_spear})
#'@param X List of omics matrices to be used for prediction. Defaults to \code{NULL} (use training data stored within SPEAR)
#'@param w Weight for SPEAR model. Defaults to \code{"best"}, choosing the best weight per response. Can also be \code{"overall"} (choosing the weight with the best overall mean cross-validated error), or one of the weights used to train SPEAR (\code{SPEARobj$params$weights})
#'@param w.method How to choose best weight? Options include \code{"min"} (lowest mean CV error, default) and \code{"sd"} (choose a higher weight within 1 standard deviation of the CV errors).
#'@param scale.x Should \code{X} be scaled (only used if new \code{X} is supplied). Defaults to \code{FALSE}
#'@return A list of predictions from the SPEAR model. One list for each response predicted, as well as \code{predictions}, \code{probabilites}, and \code{w} used.
#' @examples
#' SPEAR.get_factor_scores(SPEARobj)
#' SPEAR.get_factor_scores(SPEARobj, w = 1)
#'@export
SPEAR.get_factor_scores <- function(SPEARobj, w = "best", X = NULL, w.method = "min", scale.x = FALSE){
  # Weight ----------------
  if(w == "best"){
    w.idxs <- SPEAR.get_best_weights(SPEARobj = SPEARobj, w.method = w.method, return.overall = FALSE, return.as.index = TRUE)
  } else if(w == "overall"){
    w.idxs <- SPEAR.get_best_weights(SPEARobj = SPEARobj, w.method = w.method, return.overall = TRUE, return.as.index = TRUE)
  } else {
    if(!any(SPEARobj$params$weights == w)){
      cat("*** Warning: w = ", w, " not found among possible weights (", paste(SPEARobj$params$weights, collapse = ", "), "). Will use closest SPEAR weight...\n")
      w = SPEARobj$params$weights[which.min(abs(SPEARobj$params$weights-w))]
      w.idxs <- which(SPEARobj$params$weights == w)
      cat(paste0("*** Using SPEAR with w = ", round(SPEARobj$params$weights[w.idxs], 3), "\n"))
      w.idxs <- rep(w.idxs, ncol(SPEARobj$data$Y))
    } else {
      w.idxs <- which(SPEARobj$params$weights == w)
      cat(paste0("*** Using SPEAR with w = ", round(SPEARobj$params$weights[w.idxs], 3), "\n"))
      w.idxs <- rep(w.idxs, ncol(SPEARobj$data$Y))
    }
  } #----------------------
  if(is.null(X)){
    X <- SPEARobj$data$X
  } else {
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
    X <- do.call("cbind", X)
    # Scale X by feature:
    if(scale.x){
      X <- scale(X)
    }
  }
  results <- list()
  for(i in 1:length(w.idxs)){
    response <- colnames(SPEARobj$data$Y)[i]
    w.idx <- w.idxs[i]
    ### Group loop? [,,g,]
    Uhat <- X %*% SPEARobj$fit$results$post_betas[,,,w.idx]
    colnames(Uhat) <- paste0("Factor", 1:SPEARobj$params$num_factors)
    rownames(Uhat) <- rownames(X)
    temp <- list()
    temp$factor.scores <- Uhat
    temp$weight <- SPEARobj$params$weights[w.idx]
    ### Add to group loop?
    results[[response]] <- temp
  }
  return(results)
}


#' Plot factor loadings
#'@param SPEARobj SPEAR object (returned from run_cv_spear)
#'@param w Weight for SPEAR. Defaults to "best", choosing the best weight per response. Can also be "overall" (choosing the weight with the best overall mean cross-validated error), or one of the weights used to train SPEAR (SPEARobj$params$weights)
#'@param threshold Threshold value of contribution for a factor to be relevant. Defaults to .01.
#'@param response.name Which response variable? Check colnames(SPEARobj$data$Y) for options. Defaults to 1st column (NULL)
#'@param plot.per.omic Return an individual plot per omic (each dataset in SPEARobj$data$xlist)? Defaults to TRUE
#'@param scale.per.omic Return an individual plot per omic (each dataset in SPEARobj$data$xlist)? Defaults to TRUE
#'@param show.feature.names Show row names (features) for factor loadings? Defaults to FALSE
#'@param plot.irrelevant.factors Plot loadings for factors with a contribution lower than the threshold? Defaults to FALSE
#'@export
SPEAR.plot_factor_loadings <- function(SPEARobj, w = "best", w.method = "sd", threshold = .01, response.name = NULL, plot.per.omic = TRUE, scale.per.omic = TRUE, show.feature.names = FALSE, plot.irrelevant.factors = FALSE){
  # Weight ----------------
  if(w == "best"){
    w.idxs <- SPEAR.get_best_weights(SPEARobj = SPEARobj, w.method = w.method, return.overall = FALSE, return.as.index = TRUE)
  } else if(w == "overall"){
    w.idxs <- SPEAR.get_best_weights(SPEARobj = SPEARobj, w.method = w.method, return.overall = TRUE, return.as.index = TRUE)
  } else {
    if(!any(SPEARobj$params$weights == w)){
      cat("*** Warning: w = ", w, " not found among possible weights (", paste(SPEARobj$params$weights, collapse = ", "), "). Will use closest SPEAR weight...\n")
      w = SPEARobj$params$weights[which.min(abs(SPEARobj$params$weights-w))]
      w.idxs <- which(SPEARobj$params$weights == w)
      cat(paste0("*** Using SPEAR with w = ", round(SPEARobj$params$weights[w.idxs], 3), "\n"))
      w.idxs <- rep(w.idxs, ncol(SPEARobj$data$Y))
    } else {
      w.idxs <- which(SPEARobj$params$weights == w)
      cat(paste0("*** Using SPEAR with w = ", round(SPEARobj$params$weights[w.idxs], 3), "\n"))
      w.idxs <- rep(w.idxs, ncol(SPEARobj$data$Y))
    }
  } #----------------------
  
  # response.name: ------------
  if(is.null(response.name)){
    if(ncol(SPEARobj$data$Y) == 1){
      response.name <- colnames(SPEARobj$data$Y)
    } else {
      cat("\nParameter 'response.name' not provided. Using first response variable ('", colnames(SPEARobj$data$Y)[1], "')\n")
      response.name <- colnames(SPEARobj$data$Y)[1]
    }
  } else if(!any(response.name %in% colnames(SPEARobj$data$Y))){
    stop(paste0("ERROR: 'response.name' provided ('", response.name, "') is not found among the available response variables.\nRequested:\t'", response.name, "'\nAvailable:\t'", paste(colnames(SPEARobj$data$Y), collapse = "', '"), "'"))
  }
  # --------------------------
  
  plot.list.total <- list()
  k <- which(colnames(SPEARobj$data$Y) == response.name)
  w <- SPEARobj$params$weights[w.idxs[k]]
  relevant.factors <- SPEARobj$cv.eval$factor_contributions[, k, w.idxs[k]] >= threshold
  names(relevant.factors) <- paste0("Factor", 1:length(relevant.factors))
  w.loadings = SPEARobj$fit$results$post_selections[,,w.idxs[k]]
  if(!plot.irrelevant.factors){
    w.loadings[,!relevant.factors] <- NA
  }
  colnames(w.loadings) <- paste0("Factor", 1:ncol(w.loadings))
  rownames(w.loadings) <- colnames(SPEARobj$data$X)
  # Generate plots
  if(!plot.per.omic){
    # Plot all together
    df <- reshape2::melt(w.loadings)
    colnames(df)[3] <- "Probability"
    g.temp <- ggplot(df) +
      geom_tile(aes(x = Var2, y = Var1, fill = Probability)) +
      scale_fill_distiller(palette = "Reds", direction = 1, na.value="#F8F8F8") +
      ylab(paste0("Features\n(n = ", ncol(SPEARobj$data$X), ")")) +
      xlab(NULL) +
      cowplot::theme_map() +
      theme(axis.title.y = element_text(size = 10, angle = 90),
            legend.text = element_text(size = 8),
            legend.title = element_blank(),
            plot.margin = unit(c(0, 0, 0, 0), "cm"),
            axis.text.x = element_text(size = 10))
    
    if(!show.feature.names){
      g.temp <- g.temp + theme(axis.text.y=element_blank(),
                               axis.ticks.y=element_blank())
    }
    
    # Return g.temp
    p <- g.temp
  }
  else{
    # Split w.loadings into separate omics
    num.omics <- length(SPEARobj$data$xlist)
    colors <- rep(c("Reds", "Blues", "Greens", "Greys"), length.out = num.omics)
    s <- 1
    plot.list <- list()
    for(o in 1:num.omics){
      w.loadings.current <- w.loadings[s:(s + ncol(SPEARobj$data$xlist[[o]]) - 1),]
      s <- s + ncol(SPEARobj$data$xlist[[o]])
      # Make a separate plot for each:
      df <- reshape2::melt(w.loadings.current)
      colnames(df)[3] <- "Probability"
      g.temp <- ggplot(df) +
        geom_tile(aes(x = Var2, y = Var1, fill = Probability)) +
        scale_fill_distiller(palette = colors[o], direction = 1, na.value="#F8F8F8") +
        ylab(paste0(names(SPEARobj$data$xlist)[o], "\n(n = ", ncol(SPEARobj$data$xlist[[o]]), ")")) +
        xlab(NULL) +
        cowplot::theme_map() +
        theme(axis.title.y = element_text(size = 10, angle = 90),
              legend.text = element_text(size = 8),
              legend.title = element_blank(),
              plot.title = element_text(size = 10, hjust = .5, face = "plain"),
              plot.margin = unit(c(0, 0, 0, 0), "cm"))
      
      if(o == num.omics){
        g.temp <- g.temp + theme(axis.text.x = element_text(size = 10))
      }
      
      if(!show.feature.names){
        g.temp <- g.temp + theme(axis.text.y=element_blank(),
                                 axis.ticks.y=element_blank())
      }
      
      # Save
      plot.list[[o]] <- g.temp
    }
    # return list:
    if(scale.per.omic){
      rel.heights <- sapply(SPEARobj$data$xlist, ncol)
    } else {
      rel.heights <- rep(1, num.omics)
    }
    p <- cowplot::plot_grid(plotlist = plot.list, ncol = 1, rel_heights = rel.heights)
    title_gg <- ggplot() + 
      labs(title = paste0(response.name, " | w = ", round(w, 3))) +
      cowplot::theme_map() +
      theme(plot.title = element_text(hjust = .5, vjust = .5, face = "plain", size = 14))
    p <- cowplot::plot_grid(title_gg, p, ncol = 1, rel_heights = c(0.05, 1))
  }
  return(p)
}


#' Plot ordinal class probabilities
#'@param SPEARobj SPEAR object (returned from run_cv_spear)
#'@param X List of omics matrices to be used for prediction. Defaults to NULL (use training data stored within SPEAR)
#'@param Y Responses for subjects (rows) in X. Defaults to NULL (use training data stored within SPEAR)
#'@param w Weight for SPEAR. Defaults to "best", choosing the best weight per response. Can also be "overall" (choosing the weight with the best overall mean cross-validated error), or one of the weights used to train SPEAR (SPEARobj$params$weights)
#'@param response.name Which response variable? Check colnames(SPEARobj$data$Y) for options. Defaults to 1st column (NULL)
#'@param scale.x Should X be scaled (only used if X is supplied). Defaults to FALSE
#'@param plot.type Plot type - "boxplot" or "line". Defaults to "boxplot"
#'@export
SPEAR.plot_class_probabilities <- function(SPEARobj, X = NULL, Y = NULL, w = "best", response.name = NULL, scale.x = FALSE, plot.type = "boxplot"){
  ## Gaussian check: ----------
  if(SPEARobj$params$family == "gaussian"){
    stop("ERROR: SPEARobj must be of 'family' type 'binomial', 'ordinal', or 'categorical'. Current SPEARobj is of type 'gaussian'.")
  }
  
  # response.name: ------------
  if(is.null(response.name) & is.null(Y)){
    if(ncol(SPEARobj$data$Y) == 1){
      response.name <- colnames(SPEARobj$data$Y)
    } else {
      cat("\nParameter 'response.name' not provided. Using first response variable ('", colnames(SPEARobj$data$Y)[1], "')\n")
      response.name <- colnames(SPEARobj$data$Y)[1]
    }
  } else if(!any(response.name %in% colnames(SPEARobj$data$Y))){
    stop(paste0("ERROR: 'response.name' provided ('", response.name, "') is not found among the available response variables.\nRequested:\t'", response.name, "'\nAvailable:\t'", paste(colnames(SPEARobj$data$Y), collapse = "', '"), "'"))
  }
  # --------------------------
  
  if(!is.null(Y)){
    if(nrow(Y) != nrow(X[[1]])){
      stop(paste0("ERROR: Y provided has ", nrow(Y), " rows, and X has ", nrow(X[[1]]), " rows (need to match)."))
    }
    if(ncol(Y) > 1){
      if(is.null(colnames(Y))){
        stop(paste0("ERROR: Y provided has no colnames. Columns need to have matching colnames from any of colnames(SPEARobj$data$Y)."))
      } else if(!any(colnames(Y) %in% colnames(SPEARobj$data$Y))){
        stop(paste0("ERROR: One or more column names provided in Y are not found among the trained response.names.\ncolnames(Y) provided:\t'", paste(colnames(Y), collapse = "', '"), "'\ncolnames(SPEARobj$data$Y):\t'", paste(colnames(SPEARobj$data$Y), collapse = "', '"), "'"))
      }
      cat("*** Warning: Y provided has more than one response variable... Using first column only (", colnames(Y)[1], ")")
      response.name <- colnames(Y)[1]
      k <- which(response.name %in% colnames(SPEARobj$data$Y))
      Y <- Y[,k]
    }
  } else {
    k <- which(response.name %in% colnames(SPEARobj$data$Y))
    Y <- SPEARobj$data$Y[,k]
  }
  
  
  
  preds <- SPEAR.get_predictions(SPEARobj = SPEARobj, X = X, w = w, scale.x = scale.x)
  levels <- ncol(preds[[response.name]]$probabilities)
  df <- cbind(as.data.frame(preds[[response.name]]$probabilities), Y)
  colnames(df) <- c(paste0("ProbClass", 1:levels), "Y")
  df$Ychar <- as.character(df$Y)
  
  if(plot.type == "boxplot")
  {
    plotlist <- list()
    for(i in 1:levels){
      g <- ggplot(df) +
        geom_boxplot(aes_string(x = "Y", group = "Y", y = paste0("ProbClass", i), fill = "Ychar"), lwd = .25, outlier.size = 1.5, outlier.shape = 21, alpha = .6) +
        xlab("True Label") +
        ylab("Probability") +
        geom_segment(aes(x = -0.5, y = 0, xend = (levels-.5), yend = 0), lwd = 0) +
        scale_fill_brewer(palette = "Spectral", guide = FALSE) +
        scale_x_continuous(labels = c(0:(levels-1)), breaks = c(0:(levels-1))) +
        ggtitle(paste0("Probability Class ", i-1)) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
      
      plotlist[[length(plotlist) + 1]] <- g
    }
    
    
    
    p <- plot_grid(plotlist = plotlist, nrow = 1)
    return(p)
  } else if(plot.type == "line"){
    probabilities <- preds[[response.name]]$probabilities
    predictions <- preds[[response.name]]$predictions
    subject.probabilities <- as.data.frame(probabilities)
    subject.probabilities <- dplyr::mutate(subject.probabilities, 
                                           PredictedClass = as.factor(predictions),
                                           Actual = paste0("Class", SPEARobj$data$Y[,k]),
                                           Subject = rownames(subject.probabilities))
    df.melt <- reshape2::melt(subject.probabilities, id.vars = c("PredictedClass", "Actual", "Subject"))
    g <- ggplot(df.melt) +
      geom_line(aes(x = variable, y = value, group = Subject, color = PredictedClass)) +
      geom_point(aes(x = variable, y = value, group = Subject, color = PredictedClass)) +
      scale_x_discrete(labels = (1:ncol(subject.probabilities)-1)) +
      xlab(NULL) +
      ylab("Probabillity of Class Assignment") +
      scale_color_brewer(palette = "Spectral", direction = 1) +
      scale_fill_brewer(palette = "Spectral", direction = 1) +
      theme_classic() +
      facet_wrap(vars(Actual), nrow = 1) +
      theme(panel.background = element_rect(fill = "#181818"),
            panel.grid.major = element_line(size = 0.2, linetype = 'solid', colour = "#656464"))
    # Return plot:
    return(g)
    
  } else {
    stop("ERROR: plot.type '", plot.type, "' not recognized. Acceptable choices are 'boxplot' and 'line'.")
  }
}


#' Plot ordinal class predictions
#'@param SPEARobj SPEAR object (returned from run_cv_spear)
#'@param X List of omics matrices to be used for prediction. Defaults to NULL (use training data stored within SPEAR)
#'@param Y Responses for subjects (rows) in X. Defaults to NULL (use training data stored within SPEAR)
#'@param w Weight for SPEAR. Defaults to "best", choosing the best weight per response. Can also be "overall" (choosing the weight with the best overall mean cross-validated error), or one of the weights used to train SPEAR (SPEARobj$params$weights)
#'@param response.name Which response variable? Check colnames(SPEARobj$data$Y) for options. Defaults to 1st column (NULL)
#'@param scale.x Should X be scaled (only used if X is supplied). Defaults to FALSE
#'@export
SPEAR.plot_class_predictions <- function(SPEARobj, X = NULL, Y = NULL, w = "best", response.name = NULL, scale.x = FALSE, show.true.distribution = TRUE){
  ## Gaussian check: ----------
  if(SPEARobj$params$family == "gaussian"){
    stop("ERROR: SPEARobj must be of 'family' type 'binomial', 'ordinal', or 'categorical'. Current SPEARobj is of type 'gaussian'.")
  }
  
  # response.name: ------------
  if(is.null(response.name) & is.null(Y)){
    if(ncol(SPEARobj$data$Y) == 1){
      response.name <- colnames(SPEARobj$data$Y)
    } else {
      cat("\nParameter 'response.name' not provided. Using first response variable ('", colnames(SPEARobj$data$Y)[1], "')\n")
      response.name <- colnames(SPEARobj$data$Y)[1]
    }
  } else if(!any(response.name %in% colnames(SPEARobj$data$Y))){
    stop(paste0("ERROR: 'response.name' provided ('", response.name, "') is not found among the available response variables.\nRequested:\t'", response.name, "'\nAvailable:\t'", paste(colnames(SPEARobj$data$Y), collapse = "', '"), "'"))
  }
  # --------------------------
  
  if(!is.null(Y)){
    if(nrow(Y) != nrow(X[[1]])){
      stop(paste0("ERROR: Y provided has ", nrow(Y), " rows, and X has ", nrow(X[[1]]), " rows (need to match)."))
    }
    if(ncol(Y) > 1){
      if(is.null(colnames(Y))){
        stop(paste0("ERROR: Y provided has no colnames. Columns need to have matching colnames from any of colnames(SPEARobj$data$Y)."))
      } else if(!any(colnames(Y) %in% colnames(SPEARobj$data$Y))){
        stop(paste0("ERROR: One or more column names provided in Y are not found among the trained response.names.\ncolnames(Y) provided:\t'", paste(colnames(Y), collapse = "', '"), "'\ncolnames(SPEARobj$data$Y):\t'", paste(colnames(SPEARobj$data$Y), collapse = "', '"), "'"))
      }
      cat("*** Warning: Y provided has more than one response variable... Using first column only (", colnames(Y)[1], ")")
      response.name <- colnames(Y)[1]
      k <- which(response.name %in% colnames(SPEARobj$data$Y))
      Y <- Y[,k]
    }
  } else {
    k <- which(response.name %in% colnames(SPEARobj$data$Y))
    Y <- SPEARobj$data$Y[,k]
  }
  
  
  
  preds <- SPEAR.get_predictions(SPEARobj = SPEARobj, X = X, w = w, scale.x = scale.x)
  levels <- ncol(preds[[response.name]]$probabilities)
  temp <- tibble(pred = preds[[response.name]]$predictions, actual = unlist(Y))
  temp$correct <- temp$pred == temp$actual
  
  ymax <- max(table(temp$pred), table(temp$actual))
  
  p.true <- ggplot(temp) +
    geom_histogram(aes(x = actual, fill = as.character(actual)), stat = "count", lwd = .25, color = "black", alpha = .6) +
    xlab("True Label") +
    ylab("Count") +
    geom_segment(aes(x = -0.5, y = 0, xend = (levels-.5), yend = 0), lwd = 0) +
    scale_x_continuous(labels = c(0:(levels-1)), breaks = c(0:(levels-1))) +
    scale_fill_brewer(palette = "Spectral", guide = FALSE) +
    ggtitle("True Classes") +
    ylim(c(NA, ymax)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  p.pred <- ggplot(temp) +
    geom_histogram(aes(x = pred, fill = as.character(actual)), stat = "count", lwd = .25, color = "black", alpha = .6) +
    xlab("True Label") +
    ylab("Count") +
    geom_segment(aes(x = -0.5, y = 0, xend = (levels-.5), yend = 0), lwd = 0) +
    scale_x_continuous(labels = c(0:(levels-1)), breaks = c(0:(levels-1))) +
    scale_fill_brewer(palette = "Spectral", guide = FALSE) +
    ggtitle("SPEARordinal Predictions") +
    ylim(c(NA, ymax)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  if(show.true.distribution){
    p <- cowplot::plot_grid(p.true, p.pred, nrow = 1)
  } else {
    p <- p.pred
  }
  
  return(p)
}



#' Get factor contributions to both X and Y
#'@param SPEARobj SPEAR object (returned from run_cv_spear)
#'@param w Weight for SPEAR. Defaults to "best", choosing the best weight per response. Can also be "overall" (choosing the weight with the best overall mean cross-validated error), or one of the weights used to train SPEAR (SPEARobj$params$weights)
#'@param w.method How to choose best weight? Options include \code{"min"} (lowest mean CV error, default) and \code{"sd"} (choose a higher weight within 1 standard deviation of the CV errors).
#'@param threshold Threshold value of contribution for a factor to be relevant. Defaults to .01.
#'@export
SPEAR.get_factor_contributions <- function(SPEARobj, w = "best", w.method = "sd", threshold = .01){
  # Weight ----------------
  if(w == "best"){
    w.idxs <- SPEAR.get_best_weights(SPEARobj = SPEARobj, w.method = w.method, return.overall = FALSE, return.as.index = TRUE)
  } else if(w == "overall"){
    w.idxs <- SPEAR.get_best_weights(SPEARobj = SPEARobj, w.method = w.method, return.overall = TRUE, return.as.index = TRUE)
  } else {
    if(!any(SPEARobj$params$weights == w)){
      cat("*** Warning: w = ", w, " not found among possible weights (", paste(SPEARobj$params$weights, collapse = ", "), "). Will use closest SPEAR weight...\n")
      w = SPEARobj$params$weights[which.min(abs(SPEARobj$params$weights-w))]
      w.idxs <- which(SPEARobj$params$weights == w)
      cat(paste0("*** Using SPEAR with w = ", round(SPEARobj$params$weights[w.idxs], 3), "\n"))
      w.idxs <- rep(w.idxs, ncol(SPEARobj$data$Y))
    } else {
      w.idxs <- which(SPEARobj$params$weights == w)
      cat(paste0("*** Using SPEAR with w = ", round(SPEARobj$params$weights[w.idxs], 3), "\n"))
      w.idxs <- rep(w.idxs, ncol(SPEARobj$data$Y))
    }
  } #----------------------
  res <- list()
  for(w.idx in unique(w.idxs)){
    # Explanatory Vars:
    Xlist <- SPEARobj$data$xlist
    W <- SPEARobj$fit$results$post_betas[,,,w.idx]
    U <- SPEARobj$data$X %*% W
    Wlist <- list()
    ind <- 1
    for(d in 1:length(Xlist)){
      Wlist[[d]] <- W[ind:(ind - 1 + ncol(Xlist[[d]])),]
      ind <- ind + ncol(Xlist[[d]])
    }
    var_explained_X <- matrix(0, nrow = ncol(U), ncol = length(Xlist))
    for(i in 1:nrow(var_explained_X)){
      for(j in 1:ncol(var_explained_X)){
        a <- sum((Xlist[[j]] - U[,i] %*% t(Wlist[[j]][,i]))**2, na.rm = TRUE)
        b <- sum((Xlist[[j]])**2, na.rm = TRUE)
        var_explained_X[i,j] <- (1 - a/b)
      }
    }
    colnames(var_explained_X) <- names(Xlist)
    rownames(var_explained_X) <- paste0("Factor", 1:nrow(var_explained_X))
    relevant.factors.X = var_explained_X
    relevant.factors.X[var_explained_X >= threshold] <- 1
    relevant.factors.X[var_explained_X < threshold] <- 0
    # Response:
    relevant.factors.Y <- matrix(0,ncol = ncol(SPEARobj$data$Y), nrow = dim(SPEARobj$cv.eval$factor_contributions)[1])
    var_explained_Y <- matrix(0,ncol = ncol(SPEARobj$data$Y), nrow = dim(SPEARobj$cv.eval$factor_contributions)[1])
    for(i in 1:ncol(var_explained_Y)){
      var_explained_Y[,i] <- SPEARobj$cv.eval$factor_contributions[,i,w.idx]
    }
    relevant.factors.Y[var_explained_Y >= threshold] <- 1
    colnames(var_explained_Y) <- colnames(SPEARobj$data$Y)
    rownames(var_explained_Y) <- paste0("Factor", 1:nrow(var_explained_Y))
    colnames(relevant.factors.Y) <- colnames(var_explained_Y)
    rownames(relevant.factors.Y) <- rownames(var_explained_Y)
    res[[paste0("w", round(SPEARobj$params$weights[w.idx], 3))]] <- list(X.explained = var_explained_X, X.relevant.factors = relevant.factors.X, Y.explained = var_explained_Y, Y.relevant.factors = relevant.factors.Y, threshold = threshold)
  }
  return(res)
}

#' Get a color scheme for a SPEAR plot
#'@param SPEARobj SPEAR object (returned from run_cv_spear)
#'@param scale.name A character of "omic", "response", or "combined". Defaults to "combined".
#'@export
SPEAR.get_color_scheme <- function(SPEARobj, scale.name = "combined"){
  if(scale.name == "omic"){
    colors <- c("#9E0142", "#D53E4F", "#F46D43", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")
    colorvec <- colors[seq(1, length(colors), length.out = length(SPEARobj$data$xlist))]
    names(colorvec) <- names(SPEARobj$data$xlist)
    return(colorvec)
  } else if(scale.name == "response"){
    colors <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
    colorvec <- colors[seq(1, length(colors), length.out = ncol(SPEARobj$data$Y))]
    names(colorvec) <- colnames(SPEARobj$data$Y)
    return(colorvec)
  } else if(scale.name == "combined"){
    colors <- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")
    colorvec.omic <- colors[seq(1, length(colors), length.out = length(SPEARobj$data$xlist))]
    names(colorvec.omic) <- names(SPEARobj$data$xlist)
    colors <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
    colorvec.resp <- colors[seq(1, length(colors), length.out = ncol(SPEARobj$data$Y))]
    names(colorvec.resp) <- colnames(SPEARobj$data$Y)
    return(c(colorvec.omic, colorvec.resp))
  }
}


#' Plot factor contributions for X and Y
#'@param SPEARobj SPEAR object (returned from run_cv_spear)
#'@param w Weight for SPEAR. Defaults to "best", choosing the best weight per response. Can also be "overall" (choosing the weight with the best overall mean cross-validated error), or one of the weights used to train SPEAR (SPEARobj$params$weights)
#'@param threshold Threshold value of contribution for a factor to be relevant. Defaults to .01.
#'@param show.labels Show the contributions as geom_text on the plot? Defaults to TRUE
#'@param show.irrelevant Show all contributions, even those below the threshold? Defaults to FALSE
#'@export
SPEAR.plot_factor_contributions <- function(SPEARobj, w = "best", w.method = "sd", threshold = .01, show.labels = TRUE, show.irrelevant = FALSE){
  
  factor.contributions.total <- SPEAR.get_factor_contributions(SPEARobj = SPEARobj, w = w, w.method = w.method, threshold = threshold)
  plotlist <- list()
  for(idx in 1:length(factor.contributions.total)){
    w <- names(factor.contributions.total)[idx]
    factor.contributions <- factor.contributions.total[[w]]
    
    if(show.irrelevant){
      title <- paste0("Factor Contributions | ", w)
      factor.contributions$X.relevant.factors <- matrix(1, ncol = ncol(factor.contributions$X.relevant.factors), nrow = nrow(factor.contributions$X.relevant.factors))
      factor.contributions$Y.relevant.factors <- matrix(1, ncol = ncol(factor.contributions$Y.relevant.factors), nrow = nrow(factor.contributions$Y.relevant.factors))
    } else {
      title <- paste0("Factor Contributions > ", threshold, " | ", w)
    }
    
    # Combine DFs:
    df.X <- factor.contributions$X.explained * factor.contributions$X.relevant.factors
    df.Y <- factor.contributions$Y.explained * factor.contributions$Y.relevant.factors
    df.comb <- cbind(df.X, df.Y)
    # Var explained:
    df.var <- data.frame(omic = colnames(df.comb), var = colSums(df.comb), omic.num = 1:ncol(df.comb), label = round(ifelse(colSums(df.comb) < threshold, NA, colSums(df.comb)), 2))
    # Plot
    df.melt <- reshape2::melt(df.comb)
    df.melt$value.adj <- ifelse(df.melt$value < 0, 0, df.melt$value)
    df.melt$label <- ifelse(df.melt$value==0, NA, round(df.melt$value, 2))
    
    g <- ggplot(data = df.melt) +
      geom_tile(aes(y = Var2, x = Var1, alpha = value.adj, fill = Var2), color = "black") +
      geom_text(aes(y = Var2, x = Var1, label = label), color = ifelse(show.labels, "black", NA)) +
      geom_rect(data = df.var, aes(xmin = nrow(df.comb) + .5, xmax = nrow(df.comb) + 1.5, ymin = .5, ymax = nrow(df.var) + .5), color = "white", fill = "white") +
      geom_rect(data = df.var, aes(xmin = nrow(df.comb) + .5, xmax = nrow(df.comb) + .5 + var, ymin = omic.num - .5, ymax = omic.num + .5, fill = omic), color = "black") +
      geom_text(data = df.var, aes(y = omic, x = nrow(df.comb) + 1, label = label), color = "black", angle = 270) +
      geom_hline(yintercept = length(SPEARobj$data$xlist) + .5, size = 1) +
      scale_alpha_continuous(range = c(0, 1), guide = F) +
      scale_fill_manual(values = SPEAR.get_color_scheme(SPEARobj, "combined"), guide = FALSE) +
      ggtitle(title) +
      ylab("") +
      xlab("") +
      coord_equal() +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, size = 12))
    
    plotlist[[idx]] <- g
  }
  print(length(plotlist))
  return(cowplot::plot_grid(plotlist = plotlist, nrow = 1))
}

#' Plot factor features
#'@param SPEARobj SPEAR object (returned from run_cv_spear)
#'@param w Weight for SPEAR. Defaults to "best", choosing the best weight per response. Can also be "overall" (choosing the weight with the best overall mean cross-validated error), or one of the weights used to train SPEAR (SPEARobj$params$weights)
#'@param threshold Threshold value of contribution for a factor to be relevant. Defaults to .01.
#'@param probability.cutoff Posterior selection probability cutoff (0 - 1) for a feature to be selected. Defaults to .5
#'@param coefficient.cutoff Coefficient cutoff (positive) for a feature to be selected. Will use +/-. Defaults to .01
#'@export
SPEAR.plot_feature_overlap <- function(SPEARobj, w = "best", factor = NULL, response.name = NULL, threshold = .01, probability.cutoff = .5, coefficient.cutoff = .01){
  
  # response.name:
  if(is.null(response.name)){
    if(ncol(SPEARobj$data$Y) == 1){
      response.name <- colnames(SPEARobj$data$Y)
    } else {
      cat("\nParameter 'response.name' not provided. Using first response variable ('", colnames(SPEARobj$data$Y)[1], "')\n")
      response.name <- colnames(SPEARobj$data$Y)[1]
    }
  } else if(!any(response.name %in% colnames(SPEARobj$data$Y))){
    stop(paste0("ERROR: 'response.name' provided ('", response.name, "') is not found among the available response variables.\nRequested:\t'", response.name, "'\nAvailable:\t'", paste(colnames(SPEARobj$data$Y), collapse = "', '"), "'"))
  }
  
  features <- SPEAR.get_factor_features(SPEARobj, w, threshold, probability.cutoff)
  
  if(is.null(factor)){
    factor <- names(features[[response.name]])
  } else {
    factor <- paste0("Factor", factor)
    if(any(!factor %in% names(features[[response.name]]))){
      stop(paste0("ERROR: one or more factors requested not relevant or are missing for Y = '", response.name, "'.\nRequested:\t'", paste(factor, collapse = "', '"), "'\nAvailable:\t'", paste(names(features[[response.name]]), collapse = "', '"), "'"))
    }
  }
  # Get all features per Omic:
  v.list <- list()
  for(o in 1:length(SPEARobj$data$xlist)){
    l <- lapply(factor, function(f){
      return(features[[response.name]][[f]][[names(SPEARobj$data$xlist)[o]]]$features)
    })
    groups <- factor
    names(l) <- factor
    color_vec <- c("red", "purple", "orange")
    v.table <- gplots::venn(l, universe = colnames(SPEARobj$data$xlist[[o]]), show.plot = FALSE, intersections = FALSE)
    # Set the attributes to NULL, make data.frame (otherwise you can't melt it to plot the points underneath)
    attr(v.table, "class") <- NULL
    v.table <- as.data.frame(v.table)
    # I add a "rowSums" column for sorting purposes (i.e. 111, 110, 101, 011, 100, 010, 001...)
    v.table$rowSums <- rowSums(dplyr::select(v.table, -num))
    # Filter out the case with 000 (we didn't supply a 'universe' parameter, so '000' can be removed)
    v.table <- dplyr::filter(v.table, rowSums > 0) %>% arrange(rowSums) %>% select(-rowSums)
    # Add 'comb' to v.table (rownames are "000"..."111")
    v.table <- dplyr::mutate(v.table, comb = rownames(v.table))
    v.table$omic <- names(SPEARobj$data$xlist)[o]
    v.table$omic_comb <- paste0(v.table$omic, "_", v.table$comb)
    # Melt so we can plot the individual points (can use pivot_longer if preferred)
    v.table.melt <- reshape2::melt(v.table, id.vars = c("comb", "omic", "omic_comb"))
    # Convert group names into indices, remove the 'num' column
    v.table.melt <- dplyr::filter(v.table.melt, variable != "num")
    v.table.melt$ypos <- sapply(v.table.melt$variable, function(var){
      return(-which(groups == var))
    })
    # Shift the ypos up .3 (looks better in my opinion)
    v.table.melt$ypos <- v.table.melt$ypos + .3
    
    v.list[[o]] <- list(v.table = v.table, v.table.melt = v.table.melt)
  }
  v.table.comb <- do.call(rbind, lapply(v.list, function(v){return(v$v.table)}))
  v.table.melt.comb <- do.call(rbind, lapply(v.list, function(v){return(v$v.table.melt)}))
  
  # Plots:
  
  # Generate plot:
  v.table <- v.table.comb
  v.table.melt <- v.table.melt.comb
  
  text.scale <- max(v.table$num)/20
  y.scale <- max(v.table$num)/10
  vlines <- 0:length(unique(v.table$comb))+.5
  xmin <- min(vlines)
  xmax <- max(vlines)
  ymin <- min(v.table.melt$ypos * y.scale) + max(v.table.melt$ypos * y.scale)
  ymax <- max(v.table$num + 2*text.scale)
  vlinedf <- data.frame(x = vlines, xend = vlines, y = ymin, yend = ymax)
  v.table.melt.total <- v.table.melt
  v.table.melt <- filter(v.table.melt, value > 0)
  
  g <- ggplot(v.table) +
    geom_bar(aes(x = factor(comb, levels = rev(unique(comb))), y = num, fill = omic), color = "black", position=position_dodge(width=.95), stat = "identity") +
    geom_text(aes(x = factor(comb, levels = rev(unique(comb))), y = num + text.scale, label = num, group = omic), position=position_dodge(width=1)) +
    annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = 0, color = "black", fill="#EDEDED") +
    geom_line(data = v.table.melt.total, aes(x = factor(comb, levels = rev(unique(comb))), y = ypos * y.scale), size = 2, color = "white") +
    geom_point(data = v.table.melt.total, aes(x = factor(comb, levels = rev(unique(comb))), y = ypos * y.scale), fill = "white", color = "white", shape = 21, size = 4, stroke = 1) +
    geom_line(data = v.table.melt, aes(x = factor(comb, levels = rev(unique(comb))), y = ypos * y.scale), size = 2) +
    geom_point(data = v.table.melt, aes(x = factor(comb, levels = rev(unique(comb))), y = ypos * y.scale), fill = "#58D68D", shape = 21, size = 4, stroke = 1) +
    geom_segment(data = vlinedf, aes(x = x, xend = xend, y = y, yend = yend), color = "black") +
    geom_segment(aes(x = xmin, xend = xmax, y = ymax, yend = ymax), color = "black") +
    scale_y_continuous(breaks = unique(v.table.melt$ypos * y.scale), labels = unique(v.table.melt$variable)) +
    xlab(NULL) +
    ylab(NULL) +
    ggtitle(paste0("Features for Factors Significant for ", response.name, " | w = ", w)) +
    cowplot::theme_map() +
    theme(axis.text.y = element_text(size = 10),
          title = element_text(size = 12, face = "plain")) +
    scale_fill_brewer(palette = "RdBu") +
    guides(size = FALSE, color = FALSE)
  
  return(g)
}


#' Plot factor features
#'@param SPEARobj SPEAR object (returned from run_cv_spear)
#'@param w Weight for SPEAR. Defaults to "best", choosing the best weight per response. Can also be "overall" (choosing the weight with the best overall mean cross-validated error), or one of the weights used to train SPEAR (SPEARobj$params$weights)
#'@param factor Which factors to plot? Defaults to NULL (all relevant factors).
#'@param response.name Which response variable? Check colnames(SPEARobj$data$Y) for options. Defaults to 1st column (NULL)
#'@param threshold Threshold value of contribution for a factor to be relevant. Defaults to .01.
#'@param probability.cutoff Posterior selection probability cutoff (0 - 1) for a feature to be selected. Defaults to .5
#'@param coefficient.cutoff Coefficient cutoff (positive) for a feature to be selected. Will use +/-. Defaults to .01
#'@param show.top.cutoff How many features to show? Set to NULL to show all. Defaults to 50.
#'@param group.by.omic Group features from the same omic together? Defaults to FALSE
#'@export
SPEAR.plot_feature_summary <- function(SPEARobj, w = "best", factor = NULL, response.name = NULL, threshold = .01, probability.cutoff = .5, coefficient.cutoff = 0, show.top.cutoff = 50, group.by.omic = FALSE){
  # response.name: ------------
  if(is.null(response.name)){
    if(ncol(SPEARobj$data$Y) == 1){
      response.name <- colnames(SPEARobj$data$Y)
    } else {
      cat("\nParameter 'response.name' not provided. Using first response variable ('", colnames(SPEARobj$data$Y)[1], "')\n")
      response.name <- colnames(SPEARobj$data$Y)[1]
    }
  } else if(!any(response.name %in% colnames(SPEARobj$data$Y))){
    stop(paste0("ERROR: 'response.name' provided ('", response.name, "') is not found among the available response variables.\nRequested:\t'", response.name, "'\nAvailable:\t'", paste(colnames(SPEARobj$data$Y), collapse = "', '"), "'"))
  }
  # --------------------------
  
  features <- SPEAR.get_factor_features(SPEARobj, w, threshold, probability.cutoff)
  # make sure there are factors that contribute:
  if(length(features[[response.name]]) == 1){
    stop("ERROR: No contributing factors found for 'response.name' ", response.name, ".")
  }
  
  # For each factor...
  if(is.null(factor)){
    factor <- names(features[[response.name]])
  } else {
    factor <- paste0("Factor", factor)
    if(any(!factor %in% names(features[[response.name]]))){
      stop(paste0("ERROR: one or more factors requested not relevant or are missing for Y = '", response.name, "'.\nRequested:\t'", paste(factor, collapse = "', '"), "'\nAvailable:\t'", paste(names(features[[response.name]]), collapse = "', '"), "'"))
    }
  }
  
  # Generate Plots:
  plotlist <- list()
  for(f in 1:length(factor)){
    factor.feat <- features[[response.name]][[factor[f]]]
    df <- dplyr::bind_rows(factor.feat, .id = "omic")
    df <- dplyr::mutate(df, direction = sapply(df$coefficients, function(coeff){
      if(coeff >= coefficient.cutoff){
        return("positive")
      } else if(coeff <= -coefficient.cutoff){
        return("negative")
      } else {
        return("insignificant")
      }
    }))
    
    numpos <- nrow(dplyr::filter(df, direction == "positive"))
    numneg <- nrow(dplyr::filter(df, direction == "negative"))
    numinsig <- nrow(dplyr::filter(df, direction == "insignificant"))
    
    xmax <- max(abs(df$coefficients))
    
    df <- arrange(df, -abs(coefficients))
    
    # Coefficient cut-off:
    if(coefficient.cutoff > 0){
      before.coeff.cutoff <- nrow(df)
      df <- filter(df, abs(coefficients) >= coefficient.cutoff)
      if(nrow(df) == 0){
        stop("ERROR: no features found above coefficient.cutoff of ", coefficient.cutoff, " for ", factor[f], ".")
      }
      cat("*** Using 'coefficient.cutoff' of ", coefficient.cutoff, ". Removed ", before.coeff.cutoff - nrow(df), " features.\n")
    } else if(!is.null(show.top.cutoff)){
      if(nrow(df) > show.top.cutoff){
        df <- df[1:show.top.cutoff,]
      }
    }
    if(group.by.omic){
      df <- arrange(df, omic, -abs(coefficients))
    }

    # Plot:
    g <- ggplot(df) +
      geom_bar(aes(x = coefficients, y = factor(features, levels = rev(features)), fill = omic), stat = "identity", width = 1, color = "black", size = .2) +
      scale_fill_manual(values = SPEAR.get_color_scheme(SPEARobj, scale.name = "omic")) +
      ylab(NULL) +
      xlab("Coefficient") +
      ggtitle(paste0(factor[f], " | w = ", w)) +
      theme_classic() +
      theme(plot.title = element_text(size = 10)) +
      xlim(c(-xmax, xmax))
    
    if(nrow(df) > 100){
      g <- g + theme(axis.text.y = element_blank())
    }
    
    plotlist[[f]] <- g
  }
  
  p <- cowplot::plot_grid(plotlist = plotlist, nrow = 1)
  
  return(p)
}

#### OLD ______________________________________________________


#' Get factor features
#'@param SPEARobj SPEAR object (returned from run_cv_spear)
#'@param w Weight for SPEAR. Defaults to "best", choosing the best weight per response. Can also be "overall" (choosing the weight with the best overall mean cross-validated error), or one of the weights used to train SPEAR (SPEARobj$params$weights)
#'@param threshold Threshold value of contribution for a factor to be relevant. Defaults to .01.
#'@param cutoff Posterior selection probability cutoff (0 - 1) for a feature to be selected. Defaults to .5
#'@export
SPEAR.get_factor_features <- function(SPEARobj, w = "best", threshold = .01, cutoff = .5){
  # Weight ----------------
  if(w == "best"){
    w.idxs <- SPEAR.get_best_weights(SPEARobj = SPEARobj, w.method = "min", return.overall = FALSE, return.as.index = TRUE)
  } else if(w == "overall"){
    w.idxs <- SPEAR.get_best_weights(SPEARobj = SPEARobj, w.method = "min", return.overall = TRUE, return.as.index = TRUE)
  } else {
    if(!any(SPEARobj$params$weights == w)){
      cat("*** Warning: w = ", w, " not found among possible weights (", paste(SPEARobj$params$weights, collapse = ", "), "). Will use closest SPEAR weight...\n")
      w = SPEARobj$params$weights[which.min(abs(SPEARobj$params$weights-w))]
      w.idxs <- which(SPEARobj$params$weights == w)
      cat(paste0("*** Using SPEAR with w = ", round(SPEARobj$params$weights[w.idxs], 3), "\n"))
      w.idxs <- rep(w.idxs, ncol(SPEARobj$data$Y))
    } else {
      w.idxs <- which(SPEARobj$params$weights == w)
      cat(paste0("*** Using SPEAR with w = ", round(SPEARobj$params$weights[w.idxs], 3), "\n"))
      w.idxs <- rep(w.idxs, ncol(SPEARobj$data$Y))
    }
  } #----------------------
  feature.list <- list()
  for(k in 1:ncol(SPEARobj$data$Y)){
    features.temp <- list()
    relevant.factors <- SPEARobj$cv.eval$factor_contributions[, k, w.idxs[k]] >= threshold
    names(relevant.factors) <- paste0("Factor", 1:length(relevant.factors))
    w.loadings = SPEARobj$fit$results$post_selections[,,w.idxs[k]]
    w.coefficients = t(SPEARobj$fit$results$post_bxs[,,w.idxs[k]])
    w.loadings[,!relevant.factors] <- NA
    w.coefficients[,!relevant.factors] <- NA
    colnames(w.loadings) <- paste0("F", 1:ncol(w.loadings))
    rownames(w.loadings) <- colnames(SPEARobj$data$X)
    colnames(w.coefficients) <- paste0("F", 1:ncol(w.loadings))
    rownames(w.coefficients) <- colnames(SPEARobj$data$X)
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
          w.coefficients.current <- w.coefficients[s:(s + ncol(SPEARobj$data$xlist[[o]]) - 1),factor]
          s <- s + ncol(SPEARobj$data$xlist[[o]])
          w.coefficients.current <- w.coefficients.current[which(w.loadings.current > cutoff)]
          w.loadings.current <- w.loadings.current[which(w.loadings.current > cutoff)]
          # Round loading probabilities to 7 decimals (to rank 1's):
          w.loadings.current <- round(w.loadings.current, 7)
          feat.order <- order(w.loadings.current, abs(w.coefficients.current), decreasing = TRUE)
          w.coefficients.current <- w.coefficients.current[feat.order]
          w.loadings.current <- w.loadings.current[feat.order]
          if(length(w.loadings.current) > 0){
            temp[[names(SPEARobj$data$xlist)[o]]]$features <- names(w.loadings.current)
            temp[[names(SPEARobj$data$xlist)[o]]]$probabilities <- w.loadings.current
            temp[[names(SPEARobj$data$xlist)[o]]]$coefficients <- w.coefficients.current
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



#' Plot CV predictions
#'@param SPEARobj SPEAR object (returned from run_cv_spear)
#'@param w Weight for SPEAR. Defaults to "best", choosing the best weight per response. Can also be "overall" (choosing the weight with the best overall mean cross-validated error), or one of the weights used to train SPEAR (SPEARobj$params$weights)
#'@param ncol Number of columns to display in the final plot. Defaults to NULL (number of different response variables in Y)
#'@param nrow Number of rows to display in the final plot. Defaults to 1
#'@export
SPEAR.plot_cv_predictions <- function(SPEARobj, w = "best", ncol = NULL, nrow = 1){
  plotlist <- list()
  if(w == "best"){
    w.idxs <- apply(SPEARobj$cv.eval$cvm, 2, which.min)
  } else if(w == "overall"){
    w.idxs <- rep(which.min(apply(SPEARobj$cv.eval$cvm,1,sum)), ncol(SPEARobj$data$Y))
  } else {
    if(!any(SPEARobj$params$weights == w)){
      stop(paste0("ERROR: w = ", w, " not found among possible weights (", paste(SPEARobj$params$weights, collapse = ", "), ")."))
    } else {
      w.idxs <- which(SPEARobj$params$weights == w)
      cat(paste0("*** Using SPEAR with w = ", round(SPEARobj$params$weights[w.idxs], 3), "\n"))
      w.idxs <- rep(w.idxs, ncol(SPEARobj$data$Y))
    }
  }
  for(k in 1:ncol(SPEARobj$data$Y)){
    x <- SPEARobj$cv.eval$Yhat.keep[,k,w.idxs[k]] + SPEARobj$cv.eval$intercepts[[1]][w.idxs[k]]
    y <- SPEARobj$data$Y[,k]
    fit <- lm(y ~ x)
    r2 <- summary(fit)$r.squared
    g <- ggplot() +
      geom_point(aes(x = SPEARobj$cv.eval$Yhat.keep[,k,w.idxs[k]] + SPEARobj$cv.eval$intercepts[[1]][w.idxs[k]], y = SPEARobj$data$Y[,k])) + 
      geom_smooth(aes(x = SPEARobj$cv.eval$Yhat.keep[,k,w.idxs[k]] + SPEARobj$cv.eval$intercepts[[1]][w.idxs[k]], y = SPEARobj$data$Y[,k]), color = "#74149E", method = "lm") +
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



#' Plot grid of factor scores per subject
#'@param SPEARobj SPEAR object (returned from run_cv_spear)
#'@param w Weight for SPEAR. Defaults to "overall", choosing the weight with the best overall mean cross-validated error for ALL responses. Can also be one of the weights used to train SPEAR (SPEARobj$params$weights)
#'@param groups A named vector of a grouping variable, where names(groups) = rownames(SPEARobj$data$X). Names are required to match to ensure correct representation of the data.
#'@param factors A vector of factors to be plotted (i.e. c(1, 2)). Defaults to plotting all factors (NULL).
#'@param return.as.list Return plots as a list instead of a cowplot object? Defaults to FALSE
#'@param include.legend Include legend at the bottom of the plot? Defaults to TRUE
#'@export
SPEAR.plot_factor_grid <- function(SPEARobj, w = "overall", groups = NULL, factors = NULL, return.as.list = FALSE, include.legend = TRUE){
    if(is.null(factors)){
      factors <- 1:SPEARobj$params$num_factors
    }
    if(w == "overall"){
      w.idxs <- rep(which.min(apply(SPEARobj$cv.eval$cvm,1,sum)), ncol(SPEARobj$data$Y))
    } else {
      if(!any(SPEARobj$params$weights == w)){
        stop(paste0("ERROR: w = ", w, " not found among possible weights (", paste(SPEARobj$params$weights, collapse = ", "), ")."))
      } else {
        w.idxs <- which(SPEARobj$params$weights == w)
        cat(paste0("*** Using SPEAR with w = ", round(SPEARobj$params$weights[w.idxs], 3), "\n"))
        w.idxs <- rep(w.idxs, ncol(SPEARobj$data$Y))
      }
    }
    # Only use 1 w.idx for grid
    w.idx <- w.idxs[1]
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
              scale_color_distiller(palette = "Spectral", guide = FALSE)
          } else {
            g <- g + geom_density(aes_string(x = paste0("Factor",factors[i]), group = "Group", color = "Group", fill = "Group"), alpha = .15) +
              scale_color_brewer(palette = "Spectral", guide = FALSE) +
              scale_fill_brewer(palette = "Spectral", guide = FALSE)
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
            g <- g + scale_color_distiller(palette = "Spectral", guide = FALSE)
          } else {
            g <- g + scale_color_brewer(palette = "Spectral", guide = FALSE)
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
        scale_fill_brewer(palette = "Spectral", guide = FALSE) +
        theme_bw() +
        theme(legend.key.size = unit(1, 'cm'), legend.position = "top")
      if(any(is.numeric(Uhat$Group))){
        g <- g + scale_color_distiller(palette = "Spectral") +
          guides(color = guide_colourbar(barheight = 0.5))
      } else {
        g <- g + scale_color_brewer(palette = "Spectral")
      }
      legend <- get_legend(g)
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



#' Plot factor scores per subject
#'@param SPEARobj SPEAR object (returned from run_cv_spear)
#'@param w Weight for SPEAR. Defaults to "overall", choosing the weight with the best overall mean cross-validated error for ALL responses. Can also be one of the weights used to train SPEAR (SPEARobj$params$weights)
#'@param groups A named vector of a grouping variable, where names(groups) = rownames(SPEARobj$data$X). Names are required to match to ensure correct representation of the data.
#'@param factors A vector of factors to be plotted (i.e. c(1, 2)). Defaults to plotting all factors (NULL).
#'@param return.as.list Return plots as a list instead of a cowplot object? Defaults to FALSE
#'@param include.legend Include legend at the bottom of the plot? Defaults to TRUE
#'@export
SPEAR.plot_factor_scores <- function(SPEARobj, w = "overall", groups = NULL, factors = NULL, return.as.list = FALSE, include.legend = TRUE){
  if(is.null(factors)){
    factors <- 1:SPEARobj$params$num_factors
  }
  if(w == "overall"){
    w.idxs <- rep(which.min(apply(SPEARobj$cv.eval$cvm,1,sum)), ncol(SPEARobj$data$Y))
  } else {
    if(!any(SPEARobj$params$weights == w)){
      stop(paste0("ERROR: w = ", w, " not found among possible weights (", paste(SPEARobj$params$weights, collapse = ", "), ")."))
    } else {
      w.idxs <- which(SPEARobj$params$weights == w)
      cat(paste0("*** Using SPEAR with w = ", round(SPEARobj$params$weights[w.idxs], 3), "\n"))
      w.idxs <- rep(w.idxs, ncol(SPEARobj$data$Y))
    }
  }
  # Only use 1 w.idx for grid
  w.idx <- w.idxs[1]
  
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
      scale_color_brewer(palette = "Spectral", direction = -1) +
      xlab(NULL) +
      ylab("Factor Score") +
      theme_bw() +
      facet_wrap(vars(variable), nrow = 1)
  } else {
    # Numeric, plot with regression
    g <- ggplot(filter(Uhat.melt, variable %in% paste0("Factor", factors))) +
      geom_point(aes(x = Group, y = value, fill = Group), pch = 21, lwd = 1.5) +
      geom_smooth(aes(x = Group, y = value), color = "black", lwd = .6, alpha = .5, method = "lm") +
      scale_fill_distiller(palette = "Spectral") +
      xlab(NULL) +
      ylab("Factor Score") +
      theme_bw() +
      facet_wrap(vars(variable), nrow = 1)
  }

  return(g)
  
}



#' Plot factor scores per subject
#'@param SPEARobj SPEAR object (returned from run_cv_spear)
#'@param show.w.labels Label points with their weights? Defaults to TRUE
#'@param show.min.w.line Draw a dashed line at the lowest mean cv error? Defaults to TRUE
#'@param show.sd Plot cv error standard deviation? Defaults to TRUE
#'@param show.overall Plot the overall mean cv error (only works with more than one response)? Defaults to FALSE
#'@param plot.per.response Make a different plot for each response? Defaults to TRUE
#'@export
SPEAR.plot_cv_prediction_errors <- function(SPEARobj, show.w.labels = TRUE, show.min.w.line = TRUE, show.sd = TRUE, show.overall = FALSE, plot.per.response = TRUE){
  # Only show one line if there's only one response
  
  # Get the cv errors:
  cv.errors <- as.data.frame(SPEARobj$cv.eval$cvm)
  if(show.overall & ncol(SPEARobj$data$Y) > 1){
    overall.values <- rowSums(cv.errors)/ncol(cv.errors)
  }
  minimum.values <- apply(cv.errors, 2, min)
  cv.errors <- cbind(cv.errors, SPEARobj$params$weights)
  colnames(cv.errors) <- c(colnames(SPEARobj$data$Y), "w")

  # Melt to make graph
  cv.errors.melted <- reshape2::melt(cv.errors, id.vars = c("w"))
  cv.errors.melted$min <- rep(minimum.values, each = nrow(cv.errors))
  colnames(cv.errors.melted) <- c("w", "response", "cvm", "min")
  if(show.sd){
    # Get the cv errors:
    cv.errors.sd <- as.data.frame(SPEARobj$cv.eval$cvsd)
    if(show.overall & ncol(SPEARobj$data$Y) > 1){
      overall.values.sd <- rowSums(cv.errors.sd)/ncol(cv.errors.sd)
    }
    minimum.values.sd <- sapply(1:ncol(cv.errors.sd), function(idx){return(cv.errors.sd[apply(cv.errors, 2, which.min)[idx],idx])})
    cv.errors.sd <- cbind(cv.errors.sd, SPEARobj$params$weights)
    colnames(cv.errors.sd) <- c(colnames(SPEARobj$data$Y), "w")

    # Melt to make graph
    cv.errors.sd.melted <- reshape2::melt(cv.errors.sd, id.vars = c("w"))
    cv.errors.sd.melted$min <- rep(minimum.values.sd, each = nrow(cv.errors.sd))
    colnames(cv.errors.sd.melted) <- c("w", "response", "cvsd", "min.sd")
    
    cv.errors.melted <- dplyr::inner_join(cv.errors.melted, cv.errors.sd.melted, by = c("w", "response"))
  }
  
  # Make the plot:
  g <- ggplot(cv.errors.melted)
  
  if(show.min.w.line){
    g <- g + geom_hline(aes(yintercept = min, color = response), lwd = .5, linetype = "dashed")
    if(exists("overall.values")){
      g <- g + geom_hline(yintercept = min(overall.values), lwd = .5, linetype = "dashed")
    }
  }
  
  if(show.sd){
    g <- g + geom_errorbar(aes(x = w, ymin = cvm - cvsd, ymax = cvm + cvsd, group = response, color = response), size = .4)
    if(show.min.w.line){
      g <- g + geom_hline(aes(yintercept = min + min.sd, color = response), lwd = .2, linetype = "dashed")
      #g <- g + geom_hline(aes(yintercept = min - min.sd, color = response), lwd = .2, linetype = "dashed")
    }
  }
  
  g <- g + geom_line(aes(x = w, y = cvm, group = response, color = response)) +
    geom_point(aes(x = w, y = cvm, group = response, fill = response), size = 3, shape = 21)

  
  g <- g + scale_color_brewer(palette = "Spectral") +
    scale_fill_brewer(palette = "Spectral") +
    scale_x_continuous(breaks = round(SPEARobj$params$weights, 2)) +
    ylab("Mean CV Error") +
    ylim(c(0, NA)) +
    ggtitle("Mean CV Errors of SPEAR weights") +
    theme_bw()
  
  if(show.w.labels){
    g <- g + geom_text(aes(x = w, y = cvm-max(cvm)/10, label = paste0("w=", round(w, 2))), angle = 0)
  }
  if(exists("overall.values")){
    g <- g + geom_line(data = data.frame(w = SPEARobj$params$weights, overall = overall.values), aes(x = w, y = overall), color = "black") +
      geom_point(data = data.frame(w = SPEARobj$params$weights, overall = overall.values), aes(x = w, y = overall), color = "black", size = 3)
  }
  
  if(plot.per.response){
    g <- g + facet_wrap(vars(response))
  }
  
  return(g)
}
  

#' Get ordinal misclassification rate
#'@param SPEARobj SPEAR object (returned from run_cv_spear)
#'@param X List of omics matrices to be used for prediction. Defaults to NULL (use training data stored within SPEAR)
#'@param Y Responses for subjects (rows) in X. Defaults to NULL (use training data stored within SPEAR)
#'@param w Weight for SPEAR. Defaults to "best", choosing the best weight per response. Can also be "overall" (choosing the weight with the best overall mean cross-validated error), or one of the weights used to train SPEAR (SPEARobj$params$weights)
#'@param scale.x Should X be scaled (only used if X is supplied). Defaults to FALSE
#'@export
SPEAR.get_ordinal_misclassification <- function(SPEARobj, X = NULL, Y = NULL, w = "overall", scale.x = FALSE){
  if(!is.null(Y)){
    if(nrow(Y) != nrow(X[[1]])){
      stop(paste0("ERROR: Y provided has ", nrow(Y), " rows, and X has ", nrow(X[[1]]), " rows (need to match)."))
    }
  } else {
    Y <- SPEARobj$data$Y
  }
  preds <- SPEAR.predict_ordinal_classes(SPEARobj, X, w, FALSE, scale.x)
  misclassification.rate <- 1 - (sum(Y == preds)/length(preds))
  return(misclassification.rate)
}



#' Plot factor features
#'@param SPEARobj SPEAR object (returned from run_cv_spear)
#'@param w Weight for SPEAR. Defaults to "best", choosing the best weight per response. Can also be "overall" (choosing the weight with the best overall mean cross-validated error), or one of the weights used to train SPEAR (SPEARobj$params$weights)
#'@param factor Which factors to plot? Defaults to NULL (all relevant factors).
#'@param response.name Which response variable? Check colnames(SPEARobj$data$Y) for options. Defaults to 1st column (NULL)
#'@param threshold Threshold value of contribution for a factor to be relevant. Defaults to .01.
#'@param probability.cutoff Posterior selection probability cutoff (0 - 1) for a feature to be selected. Defaults to .5
#'@param coefficient.cutoff Coefficient cutoff (positive) for a feature to be selected. Will use +/-. Defaults to .01
#'@export
SPEAR.plot_factor_coefficients <- function(SPEARobj, w = "best", factor = NULL, response.name = NULL, threshold = .01, probability.cutoff = .5, coefficient.cutoff = 0){
  # response.name: ------------
  if(is.null(response.name)){
    if(ncol(SPEARobj$data$Y) == 1){
      response.name <- colnames(SPEARobj$data$Y)
    } else {
      cat("\nParameter 'response.name' not provided. Using first response variable ('", colnames(SPEARobj$data$Y)[1], "')\n")
      response.name <- colnames(SPEARobj$data$Y)[1]
    }
  } else if(!any(response.name %in% colnames(SPEARobj$data$Y))){
    stop(paste0("ERROR: 'response.name' provided ('", response.name, "') is not found among the available response variables.\nRequested:\t'", response.name, "'\nAvailable:\t'", paste(colnames(SPEARobj$data$Y), collapse = "', '"), "'"))
  }
  # --------------------------
  
  features <- SPEAR.get_factor_features(SPEARobj, w, threshold, probability.cutoff)
  
  # make sure there are factors that contribute:
  if(length(features[[response.name]]) == 1){
    stop("ERROR: No contributing factors found for 'response.name' ", response.name, ".")
  }
  
  # For each factor...
  if(is.null(factor)){
    factor <- names(features[[response.name]])
  } else {
    factor <- paste0("Factor", factor)
    if(any(!factor %in% names(features[[response.name]]))){
      stop(paste0("ERROR: one or more factors requested not relevant or are missing for Y = '", response.name, "'.\nRequested:\t'", paste(factor, collapse = "', '"), "'\nAvailable:\t'", paste(names(features[[response.name]]), collapse = "', '"), "'"))
    }
  }
  
  # Generate Plots:
  plotlist <- list()
  for(f in 1:length(factor)){
    factor.feat <- features[[response.name]][[factor[f]]]
    df <- dplyr::bind_rows(factor.feat, .id = "omic")
    df <- dplyr::mutate(df, direction = sapply(df$coefficients, function(coeff){
      if(coeff >= coefficient.cutoff){
        return("positive")
      } else if(coeff <= -coefficient.cutoff){
        return("negative")
      } else {
        return("insignificant")
      }
    }))
    
    numpos <- nrow(dplyr::filter(df, direction == "positive"))
    numneg <- nrow(dplyr::filter(df, direction == "negative"))
    numinsig <- nrow(dplyr::filter(df, direction == "insignificant"))

    # Plot:
    g <- ggplot(df) +
      annotate("rect", xmin = coefficient.cutoff, xmax = 100, ymin = -100, ymax = 100, alpha=0.4, fill="#CEFFD3") +
      annotate("rect", xmin = -100, xmax = -coefficient.cutoff, ymin = -100, ymax = 100, alpha=0.4, fill="#FFE3E3") +
      geom_vline(xintercept=coefficient.cutoff, color = "#464646") +
      geom_vline(xintercept=-coefficient.cutoff, color = "#464646") +
      geom_point(aes(x = coefficients, y = probabilities, fill = direction, color = direction), shape = 21, size = 2) +
      coord_cartesian(xlim = c(-max(abs(df$coefficients)), max(abs(df$coefficients))), ylim = c(probability.cutoff,1)) +
      scale_fill_manual(values = c("negative" = "red", "positive" = "green", "insignificant" = "#F3F3F3"), guide = FALSE) +
      scale_color_manual(values = c("negative" = "black", "positive" = "black", "insignificant" = "#464646"), guide = FALSE) +
      xlab("Coefficient") +
      ylab("Probability") +
      ggtitle(paste0(factor[f], " | w = ", w)) +
      theme_bw() +
      theme(title = element_text(size = 6),
            axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10))
    
    plotlist[[f]] <- g
  }
  
  p <- cowplot::plot_grid(plotlist = plotlist, nrow = 1)
  
  return(p)
}


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
