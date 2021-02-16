### SPEARkit Functions ###

#' Predict response values for new samples.
#'@param SPEARobj SPEAR object (returned from run_cv_spear)
#'@export
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


#' Get best SPEAR weights per response Y
#'@param SPEARobj SPEAR object (returned from run_cv_spear)
#'@export
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

#' Get factor contributions to Y
#'@param SPEARobj SPEAR object (returned from run_cv_spear)
#'@export
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


#' Plot factor contributions to Y
#'@param SPEARobj SPEAR object (returned from run_cv_spear)
#'@export
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


#' Get CV prediction errors
#'@param SPEARobj SPEAR object (returned from run_cv_spear)
#'@export
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


#' Plot factor loadings
#'@param SPEARobj SPEAR object (returned from run_cv_spear)
#'@export
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


#' Get factor features
#'@param SPEARobj SPEAR object (returned from run_cv_spear)
#'@export
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




#' Get CV predictions
#'@param SPEARobj SPEAR object (returned from run_cv_spear)
#'@export
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



#' Plot CV predictions
#'@param SPEARobj SPEAR object (returned from run_cv_spear)
#'@export
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



#' Get factor scores
#'@param SPEARobj SPEAR object (returned from run_cv_spear)
#'@export
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




#' Plot factor grid
#'@param SPEARobj SPEAR object (returned from run_cv_spear)
#'@export
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



#' Plot factor scores
#'@param SPEARobj SPEAR object (returned from run_cv_spear)
#'@export
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


