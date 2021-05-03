#' Color text for output in the terminal
#'@param text A string to be colored
#'@param fg Foreground color
#'@param bg Background color
#'@return Text that has been formatted.
#' @examples
#' SPEAR.color_text("Error", fg = "red", bg = NULL)
#'@export
SPEAR.color_text <- function(text, fg = "black", bg = NULL) {
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
      min.ws.idx <- which.min(rowSums(SPEARobj$cv.eval$cvm))
      min.ws.sum <- sum(SPEARobj$cv.eval$cvm[min.ws.idx,])
      min.ws.sd <- sum(SPEARobj$cv.eval$cvsd[min.ws.idx,])
      min.plus.sd.ws <- min.ws.sum + min.ws.sd
      best.idx <- min.ws.idx
      best.w <- SPEARobj$params$weights[min.ws.idx]
      for(t in 1:nrow(SPEARobj$cv.eval$cvm)){
        if(sum(SPEARobj$cv.eval$cvm[t,]) < min.plus.sd.ws & SPEARobj$params$weights[t] > best.w){
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



#' Get a color scheme for a SPEAR plot
#'@param SPEARmodel SPEAR object (returned from get_SPEAR_model)
#'@export
SPEAR.get_color_scheme <- function(SPEARmodel){
  num.omics <- length(SPEARmodel$data$xlist)
  num.responses <- ncol(SPEARmodel$data$Y)
  # Omic
  cs <- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")
  if(num.omics > length(cs)){
    warning("*** Warning: more than ", length(cs), " omics provided. Repeating colors in colorscale.\n")
    num.index <- length(cs)
  } else {
    num.index <- num.omics
  }
  colors.list <- list(
    cs[c(2)],
    cs[c(2, 10)],
    cs[c(2, 4, 10)],
    cs[c(2, 4, 9, 10)],
    cs[c(2, 4, 7, 9, 10)],
    cs[c(1, 2, 4, 7, 9, 10)],
    cs[c(1, 2, 4, 7, 9, 10, 11)],
    cs[c(1, 2, 3, 4, 7, 9, 10, 11)],
    cs[c(1, 2, 3, 4, 7, 8, 9, 10, 11)],
    cs[c(1, 2, 3, 4, 6, 7, 8, 9, 10, 11)],
    cs[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)]
  )
  colors <- colors.list[[num.index]]
  colorvec.omic <- rep(colors, length.out = num.omics)
  names(colorvec.omic) <- names(SPEARmodel$data$xlist)
  
  # Response
  cs <- c("#FC8D62", "#66C2A5", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
  if(num.responses > length(cs)){
    warning("*** Warning: more than ", length(cs), " responses provided. Repeating colors in colorscale.\n")
    num.index <- length(cs)
  } else {
    num.index <- num.responses
  }
  colors.list <- list(
    cs[c(2)],
    cs[c(2, 4)],
    cs[c(2, 4, 1)],
    cs[c(2, 4, 3, 1)],
    cs[c(2, 4, 3, 1, 6)],
    cs[c(2, 4, 3, 1, 6, 5)],
    cs[c(2, 4, 3, 1, 6, 5, 7)],
    cs[c(2, 4, 3, 1, 6, 5, 7, 8)]
  )
  colors <- colors.list[[num.index]]
  colorvec.response <- rep(colors, length.out = num.responses)
  names(colorvec.response) <- colnames(SPEARmodel$data$Y)
  return(c(colorvec.omic, colorvec.response))
}



#' Plot cv prediction loss
#'@param SPEARobj SPEAR object (returned from run_cv_spear)
#'@param show.w.labels Label points with their weights? Defaults to TRUE
#'@param show.min.w.line Draw a dashed line at the lowest mean cv error? Defaults to TRUE
#'@param show.sd Plot cv error standard deviation? Defaults to TRUE
#'@param show.overall Plot the overall mean cv error (only works with more than one response)? Defaults to FALSE
#'@param plot.per.response Make a different plot for each response? Defaults to TRUE
#'@export
SPEAR.plot_cv_loss <- function(SPEARobj, show.w.labels = TRUE, show.min.w.line = TRUE, show.sd = TRUE, show.overall = FALSE, plot.per.response = TRUE){
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
  g <- ggplot2::ggplot(cv.errors.melted)
  
  if(show.min.w.line){
    g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept = min, color = response), lwd = .5, linetype = "dashed")
    if(exists("overall.values")){
      g <- g + ggplot2::geom_hline(yintercept = min(overall.values), lwd = .5, linetype = "dashed")
    }
  }
  
  if(show.sd){
    g <- g + ggplot2::geom_errorbar(ggplot2::aes(x = w, ymin = cvm - cvsd, ymax = cvm + cvsd, group = response, color = response), size = .4)
    if(show.min.w.line){
      g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept = min + min.sd, color = response), lwd = .2, linetype = "dashed")
      #g <- g + geom_hline(aes(yintercept = min - min.sd, color = response), lwd = .2, linetype = "dashed")
    }
  }
  
  g <- g + ggplot2::geom_line(ggplot2::aes(x = w, y = cvm, group = response, color = response)) +
    ggplot2::geom_point(ggplot2::aes(x = w, y = cvm, group = response, fill = response), size = 3, shape = 21)
  
  
  g <- g + ggplot2::scale_color_manual(values = SPEAR.get_color_scheme(SPEARobj)) +
    ggplot2::scale_fill_manual(values = SPEAR.get_color_scheme(SPEARobj)) +
    ggplot2::scale_x_continuous(breaks = round(SPEARobj$params$weights, 2)) +
    ggplot2::ylab("Mean CV Error") +
    ggplot2::ylim(c(0, NA)) +
    ggplot2::ggtitle("Mean CV Errors of SPEAR weights") +
    ggplot2::theme_bw()
  
  if(show.w.labels){
    g <- g + ggplot2::geom_text(ggplot2::aes(x = w, y = cvm-max(cvm)/10, label = paste0("w=", round(w, 2))), angle = 0)
  }
  if(exists("overall.values")){
    g <- g + geom_line(data = data.frame(w = SPEARobj$params$weights, overall = overall.values), ggplot2::aes(x = w, y = overall), color = "black") +
      ggplot2::geom_point(data = data.frame(w = SPEARobj$params$weights, overall = overall.values), ggplot2::aes(x = w, y = overall), color = "black", size = 3)
  }
  
  if(plot.per.response){
    g <- g + ggplot2::facet_wrap(ggplot2::vars(response))
  }
  
  return(g)
}



#' Get a SPEAR model for a specific weight 'w'. The value of 'w' can be specified (i.e. w = 1, w = 0) or the "overall"/"best" weight from 
#'@param SPEARobj SPEAR object (returned from \code{run_cv_spear})
#'@param w.method How to choose best weight? Options include \code{"min"} (lowest mean CV error, default) and \code{"sd"} (choose a higher weight within 1 standard deviation of the CV errors).
#'@return A named vector of weights (or weight indices from \code{SPEARobj$params$weights} if \code{return.as.index} == \code{TRUE}).
#' @examples
#' SPEAR.get_best_weights(SPEARobj, w.method = "min", return.overall = TRUE)
#' SPEAR.get_best_weights(SPEARobj, w.method = "sd")
#'@export
get_SPEAR_model <- function(SPEARobj, w = "overall", w.method = "sd"){
  
  # Weight ----------------
  if(w == "overall" | w == "best"){
    w.idx <- SPEAR.get_best_weights(SPEARobj = SPEARobj, w.method = w.method, return.overall = TRUE, return.as.index = TRUE)[1]
  } else {
    if(!any(SPEARobj$params$weights == w)){
      warning("*** Warning: w = ", w, " not found among possible weights (", paste(SPEARobj$params$weights, collapse = ", "), "). Will use closest SPEAR weight...\n")
      w = SPEARobj$params$weights[which.min(abs(SPEARobj$params$weights-w))]
      w.idx <- which(SPEARobj$params$weights == w)
      cat(paste0("*** Using SPEAR with w = ", round(SPEARobj$params$weights[w.idx], 3), "\n"))
    } else {
      w.idx <- which(SPEARobj$params$weights == w)
      cat(paste0("*** Generating SPEARmodel for w = ", round(SPEARobj$params$weights[w.idx], 3), "\n"))
    }
  } #----------------------
  
  ########
  ### Get SPEAR model for w at w.idx:
  ########
  
  ## factors ----------------------------------------------------------------------
  factors <- list()
  
  # factor.scores 
  W.to_X <- t(SPEARobj$fit$results$post_bxs[,,w.idx])
  W.to_U <- SPEARobj$fit$results$post_betas[,,,w.idx]
  # Get out-of-sample factor scores:
  U.hat.cv <- matrix(0, ncol = SPEARobj$params$num_factors, nrow = nrow(SPEARobj$data$X))
  for(i in 1:nrow(SPEARobj$data$X)){
    U.hat.cv[i,] <- SPEARobj$data$X[i,] %*% SPEARobj$fit$factors_coefs[,,,SPEARobj$params$foldid[i], w.idx] 
  }
  U.hat.train <- SPEARobj$data$X %*% W.to_U
  rownames(U.hat.train) <- rownames(SPEARobj$data$Y)
  rownames(U.hat.cv) <- rownames(SPEARobj$data$Y)
  colnames(U.hat.train) <- paste0("Factor", 1:ncol(U.hat.train))
  colnames(U.hat.cv) <- paste0("Factor", 1:ncol(U.hat.cv))
  factors[['factor.scores']] <- list(in.sample = U.hat.train, out.of.sample = U.hat.cv)
  
  # contributions 
  cat("*** Calculating factor contributions (this may take a second...)\n")
  # X:
  W.to_X.list <- list()
  ind <- 1
  for(d in 1:length(SPEARobj$data$xlist)){
    W.to_X.list[[d]] <- W.to_X[ind:(ind - 1 + ncol(SPEARobj$data$xlist[[d]])),]
    ind <- ind + ncol(SPEARobj$data$xlist[[d]])
  }
  factor_contributions = array(NA,dim = c(SPEARobj$params$num_factors, length(SPEARobj$data$xlist)))
  factor_contributions_pvals = array(NA,dim = c(SPEARobj$params$num_factors, length(SPEARobj$data$xlist)))
  for(k in 1:SPEARobj$params$num_factors){
    for(o in 1:length(SPEARobj$data$xlist)){
      factor.contributions <- c()
      X.tilde <- U.hat.cv[,k] %*% t(W.to_X.list[[o]][,k])
      # Get var.explained per feature:
      for(j in 1:ncol(X.tilde)){
        X.tilde.j <- X.tilde[,j]
        X.j <- SPEARobj$data$xlist[[o]][,j]
        X.tilde.j = X.tilde.j + rnorm(n = length(X.tilde.j), mean = 0, sd = sqrt(var(X.j))*1/nrow(SPEARobj$data$X))
        tmp_pearson = cor(X.j, X.tilde.j)
        suppressWarnings( tmp_pearson <- cor.test(X.j, X.tilde.j, method = 'pearson') ) 
        if( tmp_pearson$estimate<0){
          factor.contributions <- c(factor.contributions, 0)
        }else{
          factor.contributions <- c(factor.contributions, tmp_pearson$estimate^2)
        }
      }
      # Return mean as factor_contributions:
      factor_contributions[k,o] <- mean(factor.contributions)
    }
  }
  colnames(factor_contributions) <- names(SPEARobj$data$xlist)
  rownames(factor_contributions) <- paste0("Factor", 1:nrow(factor_contributions))
  # Y:
  var_explained_Y <- matrix(NA, nrow = dim(SPEARobj$cv.eval$factor_contributions)[1], ncol = dim(SPEARobj$cv.eval$factor_contributions)[2])
  var_explained_Y.pvals <- matrix(NA, nrow = dim(SPEARobj$cv.eval$factor_contributions)[1], ncol = dim(SPEARobj$cv.eval$factor_contributions)[2])
  for(j in 1:dim(SPEARobj$cv.eval$factor_contributions)[2]){
    if(length(SPEARobj$cv.eval$factor_contributions[,j,w.idx]) != 0){
      var_explained_Y[,j] <- SPEARobj$cv.eval$factor_contributions[,j,w.idx]
    }
    if(length(SPEARobj$cv.eval$factor_contributions_pvals[,j,w.idx]) != 0){
      var_explained_Y.pvals[,j] <- SPEARobj$cv.eval$factor_contributions_pvals[,j,w.idx]
    }
  }
  colnames(var_explained_Y) <- colnames(SPEARobj$data$Y)
  rownames(var_explained_Y) <- paste0("Factor", 1:nrow(var_explained_Y))
  colnames(var_explained_Y.pvals) <- colnames(SPEARobj$data$Y)
  rownames(var_explained_Y.pvals) <- paste0("Factor", 1:nrow(var_explained_Y.pvals))
  factors[['contributions']] <- list(X = factor_contributions, Y = var_explained_Y, Y.pvals = var_explained_Y.pvals)
  
  
  #     features ----------------------------------------------------------------------
  cat("*** Getting feature coefficients and posterior selection probabilities for all factors...\n")
  features <- list()
  w.loadings = SPEARobj$fit$results$post_selections[,,w.idx]
  w.coefficients = SPEARobj$fit$results$post_betas[,,1,w.idx]
  colnames(w.loadings) <- paste0("Factor", 1:ncol(w.loadings))
  rownames(w.loadings) <- colnames(SPEARobj$data$X)
  colnames(w.coefficients) <- paste0("Factor", 1:ncol(w.loadings))
  rownames(w.coefficients) <- colnames(SPEARobj$data$X)
  # by factor
  num.omics <- length(SPEARobj$data$xlist)
  factor.features <- list()
  for(f in 1:SPEARobj$params$num_factors){
    temp <- list()
    factor <- paste0("Factor", f)
    s <- 1
    for(o in 1:num.omics){
      omic.features <- list()
      w.loadings.current <- w.loadings[s:(s + ncol(SPEARobj$data$xlist[[o]]) - 1),factor]
      w.coefficients.current <- w.coefficients[s:(s + ncol(SPEARobj$data$xlist[[o]]) - 1),factor]
      s <- s + ncol(SPEARobj$data$xlist[[o]])
      # Round loading probabilities to 7 decimals (to rank 1's):
      w.loadings.current <- round(w.loadings.current, 7)
      feat.order <- order(w.loadings.current, abs(w.coefficients.current), decreasing = TRUE)
      w.coefficients.current <- w.coefficients.current[feat.order]
      w.loadings.current <- w.loadings.current[feat.order]
      temp[[names(SPEARobj$data$xlist)[o]]]$features <- names(w.loadings.current)
      temp[[names(SPEARobj$data$xlist)[o]]]$probabilities <- w.loadings.current
      temp[[names(SPEARobj$data$xlist)[o]]]$coefficients <- w.coefficients.current
    }
    features[[factor]] <- temp
  }
  factors[["features"]] <- features
  
  
  ### Fit ----------------------------------------------------------------------
  fit <- list()
  fit[['post_bxs']] <- array(SPEARobj$fit$results$post_bxs[,,w.idx], dim = dim(SPEARobj$fit$results$post_bxs)[1:2])
  fit[['post_bys']] <- array(SPEARobj$fit$results$post_bys[,,w.idx], dim = dim(SPEARobj$fit$results$post_bys)[1:2])
  fit[['post_betas']] <- array(SPEARobj$fit$results$post_betas[,,,w.idx], dim = dim(SPEARobj$fit$results$post_betas)[1:3])
  fit[['post_betas_cv']] <- array(SPEARobj$fit$factors_coefs[,,,,w.idx], dim = dim(SPEARobj$fit$factors_coefs)[1:4])
  fit[['post_selections']] <- array(SPEARobj$fit$results$post_selections[,,w.idx], dim = dim(SPEARobj$fit$factors_coefs)[1:2])
  fit[['post_bys_cv']] <- array(SPEARobj$fit$projection_coefs[,,,w.idx], dim = dim(SPEARobj$fit$projection_coefs)[1:3])
  fit[['post_bys_scaled']] <- array(SPEARobj$cv.eval$projection_coefs_scaled[,,w.idx], dim = dim(SPEARobj$cv.eval$projection_coefs_scaled)[1:2])
  fit[['post_bys_cv_scaled']] <- array(SPEARobj$cv.eval$cv.projection_coefs_scaled[,,,w.idx], dim = dim(SPEARobj$cv.eval$cv.projection_coefs_scaled)[1:3])
  fit[['reg_coefs']] <- array(SPEARobj$cv.eval$reg_coefs[,,,w.idx], dim = dim(SPEARobj$cv.eval$reg_coefs)[1:3])
  fit[['intercepts']] <- lapply(SPEARobj$cv.eval$intercepts, function(temp){return(temp[w.idx,])})

  
  ### Predictions ----------------------------------------------------------------------
  cat("*** Getting predictions for training samples...\n")
  predictions <- list()
  # CURRENT FIX: Use reg_coefs for out of sample:, fix back go post_bys when they are scaled
  in.sample.preds <- U.hat.train %*% fit[['post_bys_scaled']]
  colnames(in.sample.preds) <- colnames(SPEARobj$data$Y)
  rownames(in.sample.preds) <- rownames(SPEARobj$data$Y)
  out.of.sample.preds <- matrix(0, ncol = ncol(SPEARobj$data$Y), nrow = nrow(SPEARobj$data$Y))
  for(i in 1:nrow(out.of.sample.preds)){
    out.of.sample.preds[i,] <- U.hat.cv[i,] %*% fit[['post_bys_cv_scaled']][,,SPEARobj$params$foldid[i]] 
  }
  colnames(out.of.sample.preds) <- colnames(SPEARobj$data$Y)
  rownames(out.of.sample.preds) <- rownames(SPEARobj$data$Y)
  if(SPEARobj$params$family != "gaussian"){
    in.sample.res <- list()
    out.of.sample.res <- list()
    for(j in 1:ncol(SPEARobj$data$Y)){
      resp.name <- colnames(SPEARobj$data$Y)[j]
      in.sample.signal <- matrix(in.sample.preds[,j], ncol = 1)
      rownames(in.sample.signal) <- rownames(SPEARobj$data$Y)
      colnames(in.sample.signal) <- c(resp.name)
      intercept = fit[['intercepts']][[j]]
      Pmat0 = matrix(0, ncol = length(intercept), nrow = length(in.sample.signal))
      Pmat = matrix(0, ncol = length(intercept) + 1, nrow = length(in.sample.signal))
      for(k in 1:ncol(Pmat0)){
        Pmat0[,k] = 1.0/(1+exp(-in.sample.signal - intercept[k]))
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
      rownames(class.predictions) <- rownames(SPEARobj$data$Y)
      colnames(class.predictions) <- c(resp.name)
      rownames(Pmat) <- rownames(SPEARobj$data$Y)
      colnames(Pmat) <- paste0("Class", 0:(ncol(Pmat)-1))
      in.sample.res[[resp.name]] <- list(class.predictions = class.predictions,
                                         class.probabilities = Pmat,
                                         signal = in.sample.signal)
      
      
      out.of.sample.signal <- matrix(out.of.sample.preds[,j], ncol=1)
      rownames(out.of.sample.signal) <- rownames(SPEARobj$data$Y)
      colnames(out.of.sample.signal) <- c(resp.name)
      intercept = fit[['intercepts']][[j]]
      Pmat0 = matrix(0, ncol = length(intercept), nrow = length(out.of.sample.signal))
      Pmat = matrix(0, ncol = length(intercept) + 1, nrow = length(out.of.sample.signal))
      for(k in 1:ncol(Pmat0)){
        Pmat0[,k] = 1.0/(1+exp(-out.of.sample.signal - intercept[k]))
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
      class.predictions <- matrix(apply(Pmat, 1, which.max) - 1, ncol = 1)
      rownames(class.predictions) <- rownames(SPEARobj$data$Y)
      colnames(class.predictions) <- c(resp.name)
      rownames(Pmat) <- rownames(SPEARobj$data$Y)
      colnames(Pmat) <- paste0("Class", 0:(ncol(Pmat)-1))
      out.of.sample.res[[resp.name]] <- list(class.predictions = class.predictions,
                                             class.probabilities = Pmat,
                                             signal = out.of.sample.signal)
    }
    if(SPEARobj$params$family == "categorical"){
      class.probabilities <- lapply(colnames(SPEARobj$data$Y), function(j){return(in.sample.res[[j]]$class.probabilities)})
      names(class.probabilities) <- colnames(SPEARobj$data$Y)
      overall.class.predictions <- matrix(0, nrow = nrow(SPEARobj$data$Y), ncol = ncol(SPEARobj$data$Y))
      for(i in 1:nrow(overall.class.predictions)){
        overall.class.predictions[i, which.max(sapply(1:ncol(SPEARobj$data$Y), function(ind){
          return(class.probabilities[[ind]][i,2]) # Return the probability of Class1 for each binomial
        }))] <- 1
      }
      colnames(overall.class.predictions) <- colnames(SPEARobj$data$Y)
      rownames(overall.class.predictions) <- rownames(SPEARobj$data$Y)
      predictions[['in.sample']] <- list(overall.class.predictions = overall.class.predictions,
                                         individual.class.predictions = do.call("cbind", lapply(in.sample.res, function(res){return(res$class.predictions)})),
                                         class.probabilities = class.probabilities,
                                         signal.predictions = do.call("cbind", lapply(in.sample.res, function(res){return(res$signal)}))
      )
      class.probabilities <- lapply(colnames(SPEARobj$data$Y), function(j){return(out.of.sample.res[[j]]$class.probabilities)})
      names(class.probabilities) <- colnames(SPEARobj$data$Y)
      overall.class.predictions <- matrix(0, nrow = nrow(SPEARobj$data$Y), ncol = ncol(SPEARobj$data$Y))
      for(i in 1:nrow(overall.class.predictions)){
        overall.class.predictions[i, which.max(sapply(1:ncol(SPEARobj$data$Y), function(ind){
          return(class.probabilities[[ind]][i,2]) # Return the probability of Class1 for each binomial
        }))] <- 1
      }
      colnames(overall.class.predictions) <- colnames(SPEARobj$data$Y)
      rownames(overall.class.predictions) <- rownames(SPEARobj$data$Y)
      predictions[['out.of.sample']] <- list(overall.class.predictions = overall.class.predictions,
                                             individual.class.predictions = do.call("cbind", lapply(out.of.sample.res, function(res){return(res$class.predictions)})),
                                             class.probabilities = class.probabilities,
                                             signal.predictions = do.call("cbind", lapply(out.of.sample.res, function(res){return(res$signal)}))
      )  
    } else { # Family = "binomial" or "ordinal"
      class.probabilities <- lapply(colnames(SPEARobj$data$Y), function(j){return(in.sample.res[[j]]$class.probabilities)})
      names(class.probabilities) <- colnames(SPEARobj$data$Y)
      predictions[['in.sample']] <- list(class.predictions = do.call("cbind", lapply(in.sample.res, function(res){return(res$class.predictions)})),
                                         class.probabilities = class.probabilities,
                                         signal.predictions = do.call("cbind", lapply(in.sample.res, function(res){return(res$signal)}))
                                         )
      class.probabilities <- lapply(colnames(SPEARobj$data$Y), function(j){return(out.of.sample.res[[j]]$class.probabilities)})
      names(class.probabilities) <- colnames(SPEARobj$data$Y)
      predictions[['out.of.sample']] <- list(class.predictions = do.call("cbind", lapply(out.of.sample.res, function(res){return(res$class.predictions)})),
                                         class.probabilities = class.probabilities,
                                         signal.predictions = do.call("cbind", lapply(out.of.sample.res, function(res){return(res$signal)}))
      )
    }
  } else { # gaussian
    predictions[['in.sample']] <- in.sample.preds
    predictions[['out.of.sample']] <- out.of.sample.preds
  }
  
  
  
  ### Params ----------------------------------------------------------------------
  params <- SPEARobj$params
  params[['w']] <- SPEARobj$params$weights[w.idx]
  params[['weights']] <- NULL
  
  
  
  
  # Assemble model:
  SPEARmodel <- list(factors = factors,
                     predictions = predictions,
                     fit = fit,
                     params = params,
                     data = SPEARobj$data)
  cat(paste0("*** Successfully generated SPEARmodel for w = ", round(SPEARobj$params$weights[w.idx], 3), " ***\n"))
  
  # Return the model
  return(SPEARmodel)
  
}



#' Use a SPEARmodel to predict response for new samples (Xlist)
#'@param SPEARmodel SPEAR model (returned from \code{get_SPEAR_model})
#'@param Xlist A list of matrices of testing data (must match the same dimensions as the data used to train the SPEARmodel (SPEARmodel$data$xlist)). Defaults to NULL (returning out-of-sample cv predictions for training data).
#'@return A matrix of predictions (out.of.sample)
#' @examples
#' SPEAR.get_predictions(SPEARmodel)
#' SPEAR.get_predictions(SPEARmodel, Xlist = Xlist)
#'@export
SPEAR.get_predictions <- function(SPEARmodel, Xlist = NULL){
  ### Xlist (test data) ----------
  if(is.null(Xlist)){
    cat("*** Xlist not provided. Returning out-of-sample CV predictions...\n*** (you can get these more easily with SPEARmodel$predictions$out.of.sample)\n")
    return(SPEARmodel$predictions$out.of.sample)
  } else {
    # Quickly check that the dimensions in X match:
    if(length(Xlist) != length(SPEARmodel$data$xlist)){
      stop(paste0("ERROR: Object 'X' passed in has ", length(Xlist), " datasets, whereas the training data for this SPEAR object contained ", length(SPEARmodel$data$xlist), " datasets. Incompatible dimensions."))
    }
    for(d in 1:length(Xlist)){
      if(ncol(Xlist[[d]]) != ncol(SPEARmodel$data$xlist[[d]])){
        stop(paste0("ERROR: Incompatible dimensions. New dataset ", d, " has ", ncol(Xlist[[d]]), " features, whereas training dataset ", d, " has ", ncol(SPEARmodel$data$xlist[[d]]), " features. Cannot predict new responses."))
      }
    }
    # Collapse X:
    X <- do.call("cbind", Xlist)
    if(is.null(rownames(X))){
      warning("*** WARNING: rownames for Xlist not provided. Renaming to sample_1 ...sample_", nrow(X), "\n")
      rownames(X) <- paste0("sample_", 1:nrow(X))
    }
  } # ---------------------
  
  preds.list <- list()
  for(j in 1:ncol(SPEARmodel$data$Y)){
    response <- colnames(SPEARmodel$data$Y)[j]
    ### Group loop? [,g,]
    U <- X %*% SPEARmodel$fit$post_betas[,,1]
    preds <- U %*% SPEARmodel$fit$post_bys_scaled[,j]
    ### Add to group loop?
    
    if(SPEARmodel$params$family != "gaussian"){
      resp.name <- colnames(SPEARmodel$data$Y)[j]
      signal <- matrix(preds, ncol = 1)
      rownames(signal) <- rownames(SPEARmodel$data$Y)
      colnames(signal) <- c(resp.name)
      intercept = SPEARmodel$fit$intercepts[[j]]
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
      rownames(class.predictions) <- rownames(SPEARmodel$data$Y)
      colnames(class.predictions) <- c(resp.name)
      rownames(Pmat) <- rownames(SPEARmodel$data$Y)
      colnames(Pmat) <- paste0("Class", 0:(ncol(Pmat)-1))
      if(SPEARmodel$params$family == "categorical"){
        overall.class.predictions <- matrix(0, nrow = nrow(SPEARmodel$data$Y), ncol = ncol(SPEARmodel$data$Y))
        for(i in 1:nrow(overall.class.predictions)){
          overall.class.predictions[i, which.max(sapply(1:ncol(SPEARmodel$data$Y), function(ind){
            return(Pmat[i,2]) # Return the probability of Class1 for each binomial
          }))] <- 1
        }
        colnames(overall.class.predictions) <- colnames(SPEARmodel$data$Y)
        rownames(overall.class.predictions) <- rownames(SPEARmodel$data$Y)
        preds.list[[response]] <- list(
          overall.class.predictions = overall.class.predictions,
          individual.class.predictions = class.predictions,
          class.probabilities = Pmat,
          signal = signal
        )
      } else { # ordinal or binomial
        preds.list[[response]] <- list(
          class.predictions = class.predictions,
          class.probabilities = Pmat,
          signal = signal
        )
      }
    } else { # gaussian
      preds.list[[response]] <- preds
    }
  }
  # Polish results:
  if(SPEARmodel$params$family == "gaussian"){
    preds <- do.call("cbind", preds.list)
    rownames(preds) <- rownames(X)
    colnames(preds) <- colnames(SPEARmodel$data$Y)
    return(preds)
  } else if(SPEARmodel$params$family == "categorical"){
    class.probabilities <- lapply(1:ncol(SPEARmodel$data$Y), function(j){return(preds.list[[j]]$class.probabilities)})
    names(class.probabilities) <- colnames(SPEARmodel$data$Y)
    overall.class.predictions <- matrix(0, nrow = nrow(SPEARmodel$data$Y), ncol = ncol(SPEARmodel$data$Y))
    for(i in 1:nrow(overall.class.predictions)){
      overall.class.predictions[i, which.max(sapply(1:ncol(SPEARmodel$data$Y), function(ind){
        return(class.probabilities[[ind]][i,2]) # Return the probability of Class1 for each binomial
      }))] <- 1
    }
    colnames(overall.class.predictions) <- colnames(SPEARmodel$data$Y)
    rownames(overall.class.predictions) <- rownames(SPEARmodel$data$Y)
    preds <- list(overall.class.predictions = overall.class.predictions,
                  individual.class.predictions = do.call("cbind", lapply(preds.list, function(res){return(res$class.predictions)})),
                  class.probabilities = class.probabilities,
                  signal.predictions = do.call("cbind", lapply(preds.list, function(res){return(res$signal)}))
    )
    return(preds)
  } else { # Family = "binomial" or "ordinal"
    class.probabilities <- lapply(colnames(SPEARmodel$data$Y), function(j){return(preds.list[[j]]$class.probabilities)})
    names(class.probabilities) <- colnames(SPEARmodel$data$Y)
    preds <- list(class.predictions = do.call("cbind", lapply(preds.list, function(res){return(res$class.predictions)})),
                  class.probabilities = class.probabilities,
                  signal.predictions = do.call("cbind", lapply(preds.list, function(res){return(res$signal)}))
    )
    return(preds)
  }
}



#' Use a SPEARmodel to get factor scores for new samples (Xlist)
#'@param SPEARmodel SPEAR model (returned from \code{get_SPEAR_model})
#'@param Xlist A list of matrices of testing data (must match the same dimensions as the data used to train the SPEARmodel (SPEARmodel$data$xlist)). Defaults to NULL (returning out-of-sample cv predictions for training data).
#'@return A matrix of factor scores (out.of.sample)
#' @examples
#' SPEAR.get_factor_scores(SPEARmodel)
#' SPEAR.get_factor_scores(SPEARmodel, Xlist = Xlist)
#'@export
SPEAR.get_factor_scores <- function(SPEARmodel, Xlist = NULL){
  ### Xlist (test data) ----------
  if(is.null(Xlist)){
    cat("*** Xlist not provided. Returning out-of-sample CV factor scores...\n*** (you can get these more easily with SPEARmodel$factors$factor.scores$out.of.sample)\n")
    return(SPEARmodel$factors$factor.scores$out.of.sample)
  } else {
    # Quickly check that the dimensions in X match:
    if(length(Xlist) != length(SPEARmodel$data$xlist)){
      stop(paste0("ERROR: Object 'X' passed in has ", length(Xlist), " datasets, whereas the training data for this SPEAR object contained ", length(SPEARmodel$data$xlist), " datasets. Incompatible dimensions."))
    }
    for(d in 1:length(Xlist)){
      if(ncol(Xlist[[d]]) != ncol(SPEARmodel$data$xlist[[d]])){
        stop(paste0("ERROR: Incompatible dimensions. New dataset ", d, " has ", ncol(Xlist[[d]]), " features, whereas training dataset ", d, " has ", ncol(SPEARmodel$data$xlist[[d]]), " features. Cannot predict new responses."))
      }
    }
    # Collapse X:
    X <- do.call("cbind", Xlist)
    if(is.null(rownames(X))){
      warning("*** WARNING: rownames for Xlist not provided. Renaming to sample_1 ...sample_", nrow(X), "\n")
      rownames(X) <- paste0("sample_", 1:nrow(X))
    }
    # factor.scores 
    # GROUP: Assuming group = 1
    W.to_U <- SPEARmodel$fit$post_betas[,,1]
    U.hat <- X %*% W.to_U
    rownames(U.hat) <- rownames(X)
    colnames(U.hat) <- paste0("Factor", 1:ncol(U.hat))
    return(U.hat)
  }
}



#' Use a SPEARmodel to predict response for new samples (Xlist)
#'@param SPEARmodel SPEAR model (returned from \code{get_SPEAR_model})
#'@param groups A named vector of a grouping variable, where names(groups) = rownames(SPEARmodel$data$X). Names are required to match to ensure correct representation of the data.
#'@return A plot showing in-sample vs. out-of-sample predictions per response
#' @examples
#' SPEAR.get_cv_predictions(SPEARmodel)
#'@export
SPEAR.plot_overfit <- function(SPEARmodel, groups = NULL){
  if(SPEARmodel$params$family == "gaussian"){
    df.list <- list()
    for(j in 1:ncol(SPEARmodel$data$Y)){
      df <- data.frame(in.sample = SPEARmodel$predictions$in.sample[,j],
                       out.of.sample = SPEARmodel$predictions$out.of.sample[,j],
                       group = NaN,
                       response = colnames(SPEARmodel$data$Y)[j])
      rownames(df) <- rownames(SPEARmodel$data$Y)
      if(!is.null(groups)){
        # Check for mapping of metadata to subjects:
        if(is.null(names(groups)) | any(!names(groups) %in% rownames(df))){
          stop("ERROR. 'groups' needs to be a named list to ensure the correct mapping to subjects. Please check that your subject names match.")
        }
        else{
          df$group <- sapply(rownames(df), function(subject){return(groups[which(names(groups) == subject)])})
        }
      }
      df$sample <- rownames(df)
      rownames(df) <- NULL
      # Add to list:
      df.list[[j]] <- df
    }
    df.total <- do.call("rbind", df.list)

    g <- ggplot2::ggplot(df.total) +
         ggplot2::geom_point(ggplot2::aes(x = in.sample, y = out.of.sample, color = group)) +
         ggplot2::xlab("in-sample predictions") +
         ggplot2::ylab("out-of-sample predictions") +
         ggplot2::theme_bw() +
         ggplot2::theme(plot.title = ggplot2::element_text(hjust = .5)) + 
         ggplot2::facet_wrap(ggplot2::vars(response))
  } else {
    df.list <- list()
    for(j in 1:ncol(SPEARmodel$data$Y)){
      df <- data.frame(in.sample = SPEARmodel$predictions$in.sample$signal.predictions[,j],
                       out.of.sample = SPEARmodel$predictions$out.of.sample$signal.predictions[,j],
                       group = NaN,
                       response = colnames(SPEARmodel$data$Y)[j])
      rownames(df) <- rownames(SPEARmodel$data$Y)
      if(!is.null(groups)){
        # Check for mapping of metadata to subjects:
        if(is.null(names(groups)) | any(!names(groups) %in% rownames(df))){
          stop("ERROR. 'groups' needs to be a named list to ensure the correct mapping to subjects. Please check that your subject names match.")
        }
        else{
          df$group <- sapply(rownames(df), function(subject){return(groups[which(names(groups) == subject)])})
        }
      }
      df$sample <- rownames(df)
      rownames(df) <- NULL
      # Add to list:
      df.list[[j]] <- df
    }
    df.total <- do.call("rbind", df.list)
    
    g <- ggplot2::ggplot(df.total) +
      ggplot2::geom_point(ggplot2::aes(x = in.sample, y = out.of.sample, color = group)) +
      ggplot2::xlab("in-sample predictions") +
      ggplot2::ylab("out-of-sample predictions") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = .5)) +
      ggplot2::facet_wrap(ggplot2::vars(response))
  }
  
  # Return plots
  return(g)
}



#' Plot class probabilities for a SPEARmodel (of family "binomial", "categorical", or "ordinal")
#'@param SPEARmodel SPEAR model (returned from \code{get_SPEAR_model})
#'@param groups A named vector of a grouping variable, where names(groups) = rownames(SPEARmodel$data$X). Names are required to match to ensure correct representation of the data.
#'@param forecast Which probabilities to use? A string with "out.of.sample" or "in.sample". Defaults to "out.of.sample".
#'@param Xlist A list of matrices of testing data (must match the same dimensions as the data used to train the SPEARmodel (SPEARmodel$data$xlist)). Defaults to NULL (use training data).
#'@return A plot showing in-sample vs. out-of-sample predictions per response
#' @examples
#' SPEAR.plot_class_probabilities(SPEARmodel)
#' 
#' groups <- SPEARmodel$params$foldid
#' names(groups) <- rownames(SPEARmodel$data$Y)
#' SPEAR.plot_class_probabilities(SPEARmodel, groups = groups, forecast = "in.sample")
#'@export
SPEAR.plot_class_probabilities <- function(SPEARmodel, groups = NULL, forecast = "out.of.sample", Xlist = NULL){
  if(SPEARmodel$params$family == "gaussian"){
    stop('*** ERROR: SPEARmodel must be of type "ordinal", "categorical", or "binomial" to plot class probabilities.')
  } else {
    if(!forecast %in% c("in.sample", "out.of.sample")){
      stop('*** ERROR: unknown forecast parameter "', forecast, '". Can be "in.sample" or "out.of.sample"')
    }
    ### Xlist (test data) ----------
    if(is.null(Xlist)){
      class.probabilities <- SPEARmodel$predictions[[forecast]]
    } else {
      class.probabilities <- SPEAR::SPEAR.get_predictions(SPEARmodel, Xlist)
    } # ------
    
    df.list <- list()
    for(j in 1:length(class.probabilities$class.probabilities)){
      df <- as.data.frame(class.probabilities$class.probabilities[[j]])
      df$group <- NaN
      df$sample <- rownames(df)

      if(!is.null(groups)){
        # Check for mapping of metadata to subjects:
        if(is.null(names(groups)) | any(!names(groups) %in% rownames(df))){
          stop("ERROR. 'groups' needs to be a named list to ensure the correct mapping to subjects. Please check that your subject names match.")
        }
        else{
          df$group <- sapply(rownames(df), function(sample){return(groups[which(names(groups) == sample)])})
        }
      }
      df.melt <- reshape2::melt(df, id.vars = c("group", "sample"))
      df.melt$response = response = colnames(SPEARmodel$data$Y)[j]
      # Add to list:
      df.list[[j]] <- df.melt
    }
    df.total <- do.call("rbind", df.list)
    
    g <- ggplot2::ggplot(df.total) +
      ggplot2::geom_line(ggplot2::aes(x = variable, y = value, color = group, group = sample)) +
      ggplot2::geom_point(ggplot2::aes(x = variable, y = value, color = group)) +
      ggplot2::xlab("Class") +
      ggplot2::ylab("Probability") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = .5)) +
      ggplot2::facet_wrap(ggplot2::vars(response), nrow = 1)
  }
  
  # Return plots
  return(g)
}



#' Plot factor contributions for X and Y
#'@param SPEARmodel SPEAR model (returned from get_SPEAR_model)
#'@param threshold Threshold value of contribution for a factor to be relevant. Defaults to .01.
#'@param show.labels Show the contributions as geom_text on the plot? Defaults to TRUE
#'@param show.irrelevant Show all contributions, even those below the threshold? Defaults to FALSE
#'@export
SPEAR.plot_factor_contributions <- function(SPEARmodel, threshold = .01, show.labels = TRUE, show.irrelevant = FALSE){
  
    w <- SPEARmodel$params$w
    factor.contributions <- SPEARmodel$factors$contributions
    
    if(show.irrelevant){
      title <- paste0("Factor Contributions | ", w)
    } else {
      title <- paste0("Factor Contributions >= ", threshold)
      SPEARmodel$factors$contributions$X[SPEARmodel$factors$contributions$X < threshold] <- 0
      SPEARmodel$factors$contributions$Y[SPEARmodel$factors$contributions$Y < threshold] <- 0
    }
    
    # Combine DFs:
    df.X <- SPEARmodel$factors$contributions$X
    df.Y <- SPEARmodel$factors$contributions$Y
    df.comb <- cbind(df.X, df.Y)
    # Var explained:
    df.var <- data.frame(omic = colnames(df.comb), var = colSums(df.comb), omic.num = 1:ncol(df.comb), label = format(round(ifelse(colSums(df.comb) < threshold, NA, colSums(df.comb)), 2), nsmall = 2))
    # Plot
    df.melt <- reshape2::melt(df.comb)
    df.melt$value.adj <- ifelse(df.melt$value < 0, 0, df.melt$value)
    if(!show.irrelevant){
      df.melt$label <- ifelse(df.melt$value==0, NA, format(round(df.melt$value, 2), nsmall = 2))
    } else {
      df.melt$label <- round(df.melt$value, 2)
    }
    
    g <- ggplot2::ggplot(data = df.melt) +
      ggplot2::geom_tile(ggplot2::aes(y = Var2, x = Var1, alpha = value.adj, fill = Var2), color = "black") +
      ggplot2::geom_text(ggplot2::aes(y = Var2, x = Var1, label = label), color = ifelse(show.labels, "black", NA)) +
      ggplot2::geom_rect(data = df.var, ggplot2::aes(xmin = nrow(df.comb) + .5, xmax = nrow(df.comb) + 1.5, ymin = .5, ymax = nrow(df.var) + .5), color = "white", fill = "white") +
      ggplot2::geom_rect(data = df.var, ggplot2::aes(xmin = nrow(df.comb) + .5, xmax = nrow(df.comb) + .5 + var, ymin = omic.num - .5, ymax = omic.num + .5, fill = omic), color = "black") +
      ggplot2::geom_text(data = df.var, ggplot2::aes(y = omic, x = nrow(df.comb) + 1, label = label), color = "black", angle = 270) +
      ggplot2::geom_hline(yintercept = length(SPEARmodel$data$xlist) + .5, size = 1) +
      ggplot2::scale_alpha_continuous(range = c(0, 1), guide = F) +
      ggplot2::scale_fill_manual(values = SPEAR.get_color_scheme(SPEARmodel), guide = FALSE) +
      ggplot2::ggtitle(title, paste0("w = ", w)) +
      ggplot2::ylab("") +
      ggplot2::xlab("") +
      ggplot2::coord_equal() +
      ggplot2::theme_classic() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10))

  return(g)
}



#' Plot factor scores for a SPEARmodel
#'@param SPEARmodel SPEAR model (returned from \code{get_SPEAR_model})
#'@param groups A named vector of a grouping variable, where names(groups) = rownames(SPEARmodel$data$X). Names are required to match to ensure correct representation of the data.
#'@param forecast Which probabilities to use? A string with "out.of.sample" or "in.sample". Defaults to "out.of.sample".
#'@param Xlist A list of matrices of testing data (must match the same dimensions as the data used to train the SPEARmodel (SPEARmodel$data$xlist)). Defaults to NULL (use training data).
#'@param jitter.points Use geom_jitter vs. geom_point? Defaults to TRUE
#'@return A plot showing in-sample vs. out-of-sample predictions per response
#' @examples
#' SPEAR.plot_factor_scores(SPEARmodel, jitter.points = FALSE)
#' 
#' groups <- SPEARmodel$params$foldid
#' names(groups) <- rownames(SPEARmodel$data$Y)
#' SPEAR.plot_factor_scores(SPEARmodel, groups = groups, forecast = "in.sample")
#'@export
SPEAR.plot_factor_scores <- function(SPEARmodel, groups = NULL, forecast = "out.of.sample", Xlist = NULL, jitter.points = TRUE){
  ### Xlist (test data) ----------
  if(is.null(Xlist)){
    factor.scores <- as.data.frame(SPEARmodel$factors$factor.scores[[forecast]])
  } else {
    factor.scores <- as.data.frame(SPEAR.get_factor_scores(SPEARmodel, Xlist))
  } # ------
  
  factor.scores$group <- NaN
  show.groups <- FALSE
  factor.scores$sample <- rownames(factor.scores)

  if(!is.null(groups)){
    # Check for mapping of metadata to subjects:
    if(is.null(names(groups)) | any(!names(groups) %in% rownames(factor.scores))){
      stop("ERROR. 'groups' needs to be a named list to ensure the correct mapping to subjects. Please check that your subject names match.")
    }
    else{
      factor.scores$group <- sapply(rownames(factor.scores), function(sample){return(groups[which(names(groups) == sample)])})
      show.groups <- TRUE
    }
  }
  factor.scores.melt <- reshape2::melt(factor.scores, id.vars = c("group", "sample"))
  
  g <- ggplot2::ggplot(factor.scores.melt) + ggplot2::theme_bw()
  
  if(show.groups){
    if(jitter.points){
      g <- g + ggplot2::geom_jitter(ggplot2::aes(x = group, y = value, fill = group), color = "grey", shape = 21, size = 1.5, stroke = .1)
    } else {
      g <- g + ggplot2::geom_point(ggplot2::aes(x = group, y = value, fill = group), color = "grey", shape = 21, size = 1.5, stroke = .1)
    }
    if(any(is.numeric(groups))){
      g <- g + ggplot2::geom_smooth(ggplot2::aes(x = group, y = value), method = "lm", color = "black")
    } else {
      g <- g + ggplot2::geom_boxplot(ggplot2::aes(x = group, y = value), alpha = 0)
    }
  } else {
    if(jitter.points){
      g <- g + ggplot2::geom_jitter(ggplot2::aes(x = 1, y = value), fill = "white", color = "black", shape = 21, size = 1.5, stroke = .5) +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank())
    } else {
      g <- g + ggplot2::geom_point(ggplot2::aes(x = 1, y = value), fill = "white", color = "black", shape = 21, size = 1.5, stroke = .5) +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank())
    }
    if(any(is.numeric(groups))){
      g <- g + ggplot2::geom_smooth(ggplot2::aes(x = group, y = value), method = "lm", color = "black")
    } else {
      g <- g + ggplot2::geom_boxplot(ggplot2::aes(x = 1, y = value), alpha = 0)
    }
  }
  
  g <- g + ggplot2::xlab(NULL) +
    ggplot2::ylab(paste0("Factor Score (", forecast, ")")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = .5)) +
    ggplot2::facet_wrap(ggplot2::vars(variable), nrow = 1)
  
  return(g)
}



#' Plot factor scores for a SPEARmodel
#'@param SPEARmodel SPEAR model (returned from \code{get_SPEAR_model})
#'@param groups A named vector of a grouping variable, where names(groups) = rownames(SPEARmodel$data$X). Names are required to match to ensure correct representation of the data.
#'@param forecast Which probabilities to use? A string with "out.of.sample" or "in.sample". Defaults to "out.of.sample".
#'@param Xlist A list of matrices of testing data (must match the same dimensions as the data used to train the SPEARmodel (SPEARmodel$data$xlist)). Defaults to NULL (use training data).
#'@param fit.line Add linear line of best fit? Defaults to TRUE
#'@return A plot showing in-sample vs. out-of-sample predictions per response
#' @examples
#' SPEAR.plot_factor_grid(SPEARmodel, jitter.points = FALSE)
#' 
#' groups <- SPEARmodel$params$foldid
#' names(groups) <- rownames(SPEARmodel$data$Y)
#' SPEAR.plot_factor_grid(SPEARmodel, groups = groups, forecast = "in.sample")
#'@export
SPEAR.plot_factor_grid <- function(SPEARmodel, groups = NULL, forecast = "out.of.sample", Xlist = NULL, fit.line = TRUE){
  ### Xlist (test data) ----------
  if(is.null(Xlist)){
    factor.scores <- as.data.frame(SPEARmodel$factors$factor.scores[[forecast]])
  } else {
    factor.scores <- as.data.frame(SPEAR.get_factor_scores(SPEARmodel, Xlist))
  } # ------
  factor.scores$group <- NaN
  show.groups <- FALSE
  factor.scores$sample <- rownames(factor.scores)
  if(!is.null(groups)){
    # Check for mapping of metadata to subjects:
    if(is.null(names(groups)) | any(!names(groups) %in% rownames(factor.scores))){
      stop("ERROR. 'groups' needs to be a named list to ensure the correct mapping to subjects. Please check that your subject names match.")
    }
    else{
      # Check if there are at least 3 labels for each group:
      if(any(table(groups) <= 2)){
        stop("*** ERROR: 1 or more groups provided in the 'groups' argument have 2 or less entries. Cannot perform correlation analyses. Groups need to have at least 3 members.")
      }
      # Make categorical if is.numeric:
      if(any(is.numeric(groups))){
        
      }
      factor.scores$group <- sapply(rownames(factor.scores), function(sample){return(groups[which(names(groups) == sample)])})
      show.groups <- TRUE
    }
  }
  if(!show.groups){
    g <- GGally::ggpairs(dplyr::select(factor.scores, dplyr::starts_with("Factor"), group)) + ggplot2::theme_bw()
  } else {
    g <- GGally::ggpairs(dplyr::select(factor.scores, dplyr::starts_with("Factor"), group), mapping=ggplot2::aes(colour = group)) + ggplot2::theme_bw()
  }
  return(g)
}



#' Plot factor loadings for a SPEARmodel
#'@param SPEARmodel SPEAR model (returned from \code{get_SPEAR_model})
#'@param plot.per.omic Should a bar with colors for each omic be plotted along side the loadings? Defaults to TRUE
#'@return A plot with the factor loadings of the SPEAR model
#' @examples
#' SPEAR.plot_factor_loadings(SPEARmodel)
#'@export
SPEAR.plot_factor_loadings <- function(SPEARmodel, plot.per.omic = TRUE){
  ### Group loop? [,g,]
  w.coefficients = SPEARmodel$fit$post_betas[,,1]
  ### Add to group loop?
  colnames(w.coefficients) <- paste0("Factor", 1:ncol(w.coefficients))
  rownames(w.coefficients) <- colnames(SPEARmodel$data$X)
  df <- reshape2::melt(w.coefficients)
  colnames(df)[3] <- "coefficient"
  # Add bars:
  omic.ends <- c(0, as.vector(sapply(SPEARmodel$data$xlist, ncol)))
  for(o in 2:length(omic.ends)){
    omic.ends[o] <- omic.ends[o] + omic.ends[o-1]
  }
  omic.df <- data.frame(x.start = .43,
                        x.end = .45,
                        y.start = omic.ends[1:(length(omic.ends)-1)],
                        y.end = omic.ends[2:length(omic.ends)],
                        omic = names(SPEARmodel$data$xlist))

  # Make plot
  g <- ggplot2::ggplot(df) +
    ggplot2::geom_tile(ggplot2::aes(x = Var2, y = factor(Var1, levels = unique(Var1)), fill = coefficient)) +
    ggplot2::ylab(paste0("Features\n(n = ", ncol(SPEARmodel$data$X), ")")) +
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
    ggplot2::scale_color_manual(values = SPEAR.get_color_scheme(SPEARmodel)) +
    ggplot2::xlab(NULL) +
    ggplot2::ggtitle("Factor Loadings", paste0("w = ",round(SPEARmodel$params$w, 3))) +
    ggplot2::theme_void() +
    ggplot2::theme(axis.title.y = ggplot2::element_text(size = 10, angle = 90),
          axis.text.x = ggplot2::element_text(size = 10),
          axis.text.y=ggplot2::element_blank(),
          axis.ticks.y=ggplot2::element_blank(),
          plot.margin = ggplot2::unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          plot.title = ggplot2::element_text(hjust = .5, face = "plain", size = 12),
          plot.subtitle = ggplot2::element_text(hjust = .5, face = "plain", size = 10),
          legend.text = ggplot2::element_text(size = 8),
          legend.title = ggplot2::element_text(size = 10))
  
  if(plot.per.omic){
    g <- g + ggplot2::geom_rect(data = omic.df, ggplot2::aes(xmin = x.start, xmax = x.end, ymin = y.start, ymax = y.end, color = omic), fill = "white", size = 1)
  }
  
  return(g)
}



#' Plot posterior selection probabilities for a SPEARmodel
#'@param SPEARmodel SPEAR model (returned from \code{get_SPEAR_model})
#'@param plot.per.omic Should a bar with colors for each omic be plotted along side the loadings? Defaults to TRUE
#'@return A plot with the factor loadings of the SPEAR model
#' @examples
#' SPEAR.plot_factor_loadings(SPEARmodel)
#'@export
SPEAR.plot_factor_probabilities <- function(SPEARmodel, plot.per.omic = TRUE){
  ### Group loop? [,g,]
  w.probabilities = SPEARmodel$fit$post_selections
  ### Add to group loop?
  colnames(w.probabilities) <- paste0("Factor", 1:ncol(w.probabilities))
  rownames(w.probabilities) <- colnames(SPEARmodel$data$X)
  df <- reshape2::melt(w.probabilities)
  colnames(df)[3] <- "probability"
  # Add bars:
  omic.ends <- c(0, as.vector(sapply(SPEARmodel$data$xlist, ncol)))
  for(o in 2:length(omic.ends)){
    omic.ends[o] <- omic.ends[o] + omic.ends[o-1]
  }
  omic.df <- data.frame(x.start = .43,
                        x.end = .45,
                        y.start = omic.ends[1:(length(omic.ends)-1)],
                        y.end = omic.ends[2:length(omic.ends)],
                        omic = names(SPEARmodel$data$xlist))
  # Make plot
  g <- ggplot2::ggplot(df) +
    ggplot2::geom_tile(ggplot2::aes(x = Var2, y = factor(Var1, levels = unique(Var1)), fill = probability)) +
    ggplot2::ylab(paste0("Features\n(n = ", ncol(SPEARmodel$data$X), ")")) +
    ggplot2::scale_fill_gradient2(mid = "white", high = "#b2182b", limits = c(0, 1)) +
    ggplot2::scale_color_manual(values = SPEAR.get_color_scheme(SPEARmodel)) +
    ggplot2::xlab(NULL) +
    ggplot2::ggtitle("Post Selection Probabilities", paste0("w = ",round(SPEARmodel$params$w, 3))) +
    ggplot2::theme_void() +
    ggplot2::theme(axis.title.y = ggplot2::element_text(size = 10, angle = 90),
                   axis.text.x = ggplot2::element_text(size = 10),
                   axis.text.y=ggplot2::element_blank(),
                   axis.ticks.y=ggplot2::element_blank(),
                   plot.margin = ggplot2::unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
                   plot.title = ggplot2::element_text(hjust = .5, face = "plain", size = 12),
                   plot.subtitle = ggplot2::element_text(hjust = .5, face = "plain", size = 10),
                   legend.text = ggplot2::element_text(size = 8),
                   legend.title = ggplot2::element_text(size = 10))
  
  if(plot.per.omic){
    g <- g + ggplot2::geom_rect(data = omic.df, ggplot2::aes(xmin = x.start, xmax = x.end, ymin = y.start, ymax = y.end, color = omic), fill = "white", size = 1)
  }
  
  return(g)
}



#' Get a data.frame of all features that pass cutoffs from the chosen model:
#'@param SPEARmodel SPEAR model (returned from \code{get_SPEAR_model})
#'@param factors A vector of integers representing factors. If NULL, will return features for all factors. Defaults to NULL (all factors)
#'@param omics A vector of names of datasets (see names(SPEARmodel$data$xlist)). Defaults to NULL (all omics)
#'@param coefficient.cutoff Cutoff for coefficient magnitude. If .1 is provided, abs(coefficients) >= .1 are returned. Defaults to .01
#'@param probability.cutoff Cutoff for posterior selection probability. Ranges from 0 - 1. Defaults to .5
#'@param method Which method for correlation? Defaults to "spearman"
#'@return A data.frame with features that passed cutoffs set by the parameters.
#' @examples
#' SPEAR.get_feature_table(SPEARmodel)
#' 
#' # Get specific features for one omic:
#' omic.name <- names(SPEARmodel$data$X)[1]
#' SPEAR.get_feature_table(SPEARmodel, factors = 1, omics = omic.name)
#' 
#' # Get features that pass cutoffs:
#' SPEAR.get_feature_table(SPEARmodel, coefficient.cutoff = .2, probability.cutoff = .5)
#'@export
SPEAR.get_feature_table <- function(SPEARmodel, factors = NULL, omics = NULL, coefficient.cutoff = .01, probability.cutoff = .5, method = "spearman"){
  # Check that factors is within range:
  if(is.null(factors)){
    factors <- 1:SPEARmodel$params$num_factors
  }
  if(any(!is.numeric(factors))){
    stop("*** ERROR: Factors provided (", factors, ") are not integers. Please choose a vector of integers from 1:", SPEARmodel$params$num_factors)
  } else if(any(!factors %in% 1:SPEARmodel$params$num_factors)){
    stop("*** ERROR: Some/all factors provided (", factors, ") are not within the possible range of 1:", SPEARmodel$params$num_factors)
  }
  # Passed factors:
  factors <- paste0("Factor", factors)
  
  # Check if omics exist:
  if(is.null(omics)){
    omics <- names(SPEARmodel$data$xlist)
  }
  if(any(!omics %in% names(SPEARmodel$data$xlist))){
    stop("*** ERROR: Some/all omics provided (", paste(omics, collapse = " ,"), ") are not found in the possible omics names: ", paste(names(SPEARmodel$data$xlist), collapse = ', '))
  }
  
  # Get relevant features:
  results.list <- list()
  for(f in factors){
    for(o in omics){
      feature.vec <- SPEARmodel$factors$features[[f]][[o]]
      passed.coeff.cutoff <- abs(feature.vec$coefficients) >= coefficient.cutoff
      passed.prob.cutoff <- feature.vec$probabilities >= probability.cutoff
      passed.total.cutoff <- passed.coeff.cutoff & passed.prob.cutoff
      if(sum(passed.total.cutoff) > 0){
        res <- as.data.frame(feature.vec)[passed.coeff.cutoff & passed.prob.cutoff,]
        res$omic <- o
        res$factor <- f
        res$cor.in.sample <- NA
        res$cor.out.of.sample <- NA
        res$pval.in.sample <- NA
        res$pval.out.of.sample <- NA
        # Correlations:
        for(r in 1:nrow(res)){
          vals <- SPEARmodel$data$xlist[[o]][,which(colnames(SPEARmodel$data$xlist[[o]]) == res$features[r])]
          res$cor.in.sample[r] <- cor(vals, SPEARmodel$factors$factor.scores$in.sample[,f], method = method)
          res$cor.out.of.sample[r] <- cor(vals, SPEARmodel$factors$factor.scores$out.of.sample[,f], method = method)
          res$pval.in.sample[r] <- cor.test(vals, SPEARmodel$factors$factor.scores$in.sample[,f], method = method)$p.value
          res$pval.out.of.sample[r] <- cor.test(vals, SPEARmodel$factors$factor.scores$out.of.sample[,f], method = method)$p.value
        }
        colnames(res) <- c("Feature", "Probability", "Coefficient", "Omic", "Factor", "in.sample.cor", "out.of.sample.cor", "in.sample.pval", "out.of.sample.pval")
        rownames(res) <- NULL
        res <- dplyr::select(res, Feature, Omic, Factor, Coefficient, Probability, in.sample.cor, out.of.sample.cor, in.sample.pval, out.of.sample.pval)
        results.list[[length(results.list) + 1]] <- res
      }
    }
  }
  if(length(results.list) == 0){
    warning("*** Warning: No features found within cutoffs for chosen Factors/Omics. Returning data.frame with 0 rows. Try broadening the criteria.")
    return(data.frame(Feature = character(0),
                      Omic = character(0),
                      Factor = character(0),
                      Probability = numeric(0),
                      Coefficient = numeric(0)))
  } else {
    total.results <- do.call("rbind", results.list)
    return(total.results)
  }
}




#' Plot posterior seleciton probabilities and coefficients for features that pass the provided cutoffs.
#'@param SPEARmodel SPEAR model (returned from \code{get_SPEAR_model})
#'@param factors A vector of integers representing factors. If NULL, will return features for all factors. Defaults to NULL (all factors)
#'@param omics A vector of names of datasets (see names(SPEARmodel$data$xlist)). Defaults to NULL (all omics)
#'@param coefficient.cutoff Cutoff for coefficient magnitude. If .1 is provided, abs(coefficients) >= .1 are returned. Defaults to .01
#'@param probability.cutoff Cutoff for posterior selection probability. Ranges from 0 - 1. Defaults to .5
#'@param show.facet Separate by Factor and Omic? Defaults to TRUE
#'@return A plot of features that pass the provided cutoffs.
#' @examples
#' SPEAR.plot_feature_distribution(SPEARmodel)
#' 
#' # Plot specific features for one omic:
#' omic.name <- names(SPEARmodel$data$X)[1]
#' SPEAR.plot_feature_distribution(SPEARmodel, factors = 1, omics = omic.name)
#' 
#' # Plot features that pass cutoffs:
#' SPEAR.plot_feature_distribution(SPEARmodel, coefficient.cutoff = .2, probability.cutoff = .5)
#'@export
SPEAR.plot_feature_distribution <- function(SPEARmodel, factors = NULL, omics = NULL, coefficient.cutoff = .01, probability.cutoff = .5, show.facet = TRUE){
  feature.table <- SPEAR.get_feature_table(SPEARmodel, factors = factors, omics = omics, coefficient.cutoff = 0, probability.cutoff = probability.cutoff)
  if(nrow(feature.table) == 0){
    stop("*** ERROR: No features found for requested cutoffs. Try broadening the cutoffs.")
  }
  
  feature.table$color <- sapply(feature.table$Coefficient, function(coeff){
    if(coeff >= coefficient.cutoff){
      return("high")
    } else if(coeff <= -coefficient.cutoff){
      return("low")
    } else {
      return("insignificant")
    }
  })
  
  # Plot results:
  g <- ggplot2::ggplot(feature.table) +
    ggplot2::annotate("rect", xmin = coefficient.cutoff, xmax = 100, ymin = -100, ymax = 100, fill = "#2ECC71", alpha = .2) +
    ggplot2::annotate("rect", xmin = -100, xmax = -coefficient.cutoff, ymin = -100, ymax = 100, fill = "#E94D3D", alpha = .2) +
    ggplot2::geom_point(ggplot2::aes(x = Coefficient, y = Probability, fill = color), shape = 21, stroke = .2, size = 2) +
    ggplot2::scale_fill_manual(values = list(high = "#2ECC71", low = "#E94D3D", insignificant = "grey"), guide = FALSE) +
    ggplot2::coord_cartesian(xlim = c(-max(abs(feature.table$Coefficient)), max(abs(feature.table$Coefficient))), ylim = c(probability.cutoff,1)) +
    ggplot2::theme_bw()
  
  if(show.facet){
    g <- g + ggplot2::facet_grid(rows = ggplot2::vars(Omic), cols = ggplot2::vars(Factor))
  }
  
  return(g)
}



#' Plot posterior seleciton probabilities and coefficients for features that pass the provided cutoffs.
#'@param SPEARmodel SPEAR model (returned from \code{get_SPEAR_model})
#'@param factors A vector of integers representing factors. If NULL, will return features for all factors. Defaults to NULL (all factors)
#'@param omics A vector of names of datasets (see names(SPEARmodel$data$xlist)). Defaults to NULL (all omics)
#'@param coefficient.cutoff Cutoff for coefficient magnitude. If .1 is provided, abs(coefficients) >= .1 are returned. Defaults to .01
#'@param probability.cutoff Cutoff for posterior selection probability. Ranges from 0 - 1. Defaults to .5
#'@param sort.by How to order the features? "probability" ranks by posterior selection probability first, then coefficient. "coefficient" ranks purely by coefficient. Defaults to "probability"
#'@param max.per.factor Maximum number of features per factor? Defaults to 10 (may be more due to overlapping features)
#'@return A plot of features that pass the provided cutoffs.
#' @examples
#' SPEAR.plot_feature_grid(SPEARmodel)
#' 
#' # Plot specific features for one omic:
#' omic.name <- names(SPEARmodel$data$X)[1]
#' SPEAR.plot_feature_grid(SPEARmodel, factors = 1, omics = omic.name)
#' 
#' # Plot features that pass cutoffs:
#' SPEAR.plot_feature_grid(SPEARmodel, coefficient.cutoff = .2, probability.cutoff = .5)
#'@export
SPEAR.plot_feature_grid <- function(SPEARmodel, factors = NULL, omics = NULL, coefficient.cutoff = .01, probability.cutoff = .5, sort.by = "probability", max.per.factor = 10){
  feature.table <- SPEAR.get_feature_table(SPEARmodel, factors = factors, omics = omics, coefficient.cutoff = .01, probability.cutoff = probability.cutoff)
  if(nrow(feature.table) == 0){
    stop("*** ERROR: No features found for requested cutoffs. Try broadening the cutoffs.")
  }
  
  if(sort.by == "probability"){
    feature.table <- dplyr::arrange(feature.table, -Probability, -abs(Coefficient))
  } else if(sort.by == "coefficient"){
    feature.table <- dplyr::arrange(feature.table, -abs(Coefficient))
  }
  
  # cut down to top max.per.factor:
  final.res <- data.frame(Feature = character(0),
                          Omic = character(0),
                          Factor = character(0),
                          Probability = numeric(0),
                          Coefficient = numeric(0))
  for(f in unique(feature.table$Factor)){
    temp <- dplyr::filter(feature.table, Factor == f)
    if(nrow(temp) > 0){
      if(nrow(temp) > max.per.factor){
        if(sort.by == "probability"){
          temp <- dplyr::arrange(temp, -Probability, -abs(Coefficient))
        } else if(sort.by == "coefficient"){
          temp <- dplyr::arrange(temp, -abs(Coefficient))
        }
        temp <- temp[1:max.per.factor,]
        final.res <- rbind(final.res, temp)
      } else {
        final.res <- rbind(final.res, temp)
      }
    }
  }
  feature.table <- final.res
  
  
  # Plot results:
  feature.table$Feature <- factor(feature.table$Feature, levels = rev(unique(feature.table$Feature)))
  g <- ggplot2::ggplot(feature.table) +
    ggplot2::geom_bar(ggplot2::aes(x = Coefficient, y = Feature, fill = Omic), stat = "identity") +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::scale_fill_manual(values = SPEAR.get_color_scheme(SPEARmodel)) +
    ggplot2::theme_bw() +
    ggplot2::facet_grid(rows = ggplot2::vars(Omic), cols = ggplot2::vars(Factor), scales="free_y")
    
  
  return(g)
}



#' Plot posterior seleciton probabilities and coefficients for features that pass the provided cutoffs.
#'@param SPEARmodel SPEAR model (returned from \code{get_SPEAR_model})
#'@param factors A vector of integers representing factors. If NULL, will return features for all factors. Defaults to NULL (all factors)
#'@param omics A vector of names of datasets (see names(SPEARmodel$data$xlist)). Defaults to NULL (all omics)
#'@param coefficient.cutoff Cutoff for coefficient magnitude. If .1 is provided, abs(coefficients) >= .1 are returned. Defaults to .01
#'@param probability.cutoff Cutoff for posterior selection probability. Ranges from 0 - 1. Defaults to .5
#'@param sort.by How to order the features? "probability" ranks by posterior selection probability first, then coefficient. "coefficient" ranks purely by coefficient. Defaults to "probability"
#'@param max.per.factor Maximum number of features per factor? Defaults to 10 (may be more due to overlapping features)
#'@return A plot of features that pass the provided cutoffs.
#' @examples
#' SPEAR.plot_feature_summary(SPEARmodel)
#' 
#' # Plot specific features for one omic:
#' omic.name <- names(SPEARmodel$data$X)[1]
#' SPEAR.plot_feature_summary(SPEARmodel, factors = 1, omics = omic.name)
#' 
#' # Plot features that pass cutoffs:
#' SPEAR.plot_feature_summary(SPEARmodel, coefficient.cutoff = .2, probability.cutoff = .5)
#'@export
SPEAR.plot_feature_summary <- function(SPEARmodel, factors = NULL, omics = NULL, coefficient.cutoff = .01, probability.cutoff = .5, sort.by = "probability", max.per.factor = 10){
  feature.table <- SPEAR.get_feature_table(SPEARmodel, factors = factors, omics = omics, coefficient.cutoff = .01, probability.cutoff = probability.cutoff)
  if(nrow(feature.table) == 0){
    stop("*** ERROR: No features found for requested cutoffs. Try broadening the cutoffs.")
  }

  if(sort.by == "probability"){
    feature.table <- dplyr::arrange(feature.table, -Probability, -abs(Coefficient))
  } else if(sort.by == "coefficient"){
    feature.table <- dplyr::arrange(feature.table, -abs(Coefficient))
  }
  
  # cut down to top max.per.factor:
  final.res <- list()
  for(f in sort(unique(feature.table$Factor))){
    temp <- dplyr::filter(feature.table, Factor == f)
    if(nrow(temp) > 0){
      if(nrow(temp) > max.per.factor){
        temp <- temp[1:max.per.factor,]
        final.res[[f]] <- temp
      } else {
        final.res[[f]] <- temp
      }
    }
  }

  plotlist <- list()
  for(i in 1:length(final.res)){
    final.res[[i]]$Feature <- factor(final.res[[i]]$Feature, levels = rev(final.res[[i]]$Feature))
    g <- ggplot2::ggplot(final.res[[i]]) +
      ggplot2::geom_bar(ggplot2::aes(y = Feature, x = Coefficient, fill = Omic), stat = "identity") +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::ggtitle(names(final.res)[i]) +
      ggplot2::scale_fill_manual(values = SPEAR.get_color_scheme(SPEARmodel), guide = FALSE) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = .5, size = 12))
    plotlist[[i]] <- g
  }
  # add legend:
  g.legend <- ggplot2::ggplot(data.frame(Omic = names(SPEARmodel$data$xlist), x = 1, y = 1)) +
    ggplot2::geom_bar(ggplot2::aes(x = x, y = y, fill = Omic), stat = "identity") +
    ggplot2::scale_fill_manual(values = SPEAR.get_color_scheme(SPEARmodel))
  plotlist[[length(plotlist) + 1]] <- cowplot::get_legend(g.legend)
  
  p <- cowplot::plot_grid(plotlist = plotlist, nrow = 1)
  return(p)
}



#' Plot ordinal class predictions
#'@param SPEARmodel A SPEAR model (returned from get_SPEAR_model)
#'@param forecast Which probabilities to use? A string with "out.of.sample" or "in.sample". Defaults to "out.of.sample".
#'@param Xlist A list of matrices of testing data (must match the same dimensions as the data used to train the SPEARmodel (SPEARmodel$data$xlist)). Defaults to NULL (use training data).
#'@param show.true.distribution Plot the true distribution as well? Defaults to TRUE
#'@export
SPEAR.plot_class_predictions <- function(SPEARmodel, forecast = "out.of.sample", Xlist = NULL, show.true.distribution = TRUE){
  ## Gaussian check: ----------
  if(SPEARmodel$params$family == "gaussian"){
    stop("ERROR: SPEARobj must be of 'family' type 'binomial', 'ordinal', or 'categorical'. Current SPEARmodel is of type 'gaussian'.")
  }
  if(SPEARmodel$params$family == "categorical"){
    # Get class predictions:
    preds <- SPEARmodel$predictions[[forecast]]$overall.class.predictions
    levels <- ncol(preds)
    class.preds <- data.frame(ClassPrediction = apply(preds, 1, which.max) - 1)
    rownames(class.preds) <- rownames(preds)
    
    # Get true classes:
    true.classes <- SPEARmodel$data$Y
    true.class.long <- data.frame(Class = apply(true.classes, 1, which.max) - 1)
    rownames(true.class.long) <- rownames(true.classes)
    
    if(!all(rownames(true.class.long) == rownames(class.preds))){
      stop("*** ERROR: Not all rownames of predicted classes match with rownames of true classes...")
    }
    
    plot.type <- data.frame(PlotType = rep(paste0("SPEAR Predictions (w = ", round(SPEARmodel$params$w, 3), ")"), nrow(class.preds)))
    rownames(plot.type) <- rownames(class.preds)

    df <- cbind(class.preds, true.class.long, plot.type)
    
    if(show.true.distribution){
      df.true <- cbind(true.class.long, true.class.long, data.frame(PlotType = rep("True Class Distribution", nrow(true.class.long))))
      colnames(df.true) <- colnames(df)
      df <- rbind(df, df.true)
    }
    
    labels <- colnames(SPEARmodel$data$Y)
    
    # Change Class to the name of the column (for the legend)
    df$Class <- colnames(SPEARmodel$data$Y)[df$Class + 1]
  } else { # Binomial / Ordinal
    ### Xlist (test data) ----------
    if(is.null(Xlist)){
      class.predictions <- SPEARmodel$predictions[[forecast]]
    } else {
      class.predictions <- SPEAR::SPEAR.get_predictions(SPEARmodel, Xlist)
    } # ------
    # Get class predictions:
    class.preds <- as.data.frame(class.predictions$class.predictions)
    levels <- ncol(class.predictions$class.probabilities[[1]])
    colnames(class.preds) <- c("ClassPrediction")
    
    # Get true classes:
    true.class.long <- SPEARmodel$data$Y
    colnames(true.class.long) <- c("Class")
    
    if(!all(rownames(true.class.long) == rownames(class.preds))){
      stop("*** ERROR: Not all rownames of predicted classes match with rownames of true classes...")
    }
    
    plot.type <- data.frame(PlotType = rep(paste0("SPEAR Predictions (w = ", round(SPEARmodel$params$w, 3), ")"), nrow(class.preds)))
    rownames(plot.type) <- rownames(class.preds)
    
    df <- cbind(class.preds, true.class.long, plot.type)
    
    if(show.true.distribution){
      df.true <- cbind(true.class.long, true.class.long, data.frame(PlotType = rep("True Class Distribution", nrow(true.class.long))))
      colnames(df.true) <- colnames(df)
      df <- rbind(df, df.true)
    }
    
    labels <- c(0:(levels-1))
    
    # Change Class to the name of the column (for the legend)
    df$Class <- factor(df$Class, levels = 0:(levels-1))
  }
  df$ClassPrediction <- factor(df$ClassPrediction, levels = 0:(levels-1))
  g <- ggplot2::ggplot(df) +
    ggplot2::geom_histogram(ggplot2::aes(x = ClassPrediction, fill = Class), color = "black", stat = "count", lwd = .25) +
    ggplot2::ylab("Count") +
    ggplot2::geom_segment(ggplot2::aes(x = 0.5, y = 0, xend = (levels+.5), yend = 0), lwd = 0) +
    ggplot2::scale_fill_brewer(palette = "RdBu", direction = 1) +
    ggplot2::scale_x_discrete(labels = labels, c(0:(levels-1))) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank()) +
    ggplot2::facet_wrap(ggplot2::vars(PlotType))
  
  return(g)
}



#' Check that data are normally distributed and scaled properly for SPEAR
#'@param Xlist A list of matrices of omics data. Defaults to NULL.
#'@param Y A data.frame of response variables. Defaults to NULL.
#'@param plot.xlist.pvalues Should the function return plots of p.values from the shapiro-wilk normality tests? Defaults to TRUE
#'@return NULL (will print out results)
#' @examples
#' SPEAR.get_data_summary(Xlist, Y)
#' SPEAR.get_data_summary(SPEARmodel$data$xlist, SPEARmodel$data$Y)
#' SPEAR.get_data_summary(SPEARobject$data$xlist, SPEARobject$data$Y)
#'@export
SPEAR.get_data_summary <- function(Xlist = NULL, Y = NULL, plot.xlist.pvalues = TRUE){
  cat("******************************************\n")
  cat("| SPEAR - Checking Data Distributions... |\n")
  cat("******************************************\n\n")
  
  # remove outliers per feature before doing shapiro wilk test
  
  # what percentage of data are outliers (3-4 SDs from mean)
  
  if(is.null(Xlist) & is.null(Y)){
    cat("*** Neither Xlist nor Y provided. Nothing to check... exiting...\n")
  } else {
    all.pvalues.x <- list()
    if(!is.null(Xlist)){
      cat(paste0("*** Xlist provided: Detected ", length(Xlist), " datasets:\n\n"))
      for(i in 1:length(Xlist)){
        cat("- - - - - - - - - - - - - - - - - - - - - - - - - -\n")
        cat(paste0(names(Xlist)[i], "\tSubjects: ", nrow(Xlist[[i]]), "\tFeatures: ", ncol(Xlist[[i]]), "\n"))
        # Check if mean is centered at 0 for all features (within 10 decimals)
        cat("~~~ Mean-centered at 0?\t", all(round(apply(Xlist[[i]], 2, mean), 10) == 0), "\n")
        cat("~~~ Unit Variance?\t", all(round(apply(Xlist[[i]], 2, var), 10) == 1), "\n")
        cat("~~~ Shapiro-Wilk Normality Test:\n")
        shapiro.tests <- apply(Xlist[[i]], 2, stats::shapiro.test)
        pvals <- sapply(shapiro.tests, function(test){return(test$p.value)})
        summary.pvals <- summary(pvals)
        cat("Stat\t p.value      Sig.? (<.05)\n----------------------------------\n")
        cat("Mean:\t", signif(summary.pvals[["Mean"]], 5), ifelse(summary.pvals[["Mean"]] < .05, "\t*", "\t"), "\n")
        cat("Min:\t", signif(summary.pvals[["Min."]], 5), ifelse(summary.pvals[["Min."]] < .05, "\t*", "\t"), "\n")
        cat("1stQ:\t", signif(summary.pvals[["1st Qu."]], 5), ifelse(summary.pvals[["1st Qu."]] < .05, "\t*", "\t"), "\n")
        cat("Median:\t", signif(summary.pvals[["Median"]], 5), ifelse(summary.pvals[["Median"]] < .05, "\t*", "\t"), "\n")
        cat("3rdQ:\t", signif(summary.pvals[["3rd Qu."]], 5), ifelse(summary.pvals[["3rd Qu."]] < .05, "\t*", "\t"), "\n")
        cat("Max:\t", signif(summary.pvals[["Max."]], 5), ifelse(summary.pvals[["Max."]] < .05, "\t*", "\t"), "\n")
        cat("\n")
        # Save p.values
        all.pvalues.x[[names(Xlist)[i]]] <- pvals
      }
    }
    if(!is.null(Y)){
      cat(paste0("*** Y provided: Detected ", ncol(Y), " variables:\n\n"))
      for(i in 1:ncol(Y)){
        if(all(Y[,i] %in% 0:1)){
          if(ncol(Y) > 1){
            family <- "categorical"
          } else {
            family <- "binomial"
          }
        } else if(all(Y[,i] %in% 0:max(Y[,i]))) {
          family <- "ordinal"
        } else {
          family <- "gaussian" 
        }
        cat("- - - - - - - - - - - - - - - - - - - - - - - - - -\n")
        cat(paste0(colnames(Y)[i], "\tSubjects: ", nrow(Y), "\tFamily: ", family, "\n"))
        if(family == "gaussian"){
          # Check if mean is centered at 0 for all features (within 10 decimals)
          cat("~~~ Mean-centered at 0?\t", round(mean(Y[,i]), 10) == 0, "\n")
          cat("~~~ Unit Variance?\t", round(mean(Y[,i]), 10) == 1, "\n")
          cat("~~~ Shapiro-Wilk Normality Test:\n")
          shapiro.test <- shapiro.test(Y[,i])
          pval <- shapiro.test$p.value
          cat("Stat\t p.value      Sig.? (<.05)\n----------------------------------\n")
          cat("Pval:\t", signif(pval, 5), ifelse(pval < .05, "\t*", "\t"), "\n")
          cat("\n")
        } else {
          t <- table(Y[,i])
          t.all <- sapply(0:max(Y[,i]), function(class.temp){
            if(!class.temp %in% names(t)){
              return(0)
            } else {
              return(t[which(names(t) == class.temp)])
            }
          })
          for(k in 1:length(t.all)){
            cat("Class ", k-1, ": ", t.all[k], "\n", sep = "") 
          }
        }
      }
    }
    
    if(plot.xlist.pvalues){
      df.list <- list()
      for(o in 1:length(all.pvalues.x)){
        df.temp <- data.frame(pval = all.pvalues.x[[o]], omic = paste0("Xlist[[", o, "]]: ", names(Xlist)[o]))
        df.list[[o]] <- df.temp
      }
      df <- do.call("rbind", df.list)

      g <- ggplot2::ggplot(data = df) +
        ggplot2::geom_histogram(ggplot2::aes(x = pval), fill = "lightgrey", color = "black", binwidth = .1) +
        ggplot2::theme_bw() +
        ggplot2::facet_wrap(ggplot2::vars(omic), scales = "free_y")
      
      return(g)
    }
  }
}



#' Plot a correlation matrix for the selected features
#'@param SPEARmodel A SPEAR model (returned from get_SPEAR_model)
#'@param feature.table A table of features from SPEAR.get_feature_table
#'@param show.feature.names Should the feature names be shown on the axes? Defaults to TRUE
#'@param show.correlation.values Should the (rounded) correlation values be shown? Defaults to FALSE
#'@param method Method for correlation. Can be 'pearson', 'kendall', or 'spearman'. Defaults to 'spearman'.
#'@param cluster.features Cluster the features in the heatmap? Defaults to TRUE
#'@export
SPEAR.plot_feature_correlation <- function(SPEARmodel, feature.table, method = "spearman", show.feature.names = TRUE, show.correlation.values = FALSE, cluster.features = TRUE){

    if(length(unique(feature.table$Omic)) > 1){
    cat("*** WARNING: More than one Omic present in 'feature.table'. Combining all features...")
  }
  
  # Get matrix of samples x features:
  mat <- matrix(NA, nrow = nrow(SPEARmodel$data$X), ncol = nrow(feature.table))
  colnames(mat) <- feature.table$Feature
  rownames(mat) <- rownames(SPEARmodel$data$X)
  for(j in 1:ncol(mat)){
    mat[,j] <- SPEARmodel$data$X[,which(colnames(SPEARmodel$data$X) == colnames(mat)[j])]
  }
  # Make correlation matrix from features provided:
  corr.mat <- matrix(NA, nrow = ncol(mat), ncol = ncol(mat))
  colnames(corr.mat) <- colnames(mat)
  rownames(corr.mat) <- colnames(mat)
  for(i in 1:ncol(mat)){
    for(j in 1:ncol(mat)){
      corr.mat[i,j] <- cor(mat[,i], mat[,j], method = method)
    }
  }
  
  if(cluster.features){
    hc = hclust(as.dist(1 - corr.mat))
    corr.mat = corr.mat[hc$order, hc$order]
  }
  # remove lower triangle
  corr.mat[lower.tri(corr.mat)] = NA
  
  mat.melt <- reshape2::melt(corr.mat)
  g <- ggplot2::ggplot(mat.melt) +
    ggplot2::geom_tile(ggplot2::aes(x = Var1, y = Var2, fill = value)) +
    ggplot2::scale_fill_gradient2(high = "dark red", low = "dark blue", limits = c(-1, 1), na.value = "white", name = method) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab(NULL) +
    ggplot2::coord_equal() +
    ggplot2::theme_void() +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 10),
                   legend.position = c(.75, .25))
    
  if(show.feature.names){
    g <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, size = 8, hjust = 0),
                   axis.text.y = ggplot2::element_text(size = 8, hjust = 1))
  }
  
  if(show.correlation.values){
    g <- g + ggplot2::geom_text(ggplot2::aes(x = Var1, y = Var2, label = round(value, 2)), size = 2.5)
  }

  
  return(g)
}

#' Get the best features from the variance cutoff parameter
#'@param SPEARmodel A SPEAR model (returned from get_SPEAR_model)
#'@param coefficient.cutoff Cutoff for coefficient magnitude. If .1 is provided, abs(coefficients) >= .1 are returned. Defaults to .01
#'@param probability.cutoff Cutoff for posterior selection probability. Ranges from 0 - 1. Defaults to .5
#'@export
SPEAR.get_best_features <- function(SPEARmodel, factor = 1, var.cutoff = .5, forecast = "out.of.sample"){
  # First step: get all features:
  feature.table = SPEAR.get_feature_table(SPEARmodel, factors = factor, coefficient.cutoff = 0, probability.cutoff = 0)
  
  # Sort by Coefficient (Beta)
  feature.table = dplyr::arrange(feature.table, -Probability, -abs(Coefficient))
  
  # Now, iterate step by step...
  var.exp = 0
  features = c()
  current.residuals = SPEARmodel$factors$factor.scores[[forecast]][,paste0('Factor', factor)]
  orig.var = var(SPEARmodel$factors$factor.scores[[forecast]][,paste0('Factor', factor)])
  while(var.exp < var.cutoff){
    # Add a feature:
    new.feature = feature.table$Feature[length(features) + 1]
    features = c(features, new.feature)
    # Regress factor by feature:
    current.residuals <- lm(current.residuals ~ SPEARmodel$data$X[, new.feature])$residuals
    # Get variance explained by all features up until now:
    var.exp = 1 - (var(current.residuals)/orig.var)
    # print out:
    cat("~~~ ", length(features), ") Adding ", SPEAR.color_text(new.feature, "light cyan"), "...\t Total Var. Exp. Factor", factor, " = ", var.exp, "\n", sep = "")
  }
  
  # Return the feature table with the features chosen:
  feature.table <- dplyr::filter(feature.table, Feature %in% features)
  
  return(feature.table)
}

#' Plot variance per feature
#'@param SPEARmodel A SPEAR model (returned from get_SPEAR_model)
#'@param coefficient.cutoff Cutoff for coefficient magnitude. If .1 is provided, abs(coefficients) >= .1 are returned. Defaults to .01
#'@param probability.cutoff Cutoff for posterior selection probability. Ranges from 0 - 1. Defaults to .5
#'@export
SPEAR.plot_variance_per_feature <- function(SPEARmodel, factor = 1, num.features = 50, coefficient.cutoff = 0, probability.cutoff = 0, forecast = "out.of.sample"){
  # First step: get all features:
  feature.table = SPEAR::SPEAR.get_feature_table(SPEARmodel, factors = factor, coefficient.cutoff = coefficient.cutoff, probability.cutoff = probability.cutoff)
  
  # Sort by Coefficient (Beta)
  feature.table = dplyr::arrange(feature.table, Probability, -abs(Coefficient))
  
  # Now, iterate step by step...
  var.exp = 0
  features = c()
  current.residuals = SPEARmodel$factors$factor.scores[[forecast]][,paste0('Factor', factor)]
  orig.var = var(SPEARmodel$factors$factor.scores[[forecast]][,paste0('Factor', factor)])
  var.exp.list <- c()
  var.added.list <- c()
  for(i in 1:num.features){
    # Add a feature:
    new.feature = feature.table$Feature[i]
    # Regress factor by feature:
    current.residuals <- lm(current.residuals ~ SPEARmodel$data$X[, new.feature])$residuals
    # Get variance explained by all features up until now:
    var.exp = 1 - (var(current.residuals)/orig.var)
    var.exp.list <- c(var.exp.list, var.exp)
    if(i != 1)
      var.added.list <- c(var.added.list, var.exp - var.exp.list[i - 1])
    else
      var.added.list <- c(var.added.list, var.exp)
    # print out:
    cat("*** Trying with ", new.feature, "...\t Var Explained for Factor", factor, " = ", var.exp, "\n")
  }
  df <- feature.table[1:length(var.exp.list),]
  df$var.exp <- var.exp.list
  df$var.added <- var.added.list

  g <- ggplot2::ggplot(df) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 1:nrow(df) - .5,
                                    xmax = 1:nrow(df) + .5,
                                    ymin = 0,
                                    ymax = var.exp, 
                                    fill = Omic), color = "black", size = .1) +
    ggplot2::scale_x_continuous(breaks = 1:nrow(df), labels = df$Feature) +
    ggplot2::scale_color_manual(values = SPEAR.get_color_scheme(SPEARobj)) +
    ggplot2::scale_fill_manual(values = SPEAR.get_color_scheme(SPEARobj)) +
    ggplot2::theme_bw() +
    #ggplot2::ylim(c(0, 1)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
  
  return(g)
}

#SPEARobj <- readRDS("~/Dropbox/shared_folder/RealData/COVID-19/SPEAR_COVID_RF_Imputed_Ordinal.rds")
#SPEARmodel <- SPEAR::get_SPEAR_model(SPEARobj, w = 2)

#SPEAR.plot_variance_per_feature(SPEARmodel, factor = 1, num.features = 50)
#SPEAR.get_best_features(SPEARmodel)
#SPEARobj = readRDS("~/Documents/Coding/IMPACC/SPEAR/TCGA_tests/041421/SPEAR_tcga_encode_cat_0417_JG.rds")
#SPEARmodel = SPEAR::get_SPEAR_model(SPEARobj, w = 2)
#feature.table <- SPEAR.get_feature_table(SPEARmodel, factors = 2, omics = "mrna", probability.cutoff = .5, coefficient.cutoff = .0001, method = "pearson")
#feature.table <- dplyr::arrange(feature.table, -abs(Coefficient))
#feature.table <- dplyr::arrange(feature.table, -Probability, -abs(Coefficient))
#feature.table <- feature.table[1:20,]
#SPEAR.plot_feature_correlation(SPEARmodel, 
#                               feature.table, 
#                               show.feature.names = TRUE, 
#                               show.correlation.values = TRUE, 
#                               cluster.features = TRUE)



