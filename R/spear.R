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


#### Functions:
#' Run SPEAR
#' @examples
#' SPEARobj <- make_spear_model(...)
#' 
#' SPEARobj$run.spear()
#' 
#'@export
run.spear <- function(){
  X = self$data$X
  Y = self$data$Y
  Z = self$data$Z
  Xobs = self$data$Xobs
  Yobs = self$data$Yobs
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
  post_mu = array(0, dim = c(ncol(Z), num_factors))
  post_sigma2 = array(0.1, dim=c(ncol(Z), num_factors));
  post_pi = array(1, dim=c(ncol(Z), num_factors));
  if(!is.null(inits_post_mu)){
    if((ncol(inits_post_mu)!=num_factors) | (nrow(inits_post_mu)!= ncol(Z))){
      stop("wrong initialization dimension for post_mu!")
    }
    post_mu = inits_post_mu;
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
  ##record the model estimated with weights all 1 for initial start of y
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
    set.seed(seed)
    
    if(print_out > 0)
      cat(paste0("\n--- ", self$color.text(paste0("Running weight.x = ", all_ws[idx_w,1], " | weight.y = ",all_ws[idx_w,2]), "green"), " ---\n"))
    
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

  self$fit <- list(regression.coefs = post_betas, 
                   projection.coefs.x = post_bxs, 
                   projection.coefs.y = post_bys,
                   post_pis = post_pis, 
                   projection.probs = post_selections, 
                   marginal.probs = post_selections_marginal,
                   joint.probs = post_selections_joint,
                   intercepts.x = interceptsX, 
                   intercepts.y = interceptsY)
  
  self$update.dimnames()
  
  self$options$current.weight.idx = 1
}

#' Update the dimension names for all SPEARobject matrices. Used internally.
#'@export
update.dimnames = function(call = "fit"){
  if(call == "fit"){
    dimnames(self$fit$regression.coefs) = list()
    dimnames(self$fit$projection.coefs.x) = list() 
    dimnames(self$fit$projection.coefs.y) = list() 
    dimnames(self$fit$post_pis) = list() 
    dimnames(self$fit$projection.probs) = list() 
    dimnames(self$fit$marginal.probs) = list() 
    dimnames(self$fit$joint.probs) = list()
  }
}

#+ post_betas ($P$ x $K$ x $G$ x $W$)
#
#+ post_bys ($K$ x $Y$ x $W$)
#
#+ post_bxs ($K$ x $P$ x $W$)
#
#+ post_pis ($P$ x $K$ x $W$)
#
#+ post_selections ($P$ x $K$ x $W$)


#### Functions:
#' Print out a variety of summary information about a SPEARobject
#' @param type Which type of summary to print? Can be "data". Defaults to NULL.
#' @param remove.formatting Remove text color/bold font face. Defaults to FALSE.
#' @param quiet Do not print anything. Defaults to FALSE.
#' @examples
#' SPEARobj <- make_spear_model(...)
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
    cat(paste0("Detected ", length(self$data$Xlist), " datasets:\n"))
    for(i in 1:length(self$data$Xlist)){
      cat(self$color.text(names(self$data$Xlist)[i], dataset.color), "\tSubjects: ", nrow(self$data$Xlist[[i]]), "\tFeatures: ", ncol(self$data$Xlist[[i]]), "\n")
    }
    cat(paste0("Detected ", ncol(self$data$Y), " response ", ifelse(ncol(self$data$Y) == 1, "variable", "variables"), ":\n"))
    for(i in 1:ncol(self$data$Y)){
      cat(self$color.text(colnames(self$data$Y)[i], response.color), "\tSubjects: ", sum(!is.na(self$data$Y[,i])), "\tType: ", self$params$family, "\n")
    }
  }
  # ---------------------------------
  if(type == "help"){
    cat("For assistance, browse the SPEAR vignettes\nType browseVignettes('SPEAR')")
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


# Define SPEARobject class:
SPEARobject <- R6::R6Class("SPEARobject",
            public = list(
              data = NULL,
              params = NULL,
              inits = NULL,
              options = NULL,
              
              # This needs to be another object...
              fit = NULL,
              # This needs to be another object...
              cv.fit = NULL, # has eval member
              
              # Functions:
              initialize = function(
                                    # data:
                                    X = NULL, 
                                    Y = NULL, 
                                    Z = NULL, 
                                    # data descriptors:
                                    family = "gaussian", 
                                    # weights:
                                    weights.case = NULL,  
                                    weights.x = NULL, 
                                    weights.y = NULL,
                                    # model parameters:
                                    num_factors = 5, 
                                    inits_type = "pca", 
                                    inits_post_mu = NULL,
                                    sparsity_upper = 0.5,
                                    warm_up = 100, 
                                    max_iter = 1000,
                                    thres_elbo = 0.01, 
                                    thres_count = 5, 
                                    thres_factor = 1e-8, 
                                    print_out = 100,
                                    seed = 123, 
                                    # coefficients:
                                    a0 = NULL, 
                                    b0 = NULL, 
                                    a1 = NULL, 
                                    b1 = NULL,
                                    a2 = NULL, 
                                    b2 = NULL, 
                                    L0 = NULL,
                                    L1 = NULL,
                                    L2 = NULL,
                                    robust_eps = NULL,
                                    remove.formatting = FALSE,
                                    quiet = FALSE
              ) {
                # called by SPEARobj$new(...)
                
                # Start with options:
                options = list()
                # remove.formatting - display text without color/bold? Defaults to FALSE
                options$remove.formatting = remove.formatting
                # quiet - should extra print statements be silenced? Defaults to FALSE
                options$quiet = quiet
                self$options = options
                
                if(!quiet){
                  cat("----------------------------------------------------------------\n")
                  cat("SPEAR version 1.1.0   Please direct all questions to Jeremy Gygi\n(", self$color.text("jeremy.gygi@yale.edu", ifelse(remove.formatting, "black", "green")), ") or Leying Guan (", self$color.text("leying.guan@yale.edu", ifelse(remove.formatting, "black", "green")), ")\n", sep = "")
                  cat("----------------------------------------------------------------\n")
                  cat("Generating SPEARobject...\n")
                }
                
                # Data:
                if(!quiet){cat("$data...\t")}
                data <- list()
                
                if(any(sapply(X, function(X.d){return(class(X.d)[1])}) != "matrix")){
                  X <- lapply(X, as.matrix)
                }
                data$Xlist = X
                data$X = do.call("cbind", X)
                data$Xobs <- array(1, dim  = dim(data$X))
                if(any(is.na(data$X))){
                  data$Xobs[which(is.na(data$X))] <- 0
                } else if(any(is.nan(data$X))){
                  data$Xobs[which(is.nan(data$X))] <- 0
                }
                
                if(is.null(dim(Y))){
                  data$Y = matrix(Y, ncol = 1)
                  colnames(data$Y) = "Y"
                  rownames(data$Y) = rownames(data$X)
                } else if(class(Y)[1] != "matrix"){
                  data$Y = as.matrix(Y)
                } else{
                  data$Y = Y
                }
                data$Yobs <- array(1, dim  = dim(data$Y))
                if(any(is.na(data$Y))){
                  data$Yobs[which(is.na(data$Y))] <- 0
                } else if(any(is.nan(data$Y))){
                  data$Yobs[which(is.nan(data$Y))] <- 0
                }
                
                # generate full matrix Z (impute if necessary...)
                if(is.null(Z)){
                  Z = data$X
                }
                # check for missing values in Z...
                if(any(is.na(Z)) | any(is.nan(Z))){
                  # if missing, IMPUTE Z!
                }
                # else...
                data$Z = Z
                
                if(!quiet){cat("Done!\n")}
                if(!quiet){cat("$params...\t")}
                # Parameters:
                params <- list()
                
                # family encoded
                if(is.numeric(family)){
                  if(!family %in% 0:3){
                    stop("ERROR: Family provided (", family, ') is not accepted. Must be:\n
                         0 | "gaussian"\n
                         1 | "binomial"\n
                         2 | "ordinal"\n
                         3 | "multinomial"\n')
                  } else if (family == 0){
                    params$family_encoded = 0
                    params$family = "gaussian"
                  } else if (family == 1){
                    params$family_encoded = 1
                    params$family = "binomial"
                  } else if (family == 2){
                    params$family_encoded = 2
                    params$family = "ordinal"
                  } else if (family == 3){
                    params$family_encoded = 3
                    params$family = "multinomial"
                    if(any(rowSums(data$Y) != 0)){
                      stop("ERROR: Values in Y do not follow a multinomial structure. All row sums must be equal to 1.")
                    }
                  }
                } else {
                  if(!family %in% c("gaussian", "binomial", "ordinal", "multinomial")){
                    stop("ERROR: Family provided (", family, ') is not accepted. Must be:\n
                         0 | "gaussian"\n
                         1 | "binomial"\n
                         2 | "ordinal"\n
                         3 | "multinomial"\n')
                  } else if (family == "gaussian"){
                    params$family_encoded = 0
                    params$family = "gaussian"
                  } else if (family == "binomial"){
                    params$family_encoded = 1
                    params$family = "binomial"
                  } else if (family == "ordinal"){
                    params$family_encoded = 2
                    params$family = "ordinal"
                  } else if (family == "multinomial"){
                    params$family_encoded = 3
                    params$family = "multinomial"
                    if(any(rowSums(data$Y) != 0)){
                      stop("ERROR: Values in Y do not follow a multinomial structure. All row sums must be equal to 1.")
                    }
                  }
                }
                
                # How to determine the number of factors by default?
                params$num_factors = num_factors
                
                params$nclasses = sapply(1:ncol(data$Y), function(j){
                  if(params$family == "gaussian"){
                    return(2)
                  } else if(params$family == "binomial"){
                    if(any(!unique(data$Y[,j] %in% 0:1))){
                      stop("ERROR: Values in Y do not follow an binomial structure. Values must be either 0 or 1.")
                    }
                    return(2)
                  } else if(params$family == "ordinal"){
                    labs <- unique(data$Y[,j])
                    if(any(!labs %in% 0:(length(labs-1)))){
                      stop("ERROR: Values in Y do not follow an ordinal structure. Values must be 0 - (num.classes-1) and not skip any classes.")
                    }
                    return(length(labs))
                  } else if(params$family == "multinomial"){
                    if(any(!unique(data$Y[,j] %in% 0:1))){
                      stop("ERROR: Values in Y do not follow a multinomial structure. Values must be either 0 or 1")
                    }
                    return(2)
                  } 
                })
                
                
                params$functional_path = list()
                start.ind = 1
                for(d in 1:length(data$Xlist)){
                  end.ind = start.ind + ncol(data$Xlist[[d]]) - 1
                  params$functional_path[[d]] <- start.ind:end.ind
                  start.ind = end.ind + 1
                }
                
                # Weights:
                if(is.null(weights.x) & is.null(weights.y)){
                  weights.x <- c(0, .01, .1, .5, 1, 2)
                  weights.y <- rep(1, length(weights.x))
                } else if(is.null(weights.x)){
                  weights.x <- rep(1, length(weights.y))
                } else if(is.null(weights.y)){
                  weights.y = rep(1, length(weights.x))
                } else if(length(weights.x) != length(weights.y)){
                  stop("ERROR: lengths of weights.x and weights.y do not match. They need to have the same length.")
                }
                params$weights = cbind(weights.x, weights.y)
                params$weights = params$weights[order(params$weights[,1], decreasing = TRUE),]
                colnames(params$weights) = c("w.x", "w.y")
                
                if(is.null(weights.case)){
                  params$weights.case = rep(1, nrow(data$X))
                }else if(length(weights.case)!=nrow(data$X)){
                  stop("ERROR: Supplied weights.case do not match the data dimension.")
                }else{
                  params$weights.case = weights.case
                }
                
                # Seed
                params$seed = seed
                
                # Misc. Parameters:
                params$inits_type = inits_type
                params$inits_post_mu = inits_post_mu
                params$sparsity_upper = sparsity_upper
                params$warm_up = warm_up
                params$max_iter = max_iter
                params$thres_elbo = thres_elbo
                params$thres_count = thres_count
                params$thres_factor =thres_factor
                params$print_out = print_out
                
                
                # Initial Coefficients:
                if(!quiet){cat("Done!\n")}
                if(!quiet){cat("$inits...\t")}
                inits = list()
                if(is.null(a0)){inits$a0 = 1e-2}else{inits$a0 = a0}
                if(is.null(b0)){inits$b0 = 1e-2}else{inits$b0 = b0}
                if(is.null(a1)){inits$a1 = sqrt(nrow(data$X))}else{inits$a1 = a1}
                if(is.null(b1)){inits$b1 = sqrt(nrow(data$X))}else{inits$b1 = b1}
                if(is.null(a2)){inits$a2 = sqrt(nrow(data$X))}else{inits$a2 = a2}
                if(is.null(b2)){inits$b2 = sqrt(nrow(data$X))}else{inits$b2 = b2}
                if(is.null(L0)){inits$L0 = 1}else{inits$L0 = L0}
                if(is.null(L1)){inits$L1 = nrow(data$X)/inits$L0}else{inits$L1 = L1}
                if(is.null(L2)){inits$L2 = nrow(data$X)/log(ncol(data$X))}else{inits$L2 = L2}
                if(is.null(robust_eps)){inits$robust_eps = 1.0/(nrow(data$X))}else{inits$robust_eps = robust_eps}
                
                # Save:
                self$data = data
                self$params = params
                self$inits = inits
                
                if(!quiet){
                  cat("Done!\n")
                  cat("SPEARobject generated!\n\n")
                }
                self$print.out(type = "data", remove.formatting = self$options$remove.formatting, quiet = self$options$quiet)
                cat("\n")
                self$print.out(type = "help", remove.formatting = self$options$remove.formatting, quiet = self$options$quiet)
              },
              
              # Method functions:
              print.out = print.out,
              color.text = color.text,
              update.dimnames = update.dimnames,
              run.spear = run.spear
              
            ) # end public
)

#' Make a SPEARobject. Will return an R6 class SPEARobject used for the "SPEAR" package.
#'@param X Assay matrix.
#'@param Y Response matrix (can be multidimensional for gaussian data).
#'@param Z Complete feature matrix (usually the features are the imputed version of X, other features are attached to the end).
#'@export
make_spear_object <- function(
  # data:
  X = NULL, 
  Y = NULL, 
  Z = NULL, 
  # data descriptors:
  family = "gaussian", 
  # weights:
  weights.case = NULL,  
  weights.x = NULL, 
  weights.y = NULL,
  # model parameters:
  num_factors = 5, 
  inits_type = "pca", 
  inits_post_mu = NULL,
  sparsity_upper = 0.5,
  warm_up = 100, 
  max_iter = 1000,
  thres_elbo = 0.01, 
  thres_count = 5, 
  thres_factor = 1e-8, 
  print_out = 100,
  seed = 123, 
  # coefficients:
  a0 = NULL, 
  b0 = NULL, 
  a1 = NULL, 
  b1 = NULL,
  a2 = NULL, 
  b2 = NULL, 
  L0 = NULL,
  L1 = NULL,
  L2 = NULL,
  robust_eps = NULL,
  # options:
  remove.formatting = FALSE,
  quiet = FALSE
){
  return(SPEARobject$new(
    # data:
    X = X, 
    Y = Y, 
    Z = Z, 
    # data descriptors:
    family = family, 
    # weights:
    weights.case = weights.case,  
    weights.x = weights.x, 
    weights.y = weights.y,
    # model parameters:
    num_factors = num_factors, 
    inits_type = inits_type, 
    inits_post_mu = inits_post_mu,
    sparsity_upper = sparsity_upper,
    warm_up = warm_up, 
    max_iter = max_iter,
    thres_elbo = thres_elbo, 
    thres_count = thres_count, 
    thres_factor = thres_factor, 
    print_out = print_out,
    seed = seed, 
    # coefficients:
    a0 = a0, 
    b0 = b0, 
    a1 = a1, 
    b1 = b1,
    a2 = a2, 
    b2 = b2, 
    L0 = L0,
    L1 = L1,
    L2 = L2,
    robust_eps = robust_eps,
    # options:
    remove.formatting = remove.formatting,
    quiet = quiet
  ))
}