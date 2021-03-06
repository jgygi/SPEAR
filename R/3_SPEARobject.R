# Define SPEARobject class:
SPEARobject <- R6::R6Class("SPEARobject",
                           public = list(
                             data = NULL,
                             params = NULL,
                             inits = NULL,
                             options = NULL,
                             
                             # Update with run.spear or run.cv.spear
                             fit = NULL,
                             # Update with get.model
                             model = NULL,
                             
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
                                 cat("SPEAR version 1.1.0   Please direct all questions to Jeremy Gygi\n(", private$color.text("jeremy.gygi@yale.edu", ifelse(remove.formatting, "black", "green")), ") or Leying Guan (", private$color.text("leying.guan@yale.edu", ifelse(remove.formatting, "black", "green")), ")\n", sep = "")
                                 cat("----------------------------------------------------------------\n")
                                 cat("Generating SPEARobject...\n")
                               }
                               
                               # Data:
                               if(!quiet){cat("$data...\t")}
                               self$data = list()
                               
                               self$add.data(X = X, Y = Y, Z = Z, name = "train")
                               
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
                                   if(any(rowSums(self$data$train$Y) != 1)){
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
                                   if(any(rowSums(self$data$train$Y) != 1)){
                                     stop("ERROR: Values in Y do not follow a multinomial structure. All row sums must be equal to 1.")
                                   }
                                 }
                               }
                               
                               # How to determine the number of factors by default?
                               params$num_factors = num_factors
                               
                               params$nclasses = sapply(1:ncol(self$data$train$Y), function(j){
                                 if(params$family == "gaussian"){
                                   return(2)
                                 } else if(params$family == "binomial"){
                                   if(any(!unique(self$data$train$Y[,j] %in% 0:1))){
                                     stop("ERROR: Values in Y do not follow an binomial structure. Values must be either 0 or 1.")
                                   }
                                   return(2)
                                 } else if(params$family == "ordinal"){
                                   labs <- unique(self$data$train$Y[,j])
                                   if(any(!labs %in% 0:(length(labs-1)))){
                                     stop("ERROR: Values in Y do not follow an ordinal structure. Values must be 0 - (num.classes-1) and not skip any classes.")
                                   }
                                   return(length(labs))
                                 } else if(params$family == "multinomial"){
                                   if(any(!unique(self$data$train$Y[,j] %in% 0:1))){
                                     stop("ERROR: Values in Y do not follow a multinomial structure. Values must be either 0 or 1")
                                   }
                                   return(2)
                                 } 
                               })
                               
                               
                               params$functional_path = list()
                               start.ind = 1
                               for(d in 1:length(self$data$train$Xlist)){
                                 end.ind = start.ind + ncol(self$data$train$Xlist[[d]]) - 1
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
                               rownames(params$weights) = 1:nrow(params$weights)
                               
                               if(is.null(weights.case)){
                                 params$weights.case = rep(1, nrow(self$data$train$X))
                               }else if(length(weights.case)!=nrow(self$data$train$X)){
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
                               if(is.null(a1)){inits$a1 = sqrt(nrow(self$data$train$X))}else{inits$a1 = a1}
                               if(is.null(b1)){inits$b1 = sqrt(nrow(self$data$train$X))}else{inits$b1 = b1}
                               if(is.null(a2)){inits$a2 = sqrt(nrow(self$data$train$X))}else{inits$a2 = a2}
                               if(is.null(b2)){inits$b2 = sqrt(nrow(self$data$train$X))}else{inits$b2 = b2}
                               if(is.null(L0)){inits$L0 = 1}else{inits$L0 = L0}
                               if(is.null(L1)){inits$L1 = nrow(self$data$train$X)/inits$L0}else{inits$L1 = L1}
                               if(is.null(L2)){inits$L2 = nrow(self$data$train$X)/log(ncol(self$data$train$X))}else{inits$L2 = L2}
                               if(is.null(robust_eps)){inits$robust_eps = 1.0/(nrow(self$data$train$X))}else{inits$robust_eps = robust_eps}
                               
                               # Save:
                               self$params = params
                               self$inits = inits
                               
                               if(!quiet){
                                 cat("Done!\n")
                                 cat("SPEARobject generated!\n\n")
                               }
                               private$print.out(type = "data", remove.formatting = self$options$remove.formatting, quiet = self$options$quiet)
                               cat("\n")
                               private$print.out(type = "help", remove.formatting = self$options$remove.formatting, quiet = self$options$quiet)
                             },
                             
                             # Method functions:
                             run.spear = run.spear,
                             run.cv.spear = run.cv.spear,
                             cv.evaluate = cv.evaluate,
                             set.weights = set.weights,
                             add.data = add.data,
                             save.model = save.model,
                             get.factor.scores = get.factor.scores,
                             get.predictions = get.predictions,
                             get.features = get.features,
                             get.contributions = get.contributions,
                             get.misclassification = get.misclassification
                             
                           ), # end public
                           private = list(
                             spear = spear,
                             cv.spear = cv.spear,
                             print.out = print.out,
                             color.text = color.text,
                             update.dimnames = update.dimnames,
                             impute.z = impute.z,
                             check.fold.ids = check.fold.ids,
                             check.fit = check.fit
                           )
)


#' Make a SPEARobject. Will return an R6 class SPEARobject used for the "SPEAR" package.
#'@param X Assay matrix.
#'@param Y Response matrix (can be multidimensional for gaussian data).
#'@param Z Complete feature matrix (usually the features are the imputed version of X, other features are attached to the end).
#'@export
make.SPEARobject <- function(
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


#' Load a SPEARobject. Will return an R6 class SPEARobject used for the "SPEAR" package.
#'@param file Where the SPEARobject Rds file is located. Defaults to NULL
#'@export
load.SPEARobject <- function(file = NULL){
  if(is.null(file)){
    stop("ERROR: File SPEARobject .rds file not found at ", file)
  }
  SPEARobj <- readRDS(file = file)
  if(is.null(SPEARobj$data)){
    stop("ERROR: SPEARobject loaded has no data. Are you sure you loaded a SPEARobject saved with saveRDS or $save.model()?")
  }
  return(SPEARobj)
}



