#' Update the dimension names for all SPEARobject matrices. Used internally.
#'@export
update.dimnames = function(call = "fit"){
  if(call == "fit"){
    dimnames(self$fit$regression.coefs) = list(colnames(self$data$train$X), paste0("Factor", 1:self$params$num_factors), paste0("w.idx.", 1:nrow(self$params$weights)))
    dimnames(self$fit$projection.coefs.x) = list(colnames(self$data$train$X), paste0("Factor", 1:self$params$num_factors), paste0("w.idx.", 1:nrow(self$params$weights))) 
    dimnames(self$fit$projection.coefs.y) = list(paste0("Factor", 1:self$params$num_factors), colnames(self$data$train$Y), paste0("w.idx.", 1:nrow(self$params$weights))) 
    dimnames(self$fit$nonzero.probs) = list(colnames(self$data$train$X), paste0("Factor", 1:self$params$num_factors), paste0("w.idx.", 1:nrow(self$params$weights))) 
    dimnames(self$fit$projection.probs) = list(colnames(self$data$train$X), paste0("Factor", 1:self$params$num_factors), paste0("w.idx.", 1:nrow(self$params$weights))) 
    dimnames(self$fit$marginal.probs) = list(colnames(self$data$train$X), paste0("Factor", 1:self$params$num_factors), paste0("w.idx.", 1:nrow(self$params$weights))) 
    dimnames(self$fit$joint.probs) = list(colnames(self$data$train$X), paste0("Factor", 1:self$params$num_factors), paste0("w.idx.", 1:nrow(self$params$weights)))
  }
}



#' Add a dataset to a SPEARobject. Needs to have the same number of datasets/columns as `SPEARobject$data$train$Xlist`. Response not required (will be recorded as NULL). Can be accessed via `SPEARobject$data$__name__`.
#' @param X List of explanatory dataset matrices. Needs to have same column names and length (number of matrices) as SPEARobject$data$Xlist.
#' @param Y Response matrix. Needs to have the same number of rows/length as the number of rows in parameter `X`. Defaults to `NULL` (not required).
#' @param Z Full matrix with all data. Cannot have missing values; must be imputed. If left blank (`NULL`) will be equivalent to `do.call("cbind", X)` after imputation.
#' @param name Name for the stored dataset in the SPEARobject. Access via `SPEARobject$data$__name__`.
#' @examples
#' SPEARobj <- make_spear_model(...)
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
      cat("Setting current weight index to ", w.idx, "\n", 
          private$color.text("w.x: ", "green"), self$params$weights[w.idx,1], "\n",
          private$color.text("w.y: ", "green"), self$params$weights[w.idx,2], "\n",
          sep = "")
      self$options$current.weight.idx = w.idx
    }
    
  } else if(!is.null(method)){
    cat("To do...")
  }
}



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