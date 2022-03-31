knot1samp <- 
  function(x, n = NULL){
    
    # check class of 'x'
    xc <- class(x)[1]
    xclass <- c("character", "factor", "ordered", "integer", "numeric", "matrix")
    if(!any(xc == xclass)) stop(paste("Input 'x' must be one of the following object classes:\n", paste0(xclass, collapse = ", ")))
    nx <- ifelse(xc == "matrix", nrow(x), length(x))
    
    # check if 'x' is matrix with one column
    if(xc == "matrix" && ncol(x) == 1L){
      x <- as.vector(x)
      xc <- class(x)[1]
      if(!any(xc == xclass)) stop(paste("Input 'x' must be one of the following object classes:\n", paste0(xclass, collapse = ", ")))
    }
    
    # if x is a factor, return K-1 levels as knots
    if(any(xc == c("character", "factor", "ordered"))){
      if(xc == "character") x <- as.factor(x)
      xlev <- levels(x)
      nlev <- length(xlev)
      if(nlev == 1L) stop("Input 'x' is a factor with 1 level.")
      ordered <- ifelse(xc == "ordered", TRUE, FALSE)
      return(factor(xlev[-1], levels = xlev, ordered = ordered))
    }
    
    # check 'n'
    if(is.null(n)) {
      n <- min(5L, nx)
    } else {
      n <- as.integer(n[1])
      if(n <= 2L) stop("Need n >= 3 knots.")
      if(n >= nx) stop("Need n < length(x).")
    }
    
    # get unique 'x'
    ux <- unique(x)
    if(xc == "matrix") {
      nux <- nrow(ux)
    } else {
      ux <- sort(ux)
      nux <- length(ux)
    }
    if(nux <= n) return(ux)
    
    # if x is a vector, return quantiles as knots
    if(any(xc == c("integer", "numeric"))){
      #probs <- seq(0, 1, length.out = n + 1)[1:n]
      #return(quantile(x, probs = probs))
      probs <- seq(0, 1, length.out = n)
      return(quantile(x, probs = probs, type = 2))
    }
    
    # if x is a matrix, return bin-sampled knots
    return(bin.sample(x, n, breaks.return = TRUE)$bx)
    
  }