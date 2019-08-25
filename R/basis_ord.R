basis_ord <-
  function(x, knots, K = NULL, intercept = FALSE){
    # Ordinal Smoothing Spline Basis
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Update: 2019-04-04
    
    if(is.null(K)) K <- length(unique(x))
    x <- as.ordered(x)
    knots <- factor(knots, levels = levels(x), ordered = TRUE)
    x <- as.integer(x)
    knots <- as.integer(knots)
    const <- (K - 1) * (2 * K - 1) / (6 * K)
    X <- outer(X = x, Y = knots, FUN = "ordkern", K = K, const = const) / K
    knot.names <- paste("knot", 1:length(knots), sep = ".")
    if(intercept) {
      X <- cbind(1, X)
      knot.names <- c("(Intercept)", knot.names)
    }
    colnames(X) <- knot.names
    X
    
  }

ordkern <- 
  function(x, y, K, const){
    1 - pmax(x, y) + (x * (x - 1) + y * (y - 1)) / (2 * K) + const
  }
