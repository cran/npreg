penalty_ord <-
  function(x, K = NULL, xlev = NULL){
    # Ordinal Smoothing Spline Penalty
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Update: 2019-04-04
    
    if(is.null(K)) K <- length(unique(x))
    if(is.null(xlev)){
      x <- as.ordered(x)
    } else {
      x <- factor(x, levels = xlev, ordered = TRUE)
    }
    x <- as.integer(x)
    const <- (K - 1) * (2 * K - 1) / (6 * K)
    X <- outer(X = x, Y = x, FUN = "ordkern", K = K, const = const) / K
    X
    
  }