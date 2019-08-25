basis_nom <- 
  function(x, knots, K = NULL, intercept = FALSE){
    # Nominal Smoothing Spline Basis
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Update: 2019-04-04
    
    if(is.null(K)) K <- length(unique(x))
    X <- outer(x, knots, FUN = "==") - 1/K
    knot.names <- paste("knot", 1:length(knots), sep = ".")
    if(intercept) {
      X <- cbind(1, X)
      knot.names <- c("(Intercept)", knot.names)
    }
    colnames(X) <- knot.names
    X
    
  }