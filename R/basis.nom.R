basis.nom <- 
  function(x, knots, K = NULL, intercept = FALSE, ridge = FALSE){
    # Nominal Smoothing Spline Basis
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Update: 2022-03-22
    
    if(is.null(K)) K <- length(unique(x))
    X <- outer(x, knots, FUN = "==") - 1/K
    if(ridge){
      Q <- penalty.nom(knots, K = K)
      X <- X %*% msqrt(Q, inverse = TRUE, checkx = FALSE)
    }
    knot.names <- paste("knot", 1:ncol(X), sep = ".")
    if(intercept) {
      X <- cbind(1, X)
      knot.names <- c("(Intercept)", knot.names)
    }
    colnames(X) <- knot.names
    X
    
  }