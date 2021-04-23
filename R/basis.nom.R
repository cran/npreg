basis.nom <- 
  function(x, knots, K = NULL, intercept = FALSE, ridge = FALSE){
    # Nominal Smoothing Spline Basis
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Update: 2021-04-09
    
    if(is.null(K)) K <- length(unique(x))
    X <- outer(x, knots, FUN = "==") - 1/K
    if(ridge){
      Q <- penalty.nom(knots, K = K)
      Qeig <- eigen(Q, symmetric = TRUE)
      Qrnk <- sum(Qeig$values > Qeig$values[1] * ncol(Q) * .Machine$double.eps)
      Qisqrt <- Qeig$vectors[,1:Qrnk,drop=FALSE] %*% diag(1/sqrt(Qeig$values[1:Qrnk]), nrow = Qrnk, ncol = Qrnk)
      X <- X %*% Qisqrt
    }
    knot.names <- paste("knot", 1:ncol(X), sep = ".")
    if(intercept) {
      X <- cbind(1, X)
      knot.names <- c("(Intercept)", knot.names)
    }
    colnames(X) <- knot.names
    X
    
  }