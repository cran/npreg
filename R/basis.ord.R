basis.ord <-
  function(x, knots, K = NULL, intercept = FALSE, ridge = FALSE){
    # Ordinal Smoothing Spline Basis
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Update: 2021-04-09
    
    if(is.null(K)) K <- length(unique(x))
    x <- as.ordered(x)
    xlev <- levels(x)
    knots <- factor(knots, levels = xlev, ordered = TRUE)
    x <- as.integer(x)
    knots <- as.integer(knots)
    const <- (K - 1) * (2 * K - 1) / (6 * K)
    X <- outer(X = x, Y = knots, FUN = "ordkern", K = K, const = const) / K
    if(ridge){
      Q <- penalty.ord(knots, K = K, xlev = 1:length(xlev))
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

ordkern <- 
  function(x, y, K, const){
    1 - pmax(x, y) + (x * (x - 1) + y * (y - 1)) / (2 * K) + const
  }
