rowkro <-
  function(X, Y){
    # rowwise Kronecker product
    
    if(is.null(X) | is.null (Y)) return(NULL)
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    xd <- dim(X)
    yd <- dim(Y)
    if(xd[1] != yd[1]) stop("Inputs 'X' and 'Y' must be matrices with the same number of rows.")
    Z <- matrix(0, xd[1], xd[2] * yd[2])
    for(j in 1:xd[2]){
      zind <- (1 + (j - 1) * yd[2]):(yd[2] + (j - 1) * yd[2])
      Z[,zind] <- X[,j] * Y
    }
    return(Z)
    
  }