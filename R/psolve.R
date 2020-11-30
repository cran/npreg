psolve <- 
  function(a, b, tol){
    # Pseudo-Solve a System of Equations
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Update: 2020-11-29
    
    
    ######  CHECK  ######
    
    ## check a
    if(missing(a)){
      stop("Input 'a' must be provided.")
    } else {
      a <- as.matrix(a)
      adim <- dim(a)
    }
    
    ## check b
    if(missing(b)){
      no.b <- TRUE
    } else {
      no.b <- FALSE
      b <- as.matrix(b)
      bdim <- dim(b)
      if(adim[1] != bdim[1]) stop("Inputs 'a' and 'b' must satisfy nrow(a) == nrow(b).")
    }
    
    ## check tol
    if(missing(tol)){
      tol <- .Machine$double.eps
    } else {
      tol <- as.numeric(tol[1])
      if(tol <= 0) stop("Input 'tol' must be a positive scalar.")
    }
    
    ## max dim
    mdim <- min(dim(a))
    
    
    ######  SOLVE  ######
    
    ## symmetric a?
    symmetric <- isSymmetric(unname(a))
    if(symmetric){
      aeig <- eigen(a, symmetric = TRUE)
      arnk <- sum(aeig$values > aeig$values[1] * tol * mdim)
      ainv <- tcrossprod(aeig$vectors[,1:arnk,drop=FALSE] %*% diag(1 / sqrt(aeig$values[1:arnk]), nrow = arnk, ncol = arnk))
    } else {
      asvd <- svd(a)
      arnk <- sum(asvd$d^2 > asvd$d[1]^2 * tol * mdim)
      ainv <- tcrossprod(asvd$v[,1:arnk,drop=FALSE] %*% diag(1 / asvd$d[1:arnk], nrow = arnk, ncol = arnk), asvd$u[,1:arnk,drop=FALSE])
    }
    
    ## return solution
    if(no.b){
      return(ainv)
    } else {
      return(ainv %*% b)
    }
    
  } # end psolve