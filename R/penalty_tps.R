penalty_tps <-
  function(x, m = 2, rk = TRUE){
    # Thin-Plate Spline Penalty
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Update: 2020-11-28
    
    # initializations
    m <- as.integer(m)
    if(m < 1L | m > 3L) stop("Input 'm' must be 1 (linear), 2 (cubic), or 3 (quintic).")
    x <- as.matrix(x)
    nx <- nrow(x)
    xdim <- ncol(x)
    if(2 * m <= xdim) stop("Need input 'm' to satisfy:  2 * m > ncol(x)")
    even <- ifelse((xdim %% 2) == 0, TRUE, FALSE)
    
    # euclidean distance
    dmat <- as.matrix(dist(x))
    
    # tps semi-kernel
    if(even){
      top <- (-1) ^ (m + 1 + xdim/2)
      bot <- (2^(2*m - 1)) * (pi^(xdim/2)) * factorial(m-1) * factorial(m - xdim/2)
    } else {
      top <- gamma((xdim/2) - m)
      bot <- (2^(2*m)) * (pi^(xdim/2)) * factorial(m-1)
    }
    const <- top / bot
    Rkern <- const * (dmat ^ (2 * m - xdim))
    if(even){
      dmat <- log(dmat)
      naid <- which(is.infinite(dmat))
      dmat[naid] <- 0
      Rkern <- Rkern * dmat
    } 
    
    # return results
    if(!rk) return(Rkern)
    
    # X null
    if(m == 1L){
      Xnull <- matrix(1, nx)
    } else if(m == 2L){
      Xnull <- cbind(1, x)
    } else if(m == 3L){
      Xnull <- cbind(1, x)
      for(j in 1:xdim){
        for(i in 1:j) Xnull <- cbind(Xnull, x[,i] * x[,j])
      }
    }
    
    # svd of Xnull
    Xsvd <- svd(Xnull)
    
    # define reproducing kernel
    RPmat <- tcrossprod(Rkern %*% Xsvd$u, Xsvd$u)
    Rkern - RPmat - t(RPmat) + Xsvd$u %*% crossprod(Xsvd$u, RPmat)
    
  }