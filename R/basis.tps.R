basis.tps <-
  function(x, knots, m = 2, rk = TRUE, intercept = FALSE, ridge = FALSE){
    # Thin-Plate Spline Basis
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Update: 2021-04-09
    
    # initializations
    m <- as.integer(m)
    if(m < 1L | m > 3L) stop("Input 'm' must be 1 (linear), 2 (cubic), or 3 (quintic).")
    x <- as.matrix(x)
    knots <- as.matrix(knots)
    nx <- nrow(x)
    nk <- nrow(knots)
    xdim <- ncol(x)
    if(ncol(knots) != xdim) stop("Need inputs 'x' and 'knots' to satisfy:  ncol(x) == ncol(knots)")
    if(2 * m <= xdim) stop("Need input 'm' to satisfy:  2 * m > ncol(x)")
    if(nk < xdim) stop("Need input 'knots' to satisfy:  nrow(knots) > ncol(knots)")
    even <- ifelse((xdim %% 2) == 0, TRUE, FALSE)
    
    # euclidean distance
    dmat <- matrix(0, nx, nk)
    for(i in 1:xdim) dmat <- dmat + outer(x[,i], knots[,i], FUN = "-")^2
    dmat <- sqrt(dmat)
    
    # column names
    knot.names <- paste("knot", 1:nk, sep = ".")
    
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
    
    # null space matrix
    xnames <- colnames(x)
    if(m == 1L){
      Xnull <- NULL
      null.names <- NULL
    } else if(m == 2L){
      Xnull <- x
      if(is.null(xnames)){
        if(xdim == 1L){
          null.names <- "x"
        } else {
          null.names <- paste0("x", 1:xdim)
        }
      } else {
        null.names <- xnames
      }
    } else if(m == 3L){
      Xnull <- cbind(x, x^2)
      if(is.null(xnames)){
        if(xdim == 1L){
          null.names <- c("x", "x^2")
        } else {
          null.names <- c(paste0("x", 1:xdim), paste0("x", 1:xdim, "^2"))
        }
      } else {
        null.names <- c(xnames, paste0(xnames, "^2"))
      }
      if(xdim > 1L){
        for(j in 2:xdim){
          for(i in 1:(j-1)) {
            Xnull <- cbind(Xnull, x[,i] * x[,j])
            if(is.null(xnames)){
              null.names <- c(null.names, paste0("x", i, ":x",j))
            } else {
              null.names <- c(null.names, paste0(xnames[i], ":", xnames[j]))
            }
          } # end for(i in 1:(j-1))
        } # end for(j in 2:xdim)
      } # end if(xdim > 1L)
    } # end if(m == 1L)
    
    # return results
    if(!rk) {
      if(intercept) {
        Xnull <- cbind(matrix(1, nx), Xnull)
        null.names <- c("(Intercept)", null.names)
      }
      X <- cbind(Xnull, Rkern)
      colnames(X) <- c(null.names, knot.names)
      return(X)
    }
    
    # knots kernel
    dmatK <- as.matrix(dist(knots))
    RkernK <- const * (dmatK ^ (2 * m - xdim))
    if(even){
      dmatK <- log(dmatK)
      naid <- which(is.infinite(dmatK))
      dmatK[naid] <- 0
      RkernK <- RkernK * dmatK
    } 
    
    # knots null
    if(m == 1L){
      Knull <- matrix(1, nk)
    } else if(m == 2L){
      Knull <- cbind(1, knots)
    } else if(m == 3L){
      Knull <- cbind(1, knots, knots^2)
      if(xdim > 1L){
        for(j in 2:xdim){
          for(i in 1:(j-1)) {
            Knull <- cbind(Knull, knots[,i] * knots[,j])
          }
        }
      }
    }
    
    # svd of Knull
    Ksvd <- svd(Knull)
    
    # define reproducing kernel
    KinvR <-  Ksvd$v %*% diag(1/Ksvd$d) %*% crossprod(Ksvd$u, RkernK)
    Rkern <- Rkern - cbind(matrix(1, nx), Xnull) %*% KinvR
    Rkern <- Rkern - tcrossprod(Rkern %*% Ksvd$u, Ksvd$u)
    
    # check ridge
    if(ridge){
      Q <- penalty.tps(knots, m = m)
      Qeig <- eigen(Q, symmetric = TRUE)
      Qrnk <- sum(Qeig$values > Qeig$values[1] * ncol(Q) * .Machine$double.eps)
      Qisqrt <- Qeig$vectors[,1:Qrnk,drop=FALSE] %*% diag(1/sqrt(Qeig$values[1:Qrnk]), nrow = Qrnk, ncol = Qrnk)
      Rkern <- Rkern %*% Qisqrt
      knot.names <- knot.names[1:ncol(Rkern)]
    }
    
    # return results
    if(intercept) {
      Xnull <- cbind(matrix(1, nx), Xnull)
      null.names <- c("(Intercept)", null.names)
    }
    X <- cbind(Xnull, Rkern)
    colnames(X) <- c(null.names, knot.names)
    return(X)
    
  }