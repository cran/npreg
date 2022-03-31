tune.deep.sm <-
  function(lambdas, spar, y, Etab, rkhs, weights, tprk = TRUE, method = "GCV"){
    # deep tuning for smooth model
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2022-03-22
    
    
    #########***#########   DESIGN AND PENALTY   #########***#########
    
    # define thetas from lambdas = log(thetas)
    thetas <- exp(lambdas)
    
    # build design and penalty
    depe <- build_depe(Etab = Etab, rkhs = rkhs, tprk = tprk, thetas = thetas)
    depe$weights <- weights
    
    
    tryCatch({
      
      #########***#########   REPARAMETERIZATION   #########***#########
      
      # info
      n <- nrow(depe$K)
      nsdim <- ncol(depe$K)
      if(!tprk){
        Nknots <- sapply(depe$J, ncol)
        depe$J <- do.call(cbind, depe$J)
      }
      nknots <- ncol(depe$J)
      nullindx <- 1:nsdim
      
      # weighted matrices
      wsqrt <- sqrt(depe$weights)
      y.w <- y * wsqrt
      yss <- sum(y.w^2)
      depe$K <- depe$K * wsqrt
      depe$J <- depe$J * wsqrt
      
      # reparameterize contrast space
      if(tprk){
        
        Qisqrt <- msqrt(depe$Q, inverse = TRUE, checkx = FALSE)
        R.w <- depe$J %*% Qisqrt
        Qrnk <- ncol(Qisqrt)
        
      } else {
        
        cknots <- c(0, cumsum(Nknots))
        R.w <- Qisqrt <- vector("list", length(depe$Q))
        for(k in 1:length(depe$Q)){
          indx <- seq(cknots[k] + 1, cknots[k+1])
          Qisqrt[[k]] <- msqrt(depe$Q[[k]], inverse = TRUE, checkx = FALSE)
          R.w[[k]] <- depe$J[,indx] %*% Qisqrt[[k]]
        }
        R.w <- do.call("cbind", R.w)
        Qrnk <- ncol(R.w)
        
      } # end if(tprk)
      
      # define X
      XsvdN <- svd(depe$K)
      X.w <- cbind(depe$K, R.w - XsvdN$u %*% crossprod(XsvdN$u, R.w))
      
      # sse for null space
      beta0 <- crossprod(XsvdN$u, y.w)
      fit0 <- as.numeric(XsvdN$u %*% beta0)
      sse0 <- yss - 2 * sum(y.w * fit0) + sum(beta0^2)
      lev0 <- rowSums(XsvdN$u^2)
      
      # SVD of weighted X
      XsvdC <- svd(X.w[,-nullindx,drop=FALSE])
      bvec <- crossprod(XsvdC$u, y.w)
      dvec <- 1/XsvdC$d^2
      
      # SVD of ML or REML
      if(method == "REML") {
        avec <- XsvdC$d^2
        nval <- sum(log(XsvdN$d^2))
        mliw <- mean(log(1/depe$weights[depe$weights > 0]))
        const <- mliw + (nval / n) + (n - nsdim) * (1 + log(2 * pi) - log(n - nsdim) ) / n
      } else if (method == "ML"){
        avec <- svd(R.w, nu = 0, nv = 0)$d^2
        mliw <- mean(log(1/depe$weights[depe$weights > 0]))
        const <- mliw + 1 + log(2*pi) - log(n)
      }
      
      # reverse transformation
      Tmat <- matrix(0, nsdim + nknots, nsdim + Qrnk)
      Tmat[nullindx,nullindx] <- diag(nsdim)
      Tmat[nullindx,-nullindx] <- (-1) * solve(crossprod(X.w[,nullindx])) %*% crossprod(X.w[,nullindx], R.w)
      if(tprk){
        Tmat[-nullindx,-nullindx] <- Qisqrt
      } else {
        row.offset <- col.offset <- nsdim
        for(k in 1:length(Qisqrt)){
          nrowk <- nrow(Qisqrt[[k]])
          ncolk <- ncol(Qisqrt[[k]])
          Tmat[row.offset + 1:nrowk, col.offset + 1:ncolk] <- Qisqrt[[k]]
          row.offset <- row.offset + nrowk
          col.offset <- col.offset + ncolk
        }
      }
      Tmat[,nullindx] <- Tmat[,nullindx] %*% XsvdN$v %*% diag(1 / XsvdN$d, nrow = nsdim, ncol = nsdim)
      #Tmat[,-nullindx] <- Tmat[,-nullindx] %*% XsvdC$v %*% diag(1 / XsvdC$d)
      Tmat <- cbind(Tmat[,nullindx], Tmat[,-nullindx] %*% XsvdC$v %*% diag(1 / XsvdC$d))
      
      # get coefficient names
      coefnames <- c(colnames(depe$K), colnames(depe$J))
      
      # remove junk
      Q <- depe$Q
      rm(depe, Qisqrt, Qrnk, X.w, R.w)
      
      
      #########***#########   EVALUATE CRITERION   #########***#########
      
      df.offset <- 0
      penalty <- 1
      if(method == "GCV"){
        crit <- tune.gcv.ss(spar = spar, bvec = bvec, dvec = dvec, n = n, 
                            yss = sse0, nsdim = nsdim, df.offset = df.offset, 
                            penalty = penalty)
      } else if(method == "OCV"){
        crit <- tune.ocv.ss(spar = spar, bvec = bvec, dvec = dvec, n = n, 
                            u = XsvdC$u, y = y.w, fit0 = fit0, lev0 = lev0)
      } else if(method == "GACV"){
        crit <- tune.gacv.ss(spar = spar, bvec = bvec, dvec = dvec, n = n, 
                             yss = sse0, nsdim = nsdim, const = yss/2)
      } else if(method == "ACV"){
        crit <- tune.acv.ss(spar = spar, bvec = bvec, dvec = dvec, n = n, 
                            yss = sse0, u = XsvdC$u, y = y.w, fit0 = fit0,
                            lev0 = lev0, const = yss/2)
      } else if(any(method == c("REML", "ML"))){
        crit <- tune.mle.ss(spar = spar, avec = avec, bvec = bvec, dvec = dvec, 
                            n = n, yss = yss - sum(y.w * fit0), 
                            m = ifelse(method == "REML", nsdim, 0), const = const)
      } else if(any(method == c("AIC", "BIC"))){
        crit <- tune.aic.ss(spar = spar, bvec = bvec, dvec = dvec, n = n, 
                            yss = sse0, nsdim = nsdim, 
                            const = ifelse(method == "AIC", 2, log(n)))
      } # end if(method == "GCV")
      return(crit)
      
    }, error = function(e) .Machine$double.xmax)
    
  } # end tune.deep.sm