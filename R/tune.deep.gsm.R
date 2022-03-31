tune.deep.gsm <-
  function(lambdas, spar, y, Etab, rkhs, weights, tprk = TRUE, method = "GCV", 
           family = check_family(gaussian), control = control){
    # deep tuning for generalized smooth model
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2022-03-22
    
    
    #########***#########   DESIGN AND PENALTY   #########***#########
    
    # define thetas from lambdas = log(thetas)
    thetas <- exp(lambdas)
    
    # build design and penalty
    depe <- build_depe(Etab = Etab, rkhs = rkhs, tprk = tprk, thetas = thetas)
    depe$weights <- weights
    
    
    tryCatch({
      
      #########***#########   INITIALIZATIONS   #########***#########
      
      # info
      nobs <- nrow(depe$K)
      nsdim <- ncol(depe$K)
      if(!tprk){
        Nknots <- sapply(depe$J, ncol)
        depe$J <- do.call(cbind, depe$J)
      }
      nknots <- ncol(depe$J)
      nullindx <- 1:nsdim
      
      # reparameterize contrast space
      if(tprk){
        
        Qisqrt <- msqrt(depe$Q, inverse = TRUE, checkx = FALSE)
        Rmat <- depe$J %*% Qisqrt
        Qrnk <- ncol(Qisqrt)
        
      } else {
        
        cknots <- c(0, cumsum(Nknots))
        Rmat <- Qisqrt <- vector("list", length(depe$Q))
        for(k in 1:length(depe$Q)){
          indx <- seq(cknots[k] + 1, cknots[k+1])
          Qisqrt[[k]] <- msqrt(depe$Q[[k]], inverse = TRUE, checkx = FALSE)
          Rmat[[k]] <- depe$J[,indx] %*% Qisqrt[[k]]
        }
        Rmat <- do.call("cbind", Rmat)
        Qrnk <- ncol(Rmat)
        
      } # end if(tprk)
      
      # reverse transformation
      Tmat <- matrix(0, nsdim + nknots, nsdim + Qrnk)
      Tmat[nullindx,nullindx] <- diag(nsdim)
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
      
      
      #########***#########   ESTIMATE COEFS   #########***#########
      
      # get initial beta0
      beta0 <- family$linkfun(mean(y))
      
      # evaluate tuning criterion
      tune.gsm(spar = spar, y = y, Kmat = depe$K, Rmat = Rmat,
               weights = depe$weights, beta0 = beta0, tprk = tprk, 
               control = control, family = family, method = method)
      
    }, error = function(e) .Machine$double.xmax)
    
  } # end tune.deep.gsm