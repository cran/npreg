pred_depe <- 
  function(Etab, rkhs, tprk, thetas){
    # build design matrix for prediction
    # Nathaniel E. Helwig (helwig@umn.edu)
    # updated: 2019-04-28
    
    ### get info
    xnames <- rownames(Etab)[-1]
    nxvar <- length(xnames)
    nterms <- ncol(Etab)
    
    ### get nobs
    nobs <- nrow(rkhs$Xn[[1]])
    if(is.null(nobs)) nobs <- nrow(rkhs$Xc[[1]])
    
    ### initializations
    Kmat <- matrix(1, nrow = nobs, ncol = 1)
    Knames <- "(Intercept)"
    
    ### check rkhs$Xc
    if(nobs == 1L){
      for(k in 1:length(rkhs$Xc)){
        if(!is.null(rkhs$Xc[[k]])) rkhs$Xc[[k]] <- matrix(rkhs$Xc[[k]], nrow = 1)
      }
    }
    
    ### tprk (smoothing spline anova)
    if(tprk){
      
      # get nknots
      nknots <- NULL
      k <- 1L
      while(is.null(nknots)){
        nknots <- ncol(rkhs$Xc[[k]])
        k <- k + 1L
      }
      
      # initialize J
      Jmat <- matrix(0, nrow = nobs, ncol = nknots)
      colnames(Jmat) <- paste("knot", 1:nknots, sep=".")
      
      # sweep through model terms
      for(k in 1:nterms){
        
        # which predictor(s)?
        cidx <- which(Etab[,k] > 0L) - 1L
        lencidx <- length(cidx)
        
        # main or interaction effect?
        if(lencidx == 1L){
          # one-way (main effect) term
          
          # add to null space matrix
          if(!is.null(rkhs$Xn[[cidx]])){
            nk <- ncol(as.matrix(rkhs$Xn[[cidx]]))
            Kmat <- cbind(Kmat, rkhs$Xn[[cidx]])
            newnames <- colnames(rkhs$Xn[[cidx]])
            if(is.null(newnames)) stop("missing names for null space component")
            Knames <- c(Knames, newnames)
          }
          
          # add to contrast space matrix
          if(!is.null(rkhs$Xc[[cidx]])){
            thetanew <- thetas[1]
            thetas <- thetas[-1]
            Jmat <- Jmat + thetanew * rkhs$Xc[[cidx]]
          }
          
        } else if(lencidx == 2L){
          # two-way interaction effect
          
          ii <- cidx[1]
          jj <- cidx[2]
          
          # add to null space matrix
          Knew <- rowkro(rkhs$Xn[[ii]], rkhs$Xn[[jj]])
          if(!is.null(Knew)){
            Kmat <- cbind(Kmat, Knew)
            newnames <- charkro(colnames(rkhs$Xn[[ii]]), colnames(rkhs$Xn[[jj]]))
            if(is.null(newnames)) stop("missing names for null space component")
            Knames <- c(Knames, newnames)
          }
          
          # add to contrast space matrix (ii.cont and jj.null)
          if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xn[[jj]])){
            thetanew <- thetas[1]
            thetas <- thetas[-1]
            Jmat <- Jmat + thetanew * rkhs$Xc[[ii]] * tcrossprod(rkhs$Xn[[jj]], rkhs$Qn[[jj]])
          }
          
          # add to contrast space matrix (ii.null and jj.cont)
          if(!is.null(rkhs$Xn[[ii]]) && !is.null(rkhs$Xc[[jj]])){
            thetanew <- thetas[1]
            thetas <- thetas[-1]
            Jmat <- Jmat + thetanew * tcrossprod(rkhs$Xn[[ii]], rkhs$Qn[[ii]]) * rkhs$Xc[[jj]]
          }
          
          # add to contrast space matrix (ii.cont and jj.cont)
          if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xc[[jj]])){
            thetanew <- thetas[1]
            thetas <- thetas[-1]
            Jmat <- Jmat + thetanew * rkhs$Xc[[ii]] * rkhs$Xc[[jj]]
          }
          
        } else if(lencidx == 3L){
          # three-way interaction effect
          
          ii <- cidx[1]
          jj <- cidx[2]
          kk <- cidx[3]
          
          # add to null space matrix
          if(!any(is.null(rkhs$Xn[[ii]]), is.null(rkhs$Xn[[jj]]), is.null(rkhs$Xn[[kk]]))){
            Kmat <- cbind(Kmat, rowkro(rowkro(rkhs$Xn[[ii]], rkhs$Xn[[jj]]), rkhs$Xn[[kk]]))
            newnames <- charkro(charkro(colnames(rkhs$Xn[[ii]]), colnames(rkhs$Xn[[jj]])), colnames(rkhs$Xn[[kk]]))
            if(is.null(newnames)) stop("missing names for null space component")
            Knames <- c(Knames, newnames)
          }
          
          # add to contrast space matrix (ii.cont and jj.null and kk.null)
          if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xn[[jj]]) && !is.null(rkhs$Xn[[kk]])){
            thetanew <- thetas[1]
            thetas <- thetas[-1]
            Jmat <- Jmat + thetanew * rkhs$Xc[[ii]] * tcrossprod(rkhs$Xn[[jj]], rkhs$Qn[[jj]]) * tcrossprod(rkhs$Xn[[kk]], rkhs$Qn[[kk]])
          }
          
          # add to contrast space matrix (ii.null and jj.cont and kk.null)
          if(!is.null(rkhs$Xn[[ii]]) && !is.null(rkhs$Xc[[jj]]) && !is.null(rkhs$Xn[[kk]])){
            thetanew <- thetas[1]
            thetas <- thetas[-1]
            Jmat <- Jmat + thetanew * tcrossprod(rkhs$Xn[[ii]], rkhs$Qn[[ii]]) * rkhs$Xc[[jj]] * tcrossprod(rkhs$Xn[[kk]], rkhs$Qn[[kk]])
          }
          
          # add to contrast space matrix (ii.null and jj.null and kk.cont)
          if(!is.null(rkhs$Xn[[ii]]) && !is.null(rkhs$Xn[[jj]]) && !is.null(rkhs$Xc[[kk]])){
            thetanew <- thetas[1]
            thetas <- thetas[-1]
            Jmat <- Jmat + thetanew * tcrossprod(rkhs$Xn[[ii]], rkhs$Qn[[ii]]) * tcrossprod(rkhs$Xn[[jj]], rkhs$Qn[[jj]]) * rkhs$Xc[[kk]]
          }
          
          # add to contrast space matrix (ii.cont and jj.cont and kk.null)
          if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xc[[jj]]) && !is.null(rkhs$Xn[[kk]])){
            thetanew <- thetas[1]
            thetas <- thetas[-1]
            Jmat <- Jmat + thetanew * rkhs$Xc[[ii]] * rkhs$Xc[[jj]] * tcrossprod(rkhs$Xn[[kk]], rkhs$Qn[[kk]])
          }
          
          # add to contrast space matrix (ii.cont and jj.null and kk.cont)
          if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xn[[jj]]) && !is.null(rkhs$Xc[[kk]])){
            thetanew <- thetas[1]
            thetas <- thetas[-1]
            Jmat <- Jmat + thetanew * rkhs$Xc[[ii]] * tcrossprod(rkhs$Xn[[jj]], rkhs$Qn[[jj]]) * rkhs$Xc[[kk]]
          }
          
          # add to contrast space matrix (ii.null and jj.cont and kk.cont)
          if(!is.null(rkhs$Xn[[ii]]) && !is.null(rkhs$Xc[[jj]]) && !is.null(rkhs$Xc[[kk]])){
            thetanew <- thetas[1]
            thetas <- thetas[-1]
            Jmat <- Jmat + thetanew * tcrossprod(rkhs$Xn[[ii]], rkhs$Qn[[ii]]) * rkhs$Xc[[jj]] * rkhs$Xc[[kk]]
          }
          
          # add to contrast space matrix (ii.cont and jj.cont and kk.cont)
          if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xc[[jj]]) && !is.null(rkhs$Xc[[kk]])){
            thetanew <- thetas[1]
            thetas <- thetas[-1]
            Jmat <- Jmat + thetanew * rkhs$Xc[[ii]] * rkhs$Xc[[jj]] * rkhs$Xc[[kk]]
          }
          
        } # end if(lencidx == 1L)
        
      } # end for(k in 1:nterms)
      
      colnames(Kmat) <- Knames
      return(list(K = Kmat, J = Jmat))
      
    } # end if(tprk)
    
    
    ### !tprk (generalized additive model)
    if(!tprk){
      
      # initialize J
      Jmat <- NULL
      Jnames <- NULL
      
      # sweep through model terms
      for(k in 1:nterms){
        
        # which predictor(s)?
        cidx <- which(Etab[,k] > 0L) - 1L
        lencidx <- length(cidx)
        
        # main or interaction effect?
        if(lencidx == 1L){
          # one-way (main effect) term
          
          # add to null space matrix
          if(!is.null(rkhs$Xn[[cidx]])){
            nk <- ncol(as.matrix(rkhs$Xn[[cidx]]))
            Kmat <- cbind(Kmat, rkhs$Xn[[cidx]])
            newnames <- colnames(rkhs$Xn[[cidx]])
            if(is.null(newnames)) stop("missing names for null space component")
            Knames <- c(Knames, newnames)
          }
          
          # add to contrast space matrix
          if(!is.null(rkhs$Xc[[cidx]])){
            thetanew <- thetas[1]
            thetas <- thetas[-1]
            Jnew <- rkhs$Xc[[cidx]]
            newnames <- paste(xnames[cidx], colnames(Jnew), sep = ".")
            if(is.null(newnames)) stop("missing names for contrast space component")
            Jnames <- c(Jnames, newnames)
            colnames(Jnew) <- newnames
            Jmat <- cbind(Jmat, thetanew * Jnew)
            rm(Jnew)
          }
          
        } else if(lencidx == 2L){
          # two-way interaction effect
          
          ii <- cidx[1]
          jj <- cidx[2]
          
          # add to null space matrix
          Knew <- rowkro(rkhs$Xn[[ii]], rkhs$Xn[[jj]])
          if(!is.null(Knew)){
            Kmat <- cbind(Kmat, Knew)
            newnames <- charkro(colnames(rkhs$Xn[[ii]]), colnames(rkhs$Xn[[jj]]))
            if(is.null(newnames)) stop("missing names for null space component")
            Knames <- c(Knames, newnames)
          }
          
          # add to contrast space matrix (ii.cont and jj.null)
          if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xn[[jj]])){
            thetanew <- thetas[1]
            thetas <- thetas[-1]
            newnames <- charkro(paste(xnames[ii], colnames(rkhs$Xc[[ii]]), sep = "."), 
                                colnames(rkhs$Xn[[jj]]))
            if(is.null(newnames)) stop("missing names for contrast space component")
            Jnames <- c(Jnames, newnames)
            Jnew <- rowkro(rkhs$Xc[[ii]], rkhs$Xn[[jj]])
            colnames(Jnew) <- newnames
            Jmat <- cbind(Jmat, thetanew * Jnew)
            rm(Jnew)
          }
          
          # add to contrast space matrix (ii.null and jj.cont)
          if(!is.null(rkhs$Xn[[ii]]) && !is.null(rkhs$Xc[[jj]])){
            thetanew <- thetas[1]
            thetas <- thetas[-1]
            newnames <- charkro(colnames(rkhs$Xn[[ii]]), 
                                paste(xnames[jj], colnames(rkhs$Xc[[jj]]), sep = "."))
            if(is.null(newnames)) stop("missing names for contrast space component")
            Jnames <- c(Jnames, newnames)
            Jnew <- rowkro(rkhs$Xn[[ii]], rkhs$Xc[[jj]])
            colnames(Jnew) <- newnames
            Jmat <- cbind(Jmat, thetanew * Jnew)
            rm(Jnew)
          }
          
          # add to contrast space matrix (ii.cont and jj.cont)
          if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xc[[jj]])){
            thetanew <- thetas[1]
            thetas <- thetas[-1]
            newnames <- charkro(paste(xnames[ii], colnames(rkhs$Xc[[ii]]), sep = "."), 
                                paste(xnames[jj], colnames(rkhs$Xc[[jj]]), sep = "."))
            if(is.null(newnames)) stop("missing names for contrast space component")
            Jnames <- c(Jnames, newnames)
            Jnew <- rowkro(rkhs$Xc[[ii]], rkhs$Xc[[jj]])
            colnames(Jnew) <- newnames
            Jmat <- cbind(Jmat, thetanew * Jnew)
            rm(Jnew)
          }
          
        } else if(lencidx == 3L){
          # three-way interaction effect
          
          ii <- cidx[1]
          jj <- cidx[2]
          kk <- cidx[3]
          
          # add to null space matrix
          if(!any(is.null(rkhs$Xn[[ii]]), is.null(rkhs$Xn[[jj]]), is.null(rkhs$Xn[[kk]]))){
            Kmat <- cbind(Kmat, rowkro(rowkro(rkhs$Xn[[ii]], rkhs$Xn[[jj]]), rkhs$Xn[[kk]]))
            newnames <- charkro(charkro(colnames(rkhs$Xn[[ii]]), colnames(rkhs$Xn[[jj]])), colnames(rkhs$Xn[[kk]]))
            if(is.null(newnames)) stop("missing names for null space component")
            Knames <- c(Knames, newnames)
          }
          
          # add to contrast space matrix (ii.cont and jj.null and kk.null)
          if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xn[[jj]]) && !is.null(rkhs$Xn[[kk]])){
            thetanew <- thetas[1]
            thetas <- thetas[-1]
            newnames <- charkro(charkro(paste(xnames[ii], colnames(rkhs$Xc[[ii]]), sep = "."), 
                                        colnames(rkhs$Xn[[jj]])), 
                                colnames(rkhs$Xn[[kk]]))
            if(is.null(newnames)) stop("missing names for contrast space component")
            Jnames <- c(Jnames, newnames)
            Jnew <- rowkro(rowkro(rkhs$Xc[[ii]], rkhs$Xn[[jj]]), rkhs$Xn[[kk]])
            colnames(Jnew) <- newnames
            Jmat <- cbind(Jmat, thetanew * Jnew)
            rm(Jnew)
          }
          
          # add to contrast space matrix (ii.null and jj.cont and kk.null)
          if(!is.null(rkhs$Xn[[ii]]) && !is.null(rkhs$Xc[[jj]]) && !is.null(rkhs$Xn[[kk]])){
            thetanew <- thetas[1]
            thetas <- thetas[-1]
            newnames <- charkro(charkro(colnames(rkhs$Xn[[ii]]), 
                                        paste(xnames[jj], colnames(rkhs$Xc[[jj]]), sep = ".")), 
                                colnames(rkhs$Xn[[kk]]))
            if(is.null(newnames)) stop("missing names for contrast space component")
            Jnames <- c(Jnames, newnames)
            Jnew <- rowkro(rowkro(rkhs$Xn[[ii]], rkhs$Xc[[jj]]), rkhs$Xn[[kk]])
            colnames(Jnew) <- newnames
            Jmat <- cbind(Jmat, thetanew * Jnew)
            rm(Jnew)
          }
          
          # add to contrast space matrix (ii.null and jj.null and kk.cont)
          if(!is.null(rkhs$Xn[[ii]]) && !is.null(rkhs$Xn[[jj]]) && !is.null(rkhs$Xc[[kk]])){
            thetanew <- thetas[1]
            thetas <- thetas[-1]
            newnames <- charkro(charkro(colnames(rkhs$Xn[[ii]]), 
                                        colnames(rkhs$Xn[[jj]])), 
                                paste(xnames[kk], colnames(rkhs$Xc[[kk]]), sep = "."))
            if(is.null(newnames)) stop("missing names for contrast space component")
            Jnames <- c(Jnames, newnames)
            Jnew <- rowkro(rowkro(rkhs$Xn[[ii]], rkhs$Xn[[jj]]), rkhs$Xc[[kk]])
            colnames(Jnew) <- newnames
            Jmat <- cbind(Jmat, thetanew * Jnew)
            rm(Jnew)
          }
          
          # add to contrast space matrix (ii.cont and jj.cont and kk.null)
          if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xc[[jj]]) && !is.null(rkhs$Xn[[kk]])){
            thetanew <- thetas[1]
            thetas <- thetas[-1]
            newnames <- charkro(charkro(paste(xnames[ii], colnames(rkhs$Xc[[ii]]), sep = "."), 
                                        paste(xnames[jj], colnames(rkhs$Xc[[jj]]), sep = ".")), 
                                colnames(rkhs$Xn[[kk]]))
            if(is.null(newnames)) stop("missing names for contrast space component")
            Jnames <- c(Jnames, newnames)
            Jnew <- rowkro(rowkro(rkhs$Xc[[ii]], rkhs$Xc[[jj]]), rkhs$Xn[[kk]])
            colnames(Jnew) <- newnames
            Jmat <- cbind(Jmat, thetanew * Jnew)
            rm(Jnew)
          }
          
          # add to contrast space matrix (ii.cont and jj.null and kk.cont)
          if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xn[[jj]]) && !is.null(rkhs$Xc[[kk]])){
            thetanew <- thetas[1]
            thetas <- thetas[-1]
            newnames <- charkro(charkro(paste(xnames[ii], colnames(rkhs$Xc[[ii]]), sep = "."), 
                                        colnames(rkhs$Xn[[jj]])), 
                                paste(xnames[kk], colnames(rkhs$Xc[[kk]]), sep = "."))
            if(is.null(newnames)) stop("missing names for contrast space component")
            Jnames <- c(Jnames, newnames)
            Jnew <- rowkro(rowkro(rkhs$Xc[[ii]], rkhs$Xn[[jj]]), rkhs$Xc[[kk]])
            colnames(Jnew) <- newnames
            Jmat <- cbind(Jmat, thetanew * Jnew)
            rm(Jnew)
          }
          
          # add to contrast space matrix (ii.null and jj.cont and kk.cont)
          if(!is.null(rkhs$Xn[[ii]]) && !is.null(rkhs$Xc[[jj]]) && !is.null(rkhs$Xc[[kk]])){
            thetanew <- thetas[1]
            thetas <- thetas[-1]
            newnames <- charkro(charkro(colnames(rkhs$Xn[[ii]]), 
                                        paste(xnames[jj], colnames(rkhs$Xc[[jj]]), sep = ".")), 
                                paste(xnames[kk], colnames(rkhs$Xc[[kk]]), sep = "."))
            if(is.null(newnames)) stop("missing names for contrast space component")
            Jnames <- c(Jnames, newnames)
            Jnew <- rowkro(rowkro(rkhs$Xn[[ii]], rkhs$Xc[[jj]]), rkhs$Xc[[kk]])
            colnames(Jnew) <- newnames
            Jmat <- cbind(Jmat, thetanew * Jnew)
            rm(Jnew)
          }
          
          # add to contrast space matrix (ii.cont and jj.cont and kk.cont)
          if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xc[[jj]]) && !is.null(rkhs$Xc[[kk]])){
            thetanew <- thetas[1]
            thetas <- thetas[-1]
            newnames <- charkro(charkro(paste(xnames[ii], colnames(rkhs$Xc[[ii]]), sep = "."), 
                                        paste(xnames[jj], colnames(rkhs$Xc[[jj]]), sep = ".")), 
                                paste(xnames[kk], colnames(rkhs$Xc[[kk]]), sep = "."))
            if(is.null(newnames)) stop("missing names for contrast space component")
            Jnames <- c(Jnames, newnames)
            Jnew <- rowkro(rowkro(rkhs$Xc[[ii]], rkhs$Xc[[jj]]), rkhs$Xc[[kk]])
            colnames(Jnew) <- newnames
            Jmat <- cbind(Jmat, thetanew * Jnew)
            rm(Jnew)
          }
          
        } # end if(lencidx == 1L)
        
      } # end for(k in 1:nterms)
      
      colnames(Kmat) <- Knames
      return(list(K = Kmat, J = Jmat))
      
    } # end if(!tprk)
    
  }