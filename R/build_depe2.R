build_depe2 <- 
  function(Etab, rkhs, thetas, Jcoef, depe, tprk = TRUE){
    # rebuild design and penalty matrix
    # Nathaniel E. Helwig (helwig@umn.edu)
    # updated: 2019-04-07
    
    ### get info
    xnames <- rownames(Etab)[-1]
    nxvar <- length(xnames)
    nterms <- ncol(Etab)
    
    ### tprk (smoothing spline anova)
    if(tprk){
      
      ## get nobs
      nobs <- nrow(rkhs$Xn[[1]])
      if(is.null(nobs)) nobs <- nrow(rkhs$Xc[[1]])
      
      ## initializations
      Kmat <- matrix(1, nrow = nobs, ncol = 1)
      Knames <- "(Intercept)"
      
      ## get nknots
      nknots <- NULL
      k <- 1L
      while(is.null(nknots)){
        nknots <- ncol(rkhs$Qc[[k]])
        k <- k + 1L
      }
      
      ## initialize J and Q
      Jmat <- matrix(0, nrow = nobs, ncol = nknots)
      Qmat <- matrix(0, nrow = nknots, ncol = nknots)
      colnames(Jmat) <- rownames(Qmat) <- colnames(Qmat) <- paste("knot", 1:nknots, sep=".")
      
      ## initialize smoothing parameters
      theta0 <- thetas
      thetas <- thetaNames <- thetaAttr <- NULL
      
      ## sweep through model terms
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
          
          # add to contrast and penalty matrices
          if(!is.null(rkhs$Xc[[cidx]])){
            ssnew <- as.numeric(t(Jcoef) %*% (rkhs$Qc[[cidx]] %*% Jcoef))
            thetanew <- theta0[1]^2 * ssnew
            theta0 <- theta0[-1]
            thetas <- c(thetas, thetanew)
            thetaAttr <- c(thetaAttr, xnames[cidx])
            thetaNames <- c(thetaNames, xnames[cidx])
            Qmat <- Qmat + thetanew * rkhs$Qc[[cidx]]
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
          
          # add to contrast and penalty matrices (ii.cont and jj.null)
          if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xn[[jj]])){
            Qnew <- rkhs$Qc[[ii]] * tcrossprod(rkhs$Qn[[jj]])
            ssnew <- as.numeric(t(Jcoef) %*% (Qnew %*% Jcoef))
            thetanew <- theta0[1]^2 * ssnew
            theta0 <- theta0[-1]
            thetas <- c(thetas, thetanew)
            thetaAttr <- c(thetaAttr, paste0(xnames[ii],".c:",xnames[jj],".n"))
            thetaNames <- c(thetaNames, paste0(xnames[ii],":",xnames[jj]))
            Qmat <- Qmat + thetanew * Qnew
            Jmat <- Jmat + thetanew * rkhs$Xc[[ii]] * tcrossprod(rkhs$Xn[[jj]], rkhs$Qn[[jj]])
            rm(Qnew)
          }
          
          # add to contrast and penalty matrices (ii.null and jj.cont)
          if(!is.null(rkhs$Xn[[ii]]) && !is.null(rkhs$Xc[[jj]])){
            Qnew <- tcrossprod(rkhs$Qn[[ii]]) * rkhs$Qc[[jj]]
            ssnew <- as.numeric(t(Jcoef) %*% (Qnew %*% Jcoef))
            thetanew <- theta0[1]^2 * ssnew
            theta0 <- theta0[-1]
            thetas <- c(thetas, thetanew)
            thetaAttr <- c(thetaAttr, paste0(xnames[ii],".n:",xnames[jj],".c"))
            thetaNames <- c(thetaNames, paste0(xnames[ii],":",xnames[jj]))
            Qmat <- Qmat + thetanew * Qnew
            Jmat <- Jmat + thetanew * tcrossprod(rkhs$Xn[[ii]], rkhs$Qn[[ii]]) * rkhs$Xc[[jj]]
            rm(Qnew)
          }
          
          # add to contrast and penalty matrices (ii.cont and jj.cont)
          if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xc[[jj]])){
            Qnew <- rkhs$Qc[[ii]] * rkhs$Qc[[jj]]
            ssnew <- as.numeric(t(Jcoef) %*% (Qnew %*% Jcoef))
            thetanew <- theta0[1]^2 * ssnew
            theta0 <- theta0[-1]
            thetas <- c(thetas, thetanew)
            thetaAttr <- c(thetaAttr, paste0(xnames[ii],".c:",xnames[jj],".c"))
            thetaNames <- c(thetaNames, paste0(xnames[ii],":",xnames[jj]))
            Qmat <- Qmat + thetanew * Qnew
            Jmat <- Jmat + thetanew * rkhs$Xc[[ii]] * rkhs$Xc[[jj]]
            rm(Qnew)
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
          
          # add to contrast and penalty matrices (ii.cont and jj.null and kk.null)
          if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xn[[jj]]) && !is.null(rkhs$Xn[[kk]])){
            Qnew <- rkhs$Qc[[ii]] * tcrossprod(rkhs$Qn[[jj]]) * tcrossprod(rkhs$Qn[[kk]])
            ssnew <- as.numeric(t(Jcoef) %*% (Qnew %*% Jcoef))
            thetanew <- theta0[1]^2 * ssnew
            theta0 <- theta0[-1]
            thetas <- c(thetas, thetanew)
            thetaAttr <- c(thetaAttr, paste0(xnames[ii],".c:", xnames[jj],".n:", xnames[kk],".n"))
            thetaNames <- c(thetaNames, paste0(xnames[ii],":", xnames[jj],":", xnames[kk]))
            Qmat <- Qmat + thetanew * Qnew
            Jmat <- Jmat + thetanew * rkhs$Xc[[ii]] * tcrossprod(rkhs$Xn[[jj]], rkhs$Qn[[jj]]) * tcrossprod(rkhs$Xn[[kk]], rkhs$Qn[[kk]])
            rm(Qnew)
          }
          
          # add to contrast and penalty matrices (ii.null and jj.cont and kk.null)
          if(!is.null(rkhs$Xn[[ii]]) && !is.null(rkhs$Xc[[jj]]) && !is.null(rkhs$Xn[[kk]])){
            Qnew <- tcrossprod(rkhs$Qn[[ii]]) * rkhs$Qc[[jj]] * tcrossprod(rkhs$Qn[[kk]])
            ssnew <- as.numeric(t(Jcoef) %*% (Qnew %*% Jcoef))
            thetanew <- theta0[1]^2 * ssnew
            theta0 <- theta0[-1]
            thetas <- c(thetas, thetanew)
            thetaAttr <- c(thetaAttr, paste0(xnames[ii],".n:", xnames[jj],".c:", xnames[kk],".n"))
            thetaNames <- c(thetaNames, paste0(xnames[ii],":", xnames[jj],":", xnames[kk]))
            Qmat <- Qmat + thetanew * Qnew
            Jmat <- Jmat + thetanew * tcrossprod(rkhs$Xn[[ii]], rkhs$Qn[[ii]]) * rkhs$Xc[[jj]] * tcrossprod(rkhs$Xn[[kk]], rkhs$Qn[[kk]])
            rm(Qnew)
          }
          
          # add to contrast and penalty matrices (ii.null and jj.null and kk.cont)
          if(!is.null(rkhs$Xn[[ii]]) && !is.null(rkhs$Xn[[jj]]) && !is.null(rkhs$Xc[[kk]])){
            Qnew <- tcrossprod(rkhs$Qn[[ii]]) * tcrossprod(rkhs$Qn[[jj]]) * rkhs$Qc[[kk]]
            ssnew <- as.numeric(t(Jcoef) %*% (Qnew %*% Jcoef))
            thetanew <- theta0[1]^2 * ssnew
            theta0 <- theta0[-1]
            thetas <- c(thetas, thetanew)
            thetaAttr <- c(thetaAttr, paste0(xnames[ii],".n:", xnames[jj],".n:", xnames[kk],".c"))
            thetaNames <- c(thetaNames, paste0(xnames[ii],":", xnames[jj],":", xnames[kk]))
            Qmat <- Qmat + thetanew * Qnew
            Jmat <- Jmat + thetanew * tcrossprod(rkhs$Xn[[ii]], rkhs$Qn[[ii]]) * tcrossprod(rkhs$Xn[[jj]], rkhs$Qn[[jj]]) * rkhs$Xc[[kk]]
            rm(Qnew)
          }
          
          # add to contrast and penalty matrices (ii.cont and jj.cont and kk.null)
          if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xc[[jj]]) && !is.null(rkhs$Xn[[kk]])){
            Qnew <- rkhs$Qc[[ii]] * rkhs$Qc[[jj]] * tcrossprod(rkhs$Qn[[kk]])
            ssnew <- as.numeric(t(Jcoef) %*% (Qnew %*% Jcoef))
            thetanew <- theta0[1]^2 * ssnew
            theta0 <- theta0[-1]
            thetas <- c(thetas, thetanew)
            thetaAttr <- c(thetaAttr, paste0(xnames[ii],".c:", xnames[jj],".c:", xnames[kk],".n"))
            thetaNames <- c(thetaNames, paste0(xnames[ii],":", xnames[jj],":", xnames[kk]))
            Qmat <- Qmat + thetanew * Qnew
            Jmat <- Jmat + thetanew * rkhs$Xc[[ii]] * rkhs$Xc[[jj]] * tcrossprod(rkhs$Xn[[kk]], rkhs$Qn[[kk]])
            rm(Qnew)
          }
          
          # add to contrast and penalty matrices (ii.cont and jj.null and kk.cont)
          if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xn[[jj]]) && !is.null(rkhs$Xc[[kk]])){
            Qnew <- rkhs$Qc[[ii]] * tcrossprod(rkhs$Qn[[jj]]) * rkhs$Qc[[kk]]
            ssnew <- as.numeric(t(Jcoef) %*% (Qnew %*% Jcoef))
            thetanew <- theta0[1]^2 * ssnew
            theta0 <- theta0[-1]
            thetas <- c(thetas, thetanew)
            thetaAttr <- c(thetaAttr, paste0(xnames[ii],".c:", xnames[jj],".n:", xnames[kk],".c"))
            thetaNames <- c(thetaNames, paste0(xnames[ii],":", xnames[jj],":", xnames[kk]))
            Qmat <- Qmat + thetanew * Qnew
            Jmat <- Jmat + thetanew * rkhs$Xc[[ii]] * tcrossprod(rkhs$Xn[[jj]], rkhs$Qn[[jj]]) * rkhs$Xc[[kk]]
            rm(Qnew)
          }
          
          # add to contrast and penalty matrices (ii.null and jj.cont and kk.cont)
          if(!is.null(rkhs$Xn[[ii]]) && !is.null(rkhs$Xc[[jj]]) && !is.null(rkhs$Xc[[kk]])){
            Qnew <- tcrossprod(rkhs$Qn[[ii]]) * rkhs$Qc[[jj]] * rkhs$Qc[[kk]]
            ssnew <- as.numeric(t(Jcoef) %*% (Qnew %*% Jcoef))
            thetanew <- theta0[1]^2 * ssnew
            theta0 <- theta0[-1]
            thetas <- c(thetas, thetanew)
            thetaAttr <- c(thetaAttr, paste0(xnames[ii],".n:", xnames[jj],".c:", xnames[kk],".c"))
            thetaNames <- c(thetaNames, paste0(xnames[ii],":", xnames[jj],":", xnames[kk]))
            Qmat <- Qmat + thetanew * Qnew
            Jmat <- Jmat + thetanew * tcrossprod(rkhs$Xn[[ii]], rkhs$Qn[[ii]]) * rkhs$Xc[[jj]] * rkhs$Xc[[kk]]
            rm(Qnew)
          }
          
          # add to contrast and penalty matrices (ii.cont and jj.cont and kk.cont)
          if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xc[[jj]]) && !is.null(rkhs$Xc[[kk]])){
            Qnew <- rkhs$Qc[[ii]] * rkhs$Qc[[jj]] * rkhs$Qc[[kk]]
            ssnew <- as.numeric(t(Jcoef) %*% (Qnew %*% Jcoef))
            thetanew <- theta0[1]^2 * ssnew
            theta0 <- theta0[-1]
            thetas <- c(thetas, thetanew)
            thetaAttr <- c(thetaAttr, paste0(xnames[ii],".c:", xnames[jj],".c:", xnames[kk],".c"))
            thetaNames <- c(thetaNames, paste0(xnames[ii],":", xnames[jj],":", xnames[kk]))
            Qmat <- Qmat + thetanew * Qnew
            Jmat <- Jmat + thetanew * rkhs$Xc[[ii]] * rkhs$Xc[[jj]] * rkhs$Xc[[kk]]
            rm(Qnew)
          }
          
        } # end if(lencidx == 1L)
        
      } # end for(k in 1:nterms)
      
      colnames(Kmat) <- Knames
      names(thetas) <- thetaNames
      return(list(K = Kmat, J = Jmat, Q = Qmat, thetas = thetas))
      
    } # end if(tprk)
    
    
    ### !tprk (generalized additive model)
    if(!tprk){
      
      ## get nobs
      nobs <- nrow(depe$K)
      
      ## get nknots
      nknots <- sapply(depe$Q, ncol)
      kindex <- c(0, cumsum(nknots))
      
      ## sweep through model terms
      for(k in 1:length(depe$Q)){
        indx <- seq(kindex[k] + 1, kindex[k + 1])
        thetas[k] <- depe$thetas[k] * t(Jcoef[indx]) %*% depe$Q[[k]] %*% Jcoef[indx]
        depe$Q[[k]] <- (thetas[k] / depe$thetas[k]) * depe$Q[[k]]
        depe$J[[k]] <- (thetas[k] / depe$thetas[k]) * depe$J[[k]]
      }
      depe$thetas <- thetas
      
      ## return results
      return(depe)
      
    } # end if(!tprk)
    
  }