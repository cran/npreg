build_depe <- 
  function(Etab, rkhs, tprk, thetas = NULL){
    # build design and penalty matrix
    # Nathaniel E. Helwig (helwig@umn.edu)
    # updated: 2022-03-01
    
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
    
    ### thetas = NULL
    if(is.null(thetas)){
      
      ### tprk (smoothing spline anova)
      if(tprk){
        
        # get nknots
        nknots <- NULL
        k <- 1L
        while(is.null(nknots)){
          nknots <- ncol(rkhs$Qc[[k]])
          k <- k + 1L
        }
        
        # initialize J and Q
        Jmat <- matrix(0, nrow = nobs, ncol = nknots)
        Qmat <- matrix(0, nrow = nknots, ncol = nknots)
        colnames(Jmat) <- rownames(Qmat) <- colnames(Qmat) <- paste("knot", 1:nknots, sep=".")
        
        # initialize smoothing parameters
        thetas <- thetaNames <- thetaAttr <- NULL
        
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
            
            # add to contrast and penalty matrices
            if(!is.null(rkhs$Xc[[cidx]])){
              thetanew <- 1 / mean(diag(rkhs$Qc[[cidx]]))
              thetas <- c(thetas, thetanew)
              thetaNames <- c(thetaNames, xnames[cidx])
              thetaAttr <- c(thetaAttr, xnames[cidx])
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
              thetanew <- 1 / mean(diag(Qnew))
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
              thetanew <- 1 / mean(diag(Qnew))
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
              thetanew <- 1 / mean(diag(Qnew))
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
              thetanew <- 1 / mean(diag(Qnew))
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
              thetanew <- 1 / mean(diag(Qnew))
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
              thetanew <- 1 / mean(diag(Qnew))
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
              thetanew <- 1 / mean(diag(Qnew))
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
              thetanew <- 1 / mean(diag(Qnew))
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
              thetanew <- 1 / mean(diag(Qnew))
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
              thetanew <- 1 / mean(diag(Qnew))
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
        
        # initialize J and Q
        Jmat <- NULL
        Qmat <- NULL
        
        # initialize smoothing parameters
        thetas <- thetaNames <- thetaAttr <- Jnames <- NULL
        
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
            
            # add to contrast and penalty matrices
            if(!is.null(rkhs$Xc[[cidx]])){
              thetanew <- 1 / mean(diag(rkhs$Qc[[cidx]]))
              thetas <- c(thetas, thetanew)
              thetaAttr <- c(thetaAttr, xnames[cidx])
              thetaNames <- c(thetaNames, xnames[cidx])
              Qnew <- thetanew * rkhs$Qc[[cidx]]
              Jnew <- thetanew * rkhs$Xc[[cidx]]
              newnames <- paste(xnames[cidx], colnames(Jnew), sep = ".")
              colnames(Jnew) <- newnames
              if(is.null(newnames)) stop("missing names for contrast space component")
              Jnames <- c(Jnames, newnames)
              Qmat <- c(Qmat, list(Qnew))
              Jmat <- c(Jmat, list(Jnew))
              rm(Qnew, Jnew)
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
              Qnew <- kronecker(rkhs$Qc[[ii]], crossprod(rkhs$Qn[[jj]]))
              thetanew <- 1 / mean(diag(Qnew))
              thetas <- c(thetas, thetanew)
              thetaAttr <- c(thetaAttr, paste0(xnames[ii],".c:",xnames[jj],".n"))
              thetaNames <- c(thetaNames, paste0(xnames[ii],":",xnames[jj]))
              newnames <- charkro(paste(xnames[ii], colnames(rkhs$Xc[[ii]]), sep = "."), 
                                  colnames(rkhs$Xn[[jj]]))
              if(is.null(newnames)) stop("missing names for contrast space component")
              Jnames <- c(Jnames, newnames)
              Jnew <- rowkro(rkhs$Xc[[ii]], rkhs$Xn[[jj]])
              colnames(Jnew) <- colnames(Qnew) <- newnames
              Qmat <- c(Qmat, list(thetanew * Qnew))
              Jmat <- c(Jmat, list(thetanew * Jnew))
              rm(Qnew, Jnew)
            }
            
            # add to contrast and penalty matrices (ii.null and jj.cont)
            if(!is.null(rkhs$Xn[[ii]]) && !is.null(rkhs$Xc[[jj]])){
              Qnew <- kronecker(crossprod(rkhs$Qn[[ii]]), rkhs$Qc[[jj]])
              thetanew <- 1 / mean(diag(Qnew))
              thetas <- c(thetas, thetanew)
              thetaAttr <- c(thetaAttr, paste0(xnames[ii],".n:",xnames[jj],".c"))
              thetaNames <- c(thetaNames, paste0(xnames[ii],":",xnames[jj]))
              newnames <- charkro(colnames(rkhs$Xn[[ii]]), 
                                  paste(xnames[jj], colnames(rkhs$Xc[[jj]]), sep = "."))
              if(is.null(newnames)) stop("missing names for contrast space component")
              Jnames <- c(Jnames, newnames)
              Jnew <- rowkro(rkhs$Xn[[ii]], rkhs$Xc[[jj]])
              colnames(Jnew) <- colnames(Qnew) <- newnames
              Qmat <- c(Qmat, list(thetanew * Qnew))
              Jmat <- c(Jmat, list(thetanew * Jnew))
              rm(Qnew, Jnew)
            }
            
            # add to contrast and penalty matrices (ii.cont and jj.cont)
            if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xc[[jj]])){
              Qnew <- kronecker(rkhs$Qc[[ii]], rkhs$Qc[[jj]])
              thetanew <- 1 / mean(diag(Qnew))
              thetas <- c(thetas, thetanew)
              thetaAttr <- c(thetaAttr, paste0(xnames[ii],".c:",xnames[jj],".c"))
              thetaNames <- c(thetaNames, paste0(xnames[ii],":",xnames[jj]))
              newnames <- charkro(paste(xnames[ii], colnames(rkhs$Xc[[ii]]), sep = "."), 
                                  paste(xnames[jj], colnames(rkhs$Xc[[jj]]), sep = "."))
              if(is.null(newnames)) stop("missing names for contrast space component")
              Jnames <- c(Jnames, newnames)
              Jnew <- rowkro(rkhs$Xc[[ii]], rkhs$Xc[[jj]])
              colnames(Jnew) <- colnames(Qnew) <- newnames
              Qmat <- c(Qmat, list(thetanew * Qnew))
              Jmat <- c(Jmat, list(thetanew * Jnew))
              rm(Qnew, Jnew)
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
              Qnew <- kronecker(kronecker(rkhs$Qc[[ii]], crossprod(rkhs$Qn[[jj]])), crossprod(rkhs$Qn[[kk]]))
              thetanew <- 1 / mean(diag(Qnew))
              thetas <- c(thetas, thetanew)
              thetaAttr <- c(thetaAttr, paste0(xnames[ii],".c:", xnames[jj],".n:", xnames[kk],".n"))
              thetaNames <- c(thetaNames, paste0(xnames[ii],":", xnames[jj],":", xnames[kk]))
              newnames <- charkro(charkro(paste(xnames[ii], colnames(rkhs$Xc[[ii]]), sep = "."), 
                                          colnames(rkhs$Xn[[jj]])), 
                                  colnames(rkhs$Xn[[kk]]))
              if(is.null(newnames)) stop("missing names for contrast space component")
              Jnames <- c(Jnames, newnames)
              Jnew <- rowkro(rowkro(rkhs$Xc[[ii]], rkhs$Xn[[jj]]), rkhs$Xn[[kk]])
              colnames(Jnew) <- colnames(Qnew) <- newnames
              Qmat <- c(Qmat, list(thetanew * Qnew))
              Jmat <- c(Jmat, list(thetanew * Jnew))
              rm(Qnew, Jnew)
            }
            
            # add to contrast and penalty matrices (ii.null and jj.cont and kk.null)
            if(!is.null(rkhs$Xn[[ii]]) && !is.null(rkhs$Xc[[jj]]) && !is.null(rkhs$Xn[[kk]])){
              Qnew <- kronecker(kronecker(crossprod(rkhs$Qn[[ii]]), rkhs$Qc[[jj]]), crossprod(rkhs$Qn[[kk]]))
              thetanew <- 1 / mean(diag(Qnew))
              thetas <- c(thetas, thetanew)
              thetaAttr <- c(thetaAttr, paste0(xnames[ii],".n:", xnames[jj],".c:", xnames[kk],".n"))
              thetaNames <- c(thetaNames, paste0(xnames[ii],":", xnames[jj],":", xnames[kk]))
              newnames <- charkro(charkro(colnames(rkhs$Xn[[ii]]), 
                                          paste(xnames[jj], colnames(rkhs$Xc[[jj]]), sep = ".")), 
                                  colnames(rkhs$Xn[[kk]]))
              if(is.null(newnames)) stop("missing names for contrast space component")
              Jnames <- c(Jnames, newnames)
              Jnew <- rowkro(rowkro(rkhs$Xn[[ii]], rkhs$Xc[[jj]]), rkhs$Xn[[kk]])
              colnames(Jnew) <- colnames(Qnew) <- newnames
              Qmat <- c(Qmat, list(thetanew * Qnew))
              Jmat <- c(Jmat, list(thetanew * Jnew))
              rm(Qnew, Jnew)
            }
            
            # add to contrast and penalty matrices (ii.null and jj.null and kk.cont)
            if(!is.null(rkhs$Xn[[ii]]) && !is.null(rkhs$Xn[[jj]]) && !is.null(rkhs$Xc[[kk]])){
              Qnew <- kronecker(kronecker(crossprod(rkhs$Qn[[ii]]), crossprod(rkhs$Qn[[jj]])), rkhs$Qc[[kk]])
              thetanew <- 1 / mean(diag(Qnew))
              thetas <- c(thetas, thetanew)
              thetaAttr <- c(thetaAttr, paste0(xnames[ii],".n:", xnames[jj],".n:", xnames[kk],".c"))
              thetaNames <- c(thetaNames, paste0(xnames[ii],":", xnames[jj],":", xnames[kk]))
              newnames <- charkro(charkro(colnames(rkhs$Xn[[ii]]), 
                                          colnames(rkhs$Xn[[jj]])), 
                                  paste(xnames[kk], colnames(rkhs$Xc[[kk]]), sep = "."))
              if(is.null(newnames)) stop("missing names for contrast space component")
              Jnames <- c(Jnames, newnames)
              Jnew <- rowkro(rowkro(rkhs$Xn[[ii]], rkhs$Xn[[jj]]), rkhs$Xc[[kk]])
              colnames(Jnew) <- colnames(Qnew) <- newnames
              Qmat <- c(Qmat, list(thetanew * Qnew))
              Jmat <- c(Jmat, list(thetanew * Jnew))
              rm(Qnew, Jnew)
            }
            
            # add to contrast and penalty matrices (ii.cont and jj.cont and kk.null)
            if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xc[[jj]]) && !is.null(rkhs$Xn[[kk]])){
              Qnew <- kronecker(kronecker(rkhs$Qc[[ii]], rkhs$Qc[[jj]]), crossprod(rkhs$Qn[[kk]]))
              thetanew <- 1 / mean(diag(Qnew))
              thetas <- c(thetas, thetanew)
              thetaAttr <- c(thetaAttr, paste0(xnames[ii],".c:", xnames[jj],".c:", xnames[kk],".n"))
              thetaNames <- c(thetaNames, paste0(xnames[ii],":", xnames[jj],":", xnames[kk]))
              newnames <- charkro(charkro(paste(xnames[ii], colnames(rkhs$Xc[[ii]]), sep = "."), 
                                          paste(xnames[jj], colnames(rkhs$Xc[[jj]]), sep = ".")), 
                                  colnames(rkhs$Xn[[kk]]))
              if(is.null(newnames)) stop("missing names for contrast space component")
              Jnames <- c(Jnames, newnames)
              Jnew <- rowkro(rowkro(rkhs$Xc[[ii]], rkhs$Xc[[jj]]), rkhs$Xn[[kk]])
              colnames(Jnew) <- colnames(Qnew) <- newnames
              Qmat <- c(Qmat, list(thetanew * Qnew))
              Jmat <- c(Jmat, list(thetanew * Jnew))
              rm(Qnew, Jnew)
            }
            
            # add to contrast and penalty matrices (ii.cont and jj.null and kk.cont)
            if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xn[[jj]]) && !is.null(rkhs$Xc[[kk]])){
              Qnew <- kronecker(kronecker(rkhs$Qc[[ii]], crossprod(rkhs$Qn[[jj]])), rkhs$Qc[[kk]])
              thetanew <- 1 / mean(diag(Qnew))
              thetas <- c(thetas, thetanew)
              thetaAttr <- c(thetaAttr, paste0(xnames[ii],".c:", xnames[jj],".n:", xnames[kk],".c"))
              thetaNames <- c(thetaNames, paste0(xnames[ii],":", xnames[jj],":", xnames[kk]))
              newnames <- charkro(charkro(paste(xnames[ii], colnames(rkhs$Xc[[ii]]), sep = "."), 
                                          colnames(rkhs$Xn[[jj]])), 
                                  paste(xnames[kk], colnames(rkhs$Xc[[kk]]), sep = "."))
              if(is.null(newnames)) stop("missing names for contrast space component")
              Jnames <- c(Jnames, newnames)
              Jnew <- rowkro(rowkro(rkhs$Xc[[ii]], rkhs$Xn[[jj]]), rkhs$Xc[[kk]])
              colnames(Jnew) <- colnames(Qnew) <- newnames
              Qmat <- c(Qmat, list(thetanew * Qnew))
              Jmat <- c(Jmat, list(thetanew * Jnew))
              rm(Qnew, Jnew)
            }
            
            # add to contrast and penalty matrices (ii.null and jj.cont and kk.cont)
            if(!is.null(rkhs$Xn[[ii]]) && !is.null(rkhs$Xc[[jj]]) && !is.null(rkhs$Xc[[kk]])){
              Qnew <- kronecker(kronecker(crossprod(rkhs$Qn[[ii]]), rkhs$Qc[[jj]]), rkhs$Qc[[kk]])
              thetanew <- 1 / mean(diag(Qnew))
              thetas <- c(thetas, thetanew)
              thetaAttr <- c(thetaAttr, paste0(xnames[ii],".n:", xnames[jj],".c:", xnames[kk],".c"))
              thetaNames <- c(thetaNames, paste0(xnames[ii],":", xnames[jj],":", xnames[kk]))
              newnames <- charkro(charkro(colnames(rkhs$Xn[[ii]]), 
                                          paste(xnames[jj], colnames(rkhs$Xc[[jj]]), sep = ".")), 
                                  paste(xnames[kk], colnames(rkhs$Xc[[kk]]), sep = "."))
              if(is.null(newnames)) stop("missing names for contrast space component")
              Jnames <- c(Jnames, newnames)
              Jnew <- rowkro(rowkro(rkhs$Xn[[ii]], rkhs$Xc[[jj]]), rkhs$Xc[[kk]])
              colnames(Jnew) <- colnames(Qnew) <- newnames
              Qmat <- c(Qmat, list(thetanew * Qnew))
              Jmat <- c(Jmat, list(thetanew * Jnew))
              rm(Qnew, Jnew)
            }
            
            # add to contrast and penalty matrices (ii.cont and jj.cont and kk.cont)
            if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xc[[jj]]) && !is.null(rkhs$Xc[[kk]])){
              Qnew <- kronecker(kronecker(rkhs$Qc[[ii]], rkhs$Qc[[jj]]), rkhs$Qc[[kk]])
              thetanew <- 1 / mean(diag(Qnew))
              thetas <- c(thetas, thetanew)
              thetaAttr <- c(thetaAttr, paste0(xnames[ii],".c:", xnames[jj],".c:", xnames[kk],".c"))
              thetaNames <- c(thetaNames, paste0(xnames[ii],":", xnames[jj],":", xnames[kk]))
              newnames <- charkro(charkro(paste(xnames[ii], colnames(rkhs$Xc[[ii]]), sep = "."), 
                                          paste(xnames[jj], colnames(rkhs$Xc[[jj]]), sep = ".")), 
                                  paste(xnames[kk], colnames(rkhs$Xc[[kk]]), sep = "."))
              if(is.null(newnames)) stop("missing names for contrast space component")
              Jnames <- c(Jnames, newnames)
              Jnew <- rowkro(rowkro(rkhs$Xc[[ii]], rkhs$Xc[[jj]]), rkhs$Xc[[kk]])
              colnames(Jnew) <- colnames(Qnew) <- newnames
              Qmat <- c(Qmat, list(thetanew * Qnew))
              Jmat <- c(Jmat, list(thetanew * Jnew))
              rm(Qnew, Jnew)
            }
            
          } # end if(lencidx == 1L)
          
        } # end for(k in 1:nterms)
        
        colnames(Kmat) <- Knames
        names(thetas) <- thetaNames
        return(list(K = Kmat, J = Jmat, Q = Qmat, thetas = thetas, weights = NULL))
        
      } # end if(!tprk)
      
    } # end if(is.null(thetas))
    
    # ssa 2-way interaction between cidx = c(ii,jj)
    #  * (1) add Xn[ii] * Xn[jj] to null space (if both exist)
    #       Knew = rkron(Xn[[ii]], Xn[[jj]])
    #  * (2) add Xn[ii] * Xc[jj] to contrast space (if both exist)
    #       Xnew = tcrossprod(Xn[[ii]], Qn[[ii]]) * Xc[[jj]]
    #       Qnew = tcrossprod(Qn[[ii]]) * Qc[[jj]]
    #  * (3) add Xc[ii] * Xn[jj] to contrast space (if both exist)
    #       Xnew = Xc[[ii]] * tcrossprod(Xn[[jj]], Qn[[jj]])
    #       Qnew = Qc[[ii]] * tcrossprod(Qn[[jj]])
    #  * (4) add Xc[ii] * Xc[jj] to contrast space (if both exist)
    #       Xnew = Xc[[ii]] * Xc[[jj]]
    #       Qnew = Qc[[ii]] * Qc[[jj]]
    
    # NOTE: if term ii is parametric, only do steps 1-2 
    #       because Xc[[ii]] and Qc[[ii]] don't exist
    
    # gam 2-way interaction between cidx = c(ii,jj)
    #  * (1) add Xn[ii] * Xn[jj] to null space (if both exist)
    #       Knew = rkron(Xn[[ii]], Xn[[jj]])
    #  * (2) add Xn[ii] * Xc[jj] to contrast space (if both exist)
    #       Xnew = rkron(Xn[[ii]], Xc[[jj]])
    #       Qnew = kronecker(crossprod(Qn[[ii]]), Qc[[jj]])
    #  * (3) add Xc[ii] * Xn[jj] to contrast space (if both exist)
    #       Xnew = rkron(Xc[[ii]], Xn[[jj]])
    #       Qnew = kronecker(Qc[[ii]], crossprod(Qn[[jj]]))
    #  * (4) add Xc[ii] * Xc[jj] to contrast space (if both exist)
    #       Xnew = rkron(Xc[[ii]], Xc[[jj]])
    #       Qnew = kronecker(Qc[[ii]], Qc[[jj]])
    
    
    Thetas <- thetas
    
    ### tprk (smoothing spline anova)
    if(tprk){
      
      # get nknots
      nknots <- NULL
      k <- 1L
      while(is.null(nknots)){
        nknots <- ncol(rkhs$Qc[[k]])
        k <- k + 1L
      }
      
      # initialize J and Q
      Jmat <- matrix(0, nrow = nobs, ncol = nknots)
      Qmat <- matrix(0, nrow = nknots, ncol = nknots)
      colnames(Jmat) <- rownames(Qmat) <- colnames(Qmat) <- paste("knot", 1:nknots, sep=".")
      
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
          
          # add to contrast and penalty matrices
          if(!is.null(rkhs$Xc[[cidx]])){
            Qmat <- Qmat + thetas[1] * rkhs$Qc[[cidx]]
            Jmat <- Jmat + thetas[1] * rkhs$Xc[[cidx]]
            thetas <- thetas[-1]
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
            Qmat <- Qmat + thetas[1] * Qnew
            Jmat <- Jmat + thetas[1] * rkhs$Xc[[ii]] * tcrossprod(rkhs$Xn[[jj]], rkhs$Qn[[jj]])
            thetas <- thetas[-1]
            rm(Qnew)
          }
          
          # add to contrast and penalty matrices (ii.null and jj.cont)
          if(!is.null(rkhs$Xn[[ii]]) && !is.null(rkhs$Xc[[jj]])){
            Qnew <- tcrossprod(rkhs$Qn[[ii]]) * rkhs$Qc[[jj]]
            Qmat <- Qmat + thetas[1] * Qnew
            Jmat <- Jmat + thetas[1] * tcrossprod(rkhs$Xn[[ii]], rkhs$Qn[[ii]]) * rkhs$Xc[[jj]]
            thetas <- thetas[-1]
            rm(Qnew)
          }
          
          # add to contrast and penalty matrices (ii.cont and jj.cont)
          if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xc[[jj]])){
            Qnew <- rkhs$Qc[[ii]] * rkhs$Qc[[jj]]
            Qmat <- Qmat + thetas[1] * Qnew
            Jmat <- Jmat + thetas[1] * rkhs$Xc[[ii]] * rkhs$Xc[[jj]]
            thetas <- thetas[-1]
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
            Qmat <- Qmat + thetas[1] * Qnew
            Jmat <- Jmat + thetas[1] * rkhs$Xc[[ii]] * tcrossprod(rkhs$Xn[[jj]], rkhs$Qn[[jj]]) * tcrossprod(rkhs$Xn[[kk]], rkhs$Qn[[kk]])
            thetas <- thetas[-1]
            rm(Qnew)
          }
          
          # add to contrast and penalty matrices (ii.null and jj.cont and kk.null)
          if(!is.null(rkhs$Xn[[ii]]) && !is.null(rkhs$Xc[[jj]]) && !is.null(rkhs$Xn[[kk]])){
            Qnew <- tcrossprod(rkhs$Qn[[ii]]) * rkhs$Qc[[jj]] * tcrossprod(rkhs$Qn[[kk]])
            Qmat <- Qmat + thetas[1] * Qnew
            Jmat <- Jmat + thetas[1] * tcrossprod(rkhs$Xn[[ii]], rkhs$Qn[[ii]]) * rkhs$Xc[[jj]] * tcrossprod(rkhs$Xn[[kk]], rkhs$Qn[[kk]])
            thetas <- thetas[-1]
            rm(Qnew)
          }
          
          # add to contrast and penalty matrices (ii.null and jj.null and kk.cont)
          if(!is.null(rkhs$Xn[[ii]]) && !is.null(rkhs$Xn[[jj]]) && !is.null(rkhs$Xc[[kk]])){
            Qnew <- tcrossprod(rkhs$Qn[[ii]]) * tcrossprod(rkhs$Qn[[jj]]) * rkhs$Qc[[kk]]
            Qmat <- Qmat + thetas[1] * Qnew
            Jmat <- Jmat + thetas[1] * tcrossprod(rkhs$Xn[[ii]], rkhs$Qn[[ii]]) * tcrossprod(rkhs$Xn[[jj]], rkhs$Qn[[jj]]) * rkhs$Xc[[kk]]
            thetas <- thetas[-1]
            rm(Qnew)
          }
          
          # add to contrast and penalty matrices (ii.cont and jj.cont and kk.null)
          if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xc[[jj]]) && !is.null(rkhs$Xn[[kk]])){
            Qnew <- rkhs$Qc[[ii]] * rkhs$Qc[[jj]] * tcrossprod(rkhs$Qn[[kk]])
            Qmat <- Qmat + thetas[1] * Qnew
            Jmat <- Jmat + thetas[1] * rkhs$Xc[[ii]] * rkhs$Xc[[jj]] * tcrossprod(rkhs$Xn[[kk]], rkhs$Qn[[kk]])
            thetas <- thetas[-1]
            rm(Qnew)
          }
          
          # add to contrast and penalty matrices (ii.cont and jj.null and kk.cont)
          if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xn[[jj]]) && !is.null(rkhs$Xc[[kk]])){
            Qnew <- rkhs$Qc[[ii]] * tcrossprod(rkhs$Qn[[jj]]) * rkhs$Qc[[kk]]
            Qmat <- Qmat + thetas[1] * Qnew
            Jmat <- Jmat + thetas[1] * rkhs$Xc[[ii]] * tcrossprod(rkhs$Xn[[jj]], rkhs$Qn[[jj]]) * rkhs$Xc[[kk]]
            thetas <- thetas[-1]
            rm(Qnew)
          }
          
          # add to contrast and penalty matrices (ii.null and jj.cont and kk.cont)
          if(!is.null(rkhs$Xn[[ii]]) && !is.null(rkhs$Xc[[jj]]) && !is.null(rkhs$Xc[[kk]])){
            Qnew <- tcrossprod(rkhs$Qn[[ii]]) * rkhs$Qc[[jj]] * rkhs$Qc[[kk]]
            Qmat <- Qmat + thetas[1] * Qnew
            Jmat <- Jmat + thetas[1] * tcrossprod(rkhs$Xn[[ii]], rkhs$Qn[[ii]]) * rkhs$Xc[[jj]] * rkhs$Xc[[kk]]
            thetas <- thetas[-1]
            rm(Qnew)
          }
          
          # add to contrast and penalty matrices (ii.cont and jj.cont and kk.cont)
          if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xc[[jj]]) && !is.null(rkhs$Xc[[kk]])){
            Qnew <- rkhs$Qc[[ii]] * rkhs$Qc[[jj]] * rkhs$Qc[[kk]]
            Qmat <- Qmat + thetas[1] * Qnew
            Jmat <- Jmat + thetas[1] * rkhs$Xc[[ii]] * rkhs$Xc[[jj]] * rkhs$Xc[[kk]]
            thetas <- thetas[-1]
            rm(Qnew)
          }
          
        } # end if(lencidx == 1L)
        
      } # end for(k in 1:nterms)
      
      colnames(Kmat) <- Knames
      return(list(K = Kmat, J = Jmat, Q = Qmat, thetas = Thetas))
      
    } # end if(tprk)
    
    
    ### !tprk (generalized additive model)
    if(!tprk){
      
      # initialize J and Q
      Jmat <- NULL
      Qmat <- NULL
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
          
          # add to contrast and penalty matrices
          if(!is.null(rkhs$Xc[[cidx]])){
            Qnew <- thetas[1] * rkhs$Qc[[cidx]]
            Jnew <- thetas[1] * rkhs$Xc[[cidx]]
            newnames <- paste(xnames[cidx], colnames(Jnew), sep = ".")
            colnames(Jnew) <- newnames
            if(is.null(newnames)) stop("missing names for contrast space component")
            Jnames <- c(Jnames, newnames)
            Qmat <- c(Qmat, list(Qnew))
            Jmat <- c(Jmat, list(Jnew))
            thetas <- thetas[-1]
            rm(Qnew, Jnew)
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
            Qnew <- kronecker(rkhs$Qc[[ii]], crossprod(rkhs$Qn[[jj]]))
            newnames <- charkro(paste(xnames[ii], colnames(rkhs$Xc[[ii]]), sep = "."), 
                                colnames(rkhs$Xn[[jj]]))
            if(is.null(newnames)) stop("missing names for contrast space component")
            Jnames <- c(Jnames, newnames)
            Jnew <- rowkro(rkhs$Xc[[ii]], rkhs$Xn[[jj]])
            colnames(Jnew) <- colnames(Qnew) <- newnames
            Qmat <- c(Qmat, list(thetas[1] * Qnew))
            Jmat <- c(Jmat, list(thetas[1] * Jnew))
            thetas <- thetas[-1]
            rm(Qnew, Jnew)
          }
          
          # add to contrast and penalty matrices (ii.null and jj.cont)
          if(!is.null(rkhs$Xn[[ii]]) && !is.null(rkhs$Xc[[jj]])){
            Qnew <- kronecker(crossprod(rkhs$Qn[[ii]]), rkhs$Qc[[jj]])
            newnames <- charkro(colnames(rkhs$Xn[[ii]]), 
                                paste(xnames[jj], colnames(rkhs$Xc[[jj]]), sep = "."))
            if(is.null(newnames)) stop("missing names for contrast space component")
            Jnames <- c(Jnames, newnames)
            Jnew <- rowkro(rkhs$Xn[[ii]], rkhs$Xc[[jj]])
            colnames(Jnew) <- colnames(Qnew) <- newnames
            Qmat <- c(Qmat, list(thetas[1] * Qnew))
            Jmat <- c(Jmat, list(thetas[1] * Jnew))
            thetas <- thetas[-1]
            rm(Qnew, Jnew)
          }
          
          # add to contrast and penalty matrices (ii.cont and jj.cont)
          if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xc[[jj]])){
            Qnew <- kronecker(rkhs$Qc[[ii]], rkhs$Qc[[jj]])
            newnames <- charkro(paste(xnames[ii], colnames(rkhs$Xc[[ii]]), sep = "."), 
                                paste(xnames[jj], colnames(rkhs$Xc[[jj]]), sep = "."))
            if(is.null(newnames)) stop("missing names for contrast space component")
            Jnames <- c(Jnames, newnames)
            Jnew <- rowkro(rkhs$Xc[[ii]], rkhs$Xc[[jj]])
            colnames(Jnew) <- colnames(Qnew) <- newnames
            Qmat <- c(Qmat, list(thetas[1] * Qnew))
            Jmat <- c(Jmat, list(thetas[1] * Jnew))
            thetas <- thetas[-1]
            rm(Qnew, Jnew)
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
            Qnew <- kronecker(kronecker(rkhs$Qc[[ii]], crossprod(rkhs$Qn[[jj]])), crossprod(rkhs$Qn[[kk]]))
            newnames <- charkro(charkro(paste(xnames[ii], colnames(rkhs$Xc[[ii]]), sep = "."), 
                                        colnames(rkhs$Xn[[jj]])), 
                                colnames(rkhs$Xn[[kk]]))
            if(is.null(newnames)) stop("missing names for contrast space component")
            Jnames <- c(Jnames, newnames)
            Jnew <- rowkro(rowkro(rkhs$Xc[[ii]], rkhs$Xn[[jj]]), rkhs$Xn[[kk]])
            colnames(Jnew) <- colnames(Qnew) <- newnames
            Qmat <- c(Qmat, list(thetas[1] * Qnew))
            Jmat <- c(Jmat, list(thetas[1] * Jnew))
            thetas <- thetas[-1]
            rm(Qnew, Jnew)
          }
          
          # add to contrast and penalty matrices (ii.null and jj.cont and kk.null)
          if(!is.null(rkhs$Xn[[ii]]) && !is.null(rkhs$Xc[[jj]]) && !is.null(rkhs$Xn[[kk]])){
            Qnew <- kronecker(kronecker(crossprod(rkhs$Qn[[ii]]), rkhs$Qc[[jj]]), crossprod(rkhs$Qn[[kk]]))
            newnames <- charkro(charkro(colnames(rkhs$Xn[[ii]]), 
                                        paste(xnames[jj], colnames(rkhs$Xc[[jj]]), sep = ".")), 
                                colnames(rkhs$Xn[[kk]]))
            if(is.null(newnames)) stop("missing names for contrast space component")
            Jnames <- c(Jnames, newnames)
            Jnew <- rowkro(rowkro(rkhs$Xn[[ii]], rkhs$Xc[[jj]]), rkhs$Xn[[kk]])
            colnames(Jnew) <- colnames(Qnew) <- newnames
            Qmat <- c(Qmat, list(thetas[1] * Qnew))
            Jmat <- c(Jmat, list(thetas[1] * Jnew))
            thetas <- thetas[-1]
            rm(Qnew, Jnew)
          }
          
          # add to contrast and penalty matrices (ii.null and jj.null and kk.cont)
          if(!is.null(rkhs$Xn[[ii]]) && !is.null(rkhs$Xn[[jj]]) && !is.null(rkhs$Xc[[kk]])){
            Qnew <- kronecker(kronecker(crossprod(rkhs$Qn[[ii]]), crossprod(rkhs$Qn[[jj]])), rkhs$Qc[[kk]])
            newnames <- charkro(charkro(colnames(rkhs$Xn[[ii]]), 
                                        colnames(rkhs$Xn[[jj]])), 
                                paste(xnames[kk], colnames(rkhs$Xc[[kk]]), sep = "."))
            if(is.null(newnames)) stop("missing names for contrast space component")
            Jnames <- c(Jnames, newnames)
            Jnew <- rowkro(rowkro(rkhs$Xn[[ii]], rkhs$Xn[[jj]]), rkhs$Xc[[kk]])
            colnames(Jnew) <- colnames(Qnew) <- newnames
            Qmat <- c(Qmat, list(thetas[1] * Qnew))
            Jmat <- c(Jmat, list(thetas[1] * Jnew))
            thetas <- thetas[-1]
            rm(Qnew, Jnew)
          }
          
          # add to contrast and penalty matrices (ii.cont and jj.cont and kk.null)
          if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xc[[jj]]) && !is.null(rkhs$Xn[[kk]])){
            Qnew <- kronecker(kronecker(rkhs$Qc[[ii]], rkhs$Qc[[jj]]), crossprod(rkhs$Qn[[kk]]))
            newnames <- charkro(charkro(paste(xnames[ii], colnames(rkhs$Xc[[ii]]), sep = "."), 
                                        paste(xnames[jj], colnames(rkhs$Xc[[jj]]), sep = ".")), 
                                colnames(rkhs$Xn[[kk]]))
            if(is.null(newnames)) stop("missing names for contrast space component")
            Jnames <- c(Jnames, newnames)
            Jnew <- rowkro(rowkro(rkhs$Xc[[ii]], rkhs$Xc[[jj]]), rkhs$Xn[[kk]])
            colnames(Jnew) <- colnames(Qnew) <- newnames
            Qmat <- c(Qmat, list(thetas[1] * Qnew))
            Jmat <- c(Jmat, list(thetas[1] * Jnew))
            thetas <- thetas[-1]
            rm(Qnew, Jnew)
          }
          
          # add to contrast and penalty matrices (ii.cont and jj.null and kk.cont)
          if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xn[[jj]]) && !is.null(rkhs$Xc[[kk]])){
            Qnew <- kronecker(kronecker(rkhs$Qc[[ii]], crossprod(rkhs$Qn[[jj]])), rkhs$Qc[[kk]])
            newnames <- charkro(charkro(paste(xnames[ii], colnames(rkhs$Xc[[ii]]), sep = "."), 
                                        colnames(rkhs$Xn[[jj]])), 
                                paste(xnames[kk], colnames(rkhs$Xc[[kk]]), sep = "."))
            if(is.null(newnames)) stop("missing names for contrast space component")
            Jnames <- c(Jnames, newnames)
            Jnew <- rowkro(rowkro(rkhs$Xc[[ii]], rkhs$Xn[[jj]]), rkhs$Xc[[kk]])
            colnames(Jnew) <- colnames(Qnew) <- newnames
            Qmat <- c(Qmat, list(thetas[1] * Qnew))
            Jmat <- c(Jmat, list(thetas[1] * Jnew))
            thetas <- thetas[-1]
            rm(Qnew, Jnew)
          }
          
          # add to contrast and penalty matrices (ii.null and jj.cont and kk.cont)
          if(!is.null(rkhs$Xn[[ii]]) && !is.null(rkhs$Xc[[jj]]) && !is.null(rkhs$Xc[[kk]])){
            Qnew <- kronecker(kronecker(crossprod(rkhs$Qn[[ii]]), rkhs$Qc[[jj]]), rkhs$Qc[[kk]])
            newnames <- charkro(charkro(colnames(rkhs$Xn[[ii]]), 
                                        paste(xnames[jj], colnames(rkhs$Xc[[jj]]), sep = ".")), 
                                paste(xnames[kk], colnames(rkhs$Xc[[kk]]), sep = "."))
            if(is.null(newnames)) stop("missing names for contrast space component")
            Jnames <- c(Jnames, newnames)
            Jnew <- rowkro(rowkro(rkhs$Xn[[ii]], rkhs$Xc[[jj]]), rkhs$Xc[[kk]])
            colnames(Jnew) <- colnames(Qnew) <- newnames
            Qmat <- c(Qmat, list(thetas[1] * Qnew))
            Jmat <- c(Jmat, list(thetas[1] * Jnew))
            thetas <- thetas[-1]
            rm(Qnew, Jnew)
          }
          
          # add to contrast and penalty matrices (ii.cont and jj.cont and kk.cont)
          if(!is.null(rkhs$Xc[[ii]]) && !is.null(rkhs$Xc[[jj]]) && !is.null(rkhs$Xc[[kk]])){
            Qnew <- kronecker(kronecker(rkhs$Qc[[ii]], rkhs$Qc[[jj]]), rkhs$Qc[[kk]])
            newnames <- charkro(charkro(paste(xnames[ii], colnames(rkhs$Xc[[ii]]), sep = "."), 
                                        paste(xnames[jj], colnames(rkhs$Xc[[jj]]), sep = ".")), 
                                paste(xnames[kk], colnames(rkhs$Xc[[kk]]), sep = "."))
            if(is.null(newnames)) stop("missing names for contrast space component")
            Jnames <- c(Jnames, newnames)
            Jnew <- rowkro(rowkro(rkhs$Xc[[ii]], rkhs$Xc[[jj]]), rkhs$Xc[[kk]])
            colnames(Jnew) <- colnames(Qnew) <- newnames
            Qmat <- c(Qmat, list(thetas[1] * Qnew))
            Jmat <- c(Jmat, list(thetas[1] * Jnew))
            thetas <- thetas[-1]
            rm(Qnew, Jnew)
          }
          
        } # end if(lencidx == 1L)
        
      } # end for(k in 1:nterms)
      
      colnames(Kmat) <- Knames
      return(list(K = Kmat, J = Jmat, Q = Qmat, thetas = Thetas, weights = NULL))
      
    } # end if(!tprk)
    
    
  } # end build_depe