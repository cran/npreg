build_rkhs <-
  function(x, type, knots, xrng){
    # build reproducing kernel Hilbert space components
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: 2019-10-30
    
    ### initializations
    nxvar <- length(x)
    xnames <- names(x)
    
    ### kernel types
    alltypes <- c("par", "nom", "ord", "lin", "cub", "qui", "per.lin", "per.cub", "per.qui",
                  "sph.lin", "sph.cub", "sph.qui", "tps.lin", "tps.cub", "tps.qui", "sph", "per", "tps")
    ss.types <- c("lin", "cub", "qui", "per.lin", "per.cub", "per.qui", "per")
    sp.types <- c("sph.lin", "sph.cub", "sph.qui", "sph")
    
    ### initialize components
    Xnull <- Xcont <- Qnull <- Qcont <- vector("list", nxvar)
    names(Xnull) <- names(Xcont) <- names(Qnull) <- names(Qcont) <- xnames
    
    
    ### sweep through variables
    for(j in 1:nxvar){
      
      tj <- type[j]
      
      if(tj == "par"){
        
        # parametric effect
        if(any(class(x[[j]]) == c("factor", "ordered"))){
          jform <- as.formula(paste("~", xnames[j]))
          Xnull[[j]] <- model.matrix(jform, data = x)[,-1,drop=FALSE]
          Qnull[[j]] <- model.matrix(jform, data = knots)[,-1,drop=FALSE]  
        } else {
          Xnull[[j]] <- as.matrix(x[[j]])
          Qnull[[j]] <- as.matrix(knots[[j]])
        }
        
      } else if(tj == "nom"){
        
        # nominal smoothing spline
        Xcont[[j]] <- basis_nom(x = x[[j]], knots = knots[[j]], K = xrng[[j]][2])
        Qcont[[j]] <- penalty_nom(x = knots[[j]], K = xrng[[j]][2])
        
      } else if(tj == "ord"){
        
        # ordinal smoothing spline
        Xcont[[j]] <- basis_ord(x = x[[j]], knots = knots[[j]], K = xrng[[j]][2])
        Qcont[[j]] <- penalty_ord(x = knots[[j]], K = xrng[[j]][2])
        
      } else if(any(tj == ss.types)){
        
        # polynomial smoothing spline
        m <- ifelse(tj == "lin" | tj == "per.lin", 1L, 
                    ifelse(tj == "cub" | tj == "per.cub" | tj == "per", 2L, 3L))
        periodic <- ifelse(any(tj == c("per.lin", "per.cub", "per.qui", "per")), TRUE, FALSE)
        if(periodic | tj == "lin"){
          Xcont[[j]] <- basis_poly(x = x[[j]], knots = knots[[j]], m = m, xmin = xrng[[j]][1], xmax = xrng[[j]][2], periodic = periodic)
          Qcont[[j]] <- penalty_poly(x = knots[[j]], m = m, xmin = xrng[[j]][1], xmax = xrng[[j]][2], periodic = periodic)
        } else {
          X <- basis_poly(x = x[[j]], knots = knots[[j]], m = m, xmin = xrng[[j]][1], xmax = xrng[[j]][2], periodic = periodic)
          Q <- basis_poly(x = knots[[j]], knots = knots[[j]], m = m, xmin = xrng[[j]][1], xmax = xrng[[j]][2], periodic = periodic)
          nulldim <- ncol(Q) - nrow(Q)
          nullnames <- xnames[j]
          if(m == 3L) nullnames <- c(nullnames, paste0(xnames[j], "^2"))
          Xnull[[j]] <- X[,1:nulldim,drop=FALSE]
          Xcont[[j]] <- X[,(nulldim+1):ncol(X)]
          Qnull[[j]] <- Q[,1:nulldim,drop=FALSE]
          Qcont[[j]] <- Q[,(nulldim+1):ncol(Q)]
          colnames(Xnull[[j]]) <- nullnames
          rm(X, Q)
        }
        
      } else if(any(tj == sp.types)){
        
        # spherical spline
        m <- ifelse(tj == "sph.lin", 1L, ifelse(tj == "sph.cub" | tj == "sph", 2L, 3L))
        Xcont[[j]] <- basis_sph(x = x[[j]], knots = knots[[j]], m = m)
        Qcont[[j]] <- penalty_sph(x = knots[[j]], m = m)
        
      } else {
        
        # thin-plate spline
        xdim <- ncol(as.matrix(x[[j]]))
        m <- ifelse(tj == "tps.lin", 1L, ifelse(tj == "tps.cub" | tj == "tps", 2L, 3L))
        if(tj == "tps.lin"){
          Xcont[[j]] <- basis_tps(x = x[[j]], knots = knots[[j]], m = m)
          Qcont[[j]] <- penalty_tps(x = knots[[j]], m = m)
        } else {
          X <- basis_tps(x = x[[j]], knots = knots[[j]], m = m)
          Q <- basis_tps(x = knots[[j]], knots = knots[[j]], m = m)
          nulldim <- ncol(Q) - nrow(Q)
          Xnull[[j]] <- X[,1:nulldim,drop=FALSE]
          Xcont[[j]] <- X[,(nulldim+1):ncol(X)]
          Qnull[[j]] <- Q[,1:nulldim,drop=FALSE]
          Qcont[[j]] <- Q[,(nulldim+1):ncol(Q)]
          rm(X, Q)
        }
        
      } # end if(tj == "par")
      
    } # end for(j in 1:nxvar)
    
    return(list(Xn = Xnull, Xc = Xcont, Qn = Qnull, Qc = Qcont))
    
  }