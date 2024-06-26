pred_rkhs <-
  function(x, type, knots, xrng, xlev){
    # build reproducing kernel Hilbert space components
    # for predicting (i.e., w/o penalty components)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2023-04-06
    
    ### initializations
    nxvar <- length(x)
    xnames <- names(x)
    
    ### kernel types'
    alltypes <- c("par", "nom", "ord", "lin", "cub", "qui", "per.lin", "per.cub", "per.qui", "ran",
                  "sph.2", "sph.3", "sph.4", "tps.lin", "tps.cub", "tps.qui", "sph", "per", "tps")
    ss.types <- c("lin", "cub", "qui", "per.lin", "per.cub", "per.qui", "per")
    sp.types <- c("sph.2", "sph.3", "sph.4", "sph")
    
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
          colnames(Xnull[[j]]) <- xnames[j]
        }
        
      } else if(tj == "ran"){
        
        # nominal smoothing spline
        Xcont[[j]] <- outer(X = factor(x[[j]], levels = xlev[[j]]), Y = knots[[j]], FUN = "==") + 0.0
        colnames(Xcont[[j]]) <- as.character(knots[[j]])
        
      } else if(tj == "nom"){
        
        # nominal smoothing spline
        Xcont[[j]] <- basis.nom(x = factor(x[[j]], levels = xlev[[j]]), 
                                knots = knots[[j]], K = xrng[[j]][2])
        
      } else if(tj == "ord"){
        
        # ordinal smoothing spline
        Xcont[[j]] <- basis.ord(x = factor(x[[j]], levels = xlev[[j]], ordered = TRUE), 
                                knots = knots[[j]], K = xrng[[j]][2])
        
      } else if(any(tj == ss.types)){
        
        # polynomial smoothing spline
        m <- ifelse(tj == "lin" | tj == "per.lin", 1L, 
                    ifelse(tj == "cub" | tj == "per.cub" | tj == "per", 2L, 3L))
        periodic <- ifelse(any(tj == c("per.lin", "per.cub", "per.qui", "per")), TRUE, FALSE)
        if(periodic | tj == "lin"){
          Xcont[[j]] <- basis.poly(x = x[[j]], knots = knots[[j]], m = m, xmin = xrng[[j]][1], xmax = xrng[[j]][2], periodic = periodic)
        } else {
          X <- basis.poly(x = x[[j]], knots = knots[[j]], m = m, xmin = xrng[[j]][1], xmax = xrng[[j]][2], periodic = periodic)
          Q <- basis.poly(x = knots[[j]], knots = knots[[j]][1], m = m, xmin = xrng[[j]][1], xmax = xrng[[j]][2], periodic = periodic)
          nulldim <- m - 1L
          nullnames <- xnames[j]
          if(m == 3L) nullnames <- c(nullnames, paste0(xnames[j], "^2"))
          Xnull[[j]] <- X[,1:nulldim,drop=FALSE]
          Xcont[[j]] <- X[,(nulldim+1):ncol(X)]
          Qnull[[j]] <- Q[,1:nulldim,drop=FALSE]
          colnames(Xnull[[j]]) <- nullnames
          rm(X, Q)
        }
        
      } else if(any(tj == sp.types)){
        
        # spherical spline
        m <- ifelse(tj == "sph.2" | tj == "sph", 2L, ifelse(tj == "sph.3", 3L, 4L))
        Xcont[[j]] <- basis.sph(x = x[[j]], knots = knots[[j]], m = m)
        
      } else {
        
        # thin-plate spline
        xdim <- ncol(as.matrix(x[[j]]))
        m <- ifelse(tj == "tps.lin", 1L, ifelse(tj == "tps.cub" | tj == "tps", 2L, 3L))
        if(tj == "tps.lin"){
          Xcont[[j]] <- basis.tps(x = x[[j]], knots = knots[[j]], m = m)
        } else {
          X <- basis.tps(x = x[[j]], knots = knots[[j]], m = m)
          Q <- basis.tps(x = knots[[j]], knots = knots[[j]], m = m)
          nulldim <- ncol(Q) - nrow(Q)
          Xnull[[j]] <- X[,1:nulldim,drop=FALSE]
          Xcont[[j]] <- X[,(nulldim+1):ncol(X)]
          Qnull[[j]] <- Q[,1:nulldim,drop=FALSE]
          Qcont[[j]] <- Q[,(nulldim+1):ncol(Q)]
        }
        
      } # end if(tj == "par")
      
    } # end for(j in 1:nxvar)
    
    return(list(Xn = Xnull, Xc = Xcont, Qn = Qnull))
    
  }