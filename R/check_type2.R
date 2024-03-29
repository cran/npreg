check_type2 <- 
  function(mf, type, xdim, xrng, xlev){
    # check predictor types for new data
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: 2024-03-27
    
    # predictor types
    # par = parameteric (unpenalized)
    # nom = nominal (unordered factor)
    # ord = ordinal (ordered factor)
    # lin = linear smoothing spline
    # cub = cubic smoothing spline
    # qui = quintic smoothing spline
    # ran = random intercept (unordered factor)
    # per.lin = periodic linear smoothing spline
    # per.cub = periodic cubic smoothing spline
    # per.qui = periodic quintic smoothing spline
    # sph.2 = spherical spline (m = 2)
    # sph.3 = spherical spline (m = 3)
    # sph.4 = spherical spline (m = 4)
    # tps.lin = thin plate linear spline
    # tps.cub = thin plate cubic spline
    # tps.qui = thin plate quintic spline
    # sph = spherical spline (m = 2)
    # per = periodic (cubic) smoothing spline
    # tps = thin plate (cubic) spline
    
    # initializations
    mt <- attr(mf, "terms")                 # mt contains model info and terms 
    et <- attr(mt,"factors")                # et is effects table
    mfdim <- dim(et)                        # dim of effects table
    nobs <- dim(mf)[1]                      # total number of data points
    nxvar <- mfdim[1] - 1L                  # number of predictors
    nterm <- mfdim[2]                       # number of model terms
    xynames <- row.names(et)
    xnames <- xynames[2:(nxvar+1L)]
    
    # all possible types
    alltypes <- c("par", "nom", "ord", "lin", "cub", "qui", "per.lin", "per.cub", "per.qui", "ran",
                  "sph.2", "sph.3", "sph.4", "tps.lin", "tps.cub", "tps.qui", "sph", "per", "tps")
    ss.types <- c("lin", "cub", "qui", "per.lin", "per.cub", "per.qui", "per")
    sp.types <- c("sph.2", "sph.3", "sph.4", "sph")
    
    # check types
    for(k in 1:nxvar){
      
      # get k+1
      kp1 <- k + 1L
      
      # get class of new data
      xc <- class(mf[,k+1])[1]
      if(!any(xc == c("character", "factor", "ordered", "integer", "numeric", "matrix"))) stop(paste("Input",xynames[kp1],"needs to be of class factor, ordered, integer, numeric, or matrix."))
      
      # compatibility between variable and rktype
      if(type[[k]] == "par"){
        
        # parametric effect
        if(xc == "matrix"){
          if(ncol(mf[,kp1]) > 1) stop(paste("Input 'type' for",xynames[kp1],"is 'par' but",xynames[kp1],"is a matrix with more than 1 column.\n  Parametric effects should be of class 'factor', 'integer', or 'numeric'."))
          mf[,kp1] <- as.vector(mf[,kp1])
          xc <- class(mf[,kp1])[1]
        }
        xattr <- attr(xlev[[k]], "ordered")
        if(is.null(xattr)){
          # numeric effect
          if(!any(xc == c("integer", "numeric"))) stop(paste("Input 'type' for",xynames[kp1],"is a parametric numeric effect ('par')\n but newdata for",xynames[kp1],"is not 'integer' or 'numeric'."))
          newrng <- range(mf[,kp1])
          if(newrng[1] < xrng[[k]][1] | newrng[2] > xrng[[k]][2]) warning(paste("Input 'newdata' for",xynames[kp1],"is out of the range of the observed data."))
        } else if(xattr) {
          # ordered factor
          if(xc != "ordered") stop(paste("Input 'type' for",xynames[kp1],"is a parameteric ordered factor ('par')\n but newdata for",xynames[kp1],"is not an 'ordered' factor."))
          newlev <- levels(mf[,kp1])
          newidx <- match(newlev, xlev[[k]])
          if(any(is.na(newidx))) stop(paste("Input 'newdata' for",xynames[kp1],"contains factor levels not in the observed data."))
        } else {
          # unordered factor
          if(xc != "factor") stop(paste("Input 'type' for",xynames[kp1],"is a parameteric unordered factor ('par') but newdata for",xynames[kp1],"is not an unordered 'factor'."))
          newlev <- levels(mf[,kp1])
          newidx <- match(newlev, xlev[[k]])
          if(any(is.na(newidx))) stop(paste("Input 'newdata' for",xynames[kp1],"contains factor levels not in the observed data."))
        }
        
      } else if(type[[k]] == "ran"){
        
        # random intercept
        if(xc != "factor") {
          warning(paste("Input 'type' for",xynames[kp1],"is 'ran' but",xynames[kp1],"is not a factor.\n  Using as.factor() to coerce",xynames[kp1],"into a factor."))
          mf[,kp1] <- as.factor(mf[,kp1])
        }
        newlev <- levels(mf[,kp1])
        newidx <- match(newlev, xlev[[k]])
        if(any(is.na(newidx))) stop(paste("Input 'newdata' for",xynames[kp1],"contains factor levels not in the observed data."))
        
      } else if(type[[k]] == "nom"){
        
        # nominal spline
        if(xc != "factor") {
          warning(paste("Input 'type' for",xynames[kp1],"is 'nom' but",xynames[kp1],"is not a factor.\n  Using as.factor() to coerce",xynames[kp1],"into a factor."))
          mf[,kp1] <- as.factor(mf[,kp1])
        }
        newlev <- levels(mf[,kp1])
        newidx <- match(newlev, xlev[[k]])
        if(any(is.na(newidx))) stop(paste("Input 'newdata' for",xynames[kp1],"contains factor levels not in the observed data."))
        
      } else if(type[[k]] == "ord"){
        
        # ordinal spline 
        if(xc != "ordered"){
          warning(paste("Input 'type' for",xynames[kp1],"is 'ord' but",xynames[kp1],"is not an ordered factor.\n  Using as.ordered() to coerce",xynames[kp1],"into an ordered factor."))
          mf[,kp1] <- as.ordered(mf[,kp1])
        }
        newlev <- levels(mf[,kp1])
        newidx <- match(newlev, xlev[[k]])
        if(any(is.na(newidx))) stop(paste("Input 'newdata' for",xynames[kp1],"contains factor levels not in the observed data."))
        
      } else if(any(type[[k]] == ss.types)){
        
        # polynomial smoothing spline
        if(xc == "matrix"){
          if(ncol(mf[,kp1]) > 1) stop(paste("Input 'type' for",xynames[kp1],"is '",type[[k]],"' but",xynames[kp1],"is a matrix with more than 1 column.\n  Use a thin-plate spline ('tps') for multidimensional predictors."))
          mf[,kp1] <- as.vector(mf[,kp1])
          xc <- class(mf[,kp1])[1]
        }
        if(!any(xc == c("integer", "numeric"))) stop(paste("Input 'type' for",xynames[kp1],"is '",type[[k]],"' but",xynames[kp1],"is not of class 'integer' or 'numeric'.\n  Polynomial smoothing splines require numeric predictors."))
        newrng <- range(mf[,kp1])
        if(newrng[1] < xrng[[k]][1] | newrng[2] > xrng[[k]][2]) warning(paste("Input 'newdata' for",xynames[kp1],"is out of the range of the observed data."))
        
      } else if(any(type[[k]] == sp.types)){
        
        # spherical spline
        if(xc != "matrix") stop(paste("Input 'type' for",xynames[kp1],"is '",type[[k]],"' but",xynames[kp1],"is not of class 'matrix'.\n  Spherical splines require a 2-dimensional predictor (i.e., n x 2 matrix) containing lat/long."))
        if(ncol(mf[,kp1]) != 2L) stop(paste("Input 'type' for",xynames[kp1],"is '",type[[k]],"' but",xynames[kp1],"is not a matrix with 2 columns.\n  Spherical splines require a 2-dimensional predictor (i.e., n x 2 matrix) containing lat/long."))
        for(d in 1:2){
          newrng <- range(mf[,kp1][,d])
          if(newrng[1] < xrng[[k]][1,d] | newrng[2] > xrng[[k]][2,d]) warning(paste("Input 'newdata' for",xynames[kp1],"is out of the range of the observed data."))
        }
        
      } else {
        
        # thin-plate spline
        if(!any(class((mf[,kp1])[1])[1] == c("integer", "numeric"))) stop(paste("Input",xynames[kp1],"should be a matrix with numeric (or integer) values."))
        mf[,kp1] <- as.matrix(mf[,kp1])
        newdim <- ncol(mf[,kp1])
        if(newdim != xdim[k]) stop(paste("Input 'newdata' for",xynames[kp1],"has a different dimension", newdim ,"than original data",xdim[k],"."))
        for(d in 1:newdim){
          newrng <- range(mf[,kp1][,d])
          if(newrng[1] < xrng[[k]][1,d] | newrng[2] > xrng[[k]][2,d]) warning(paste("Input 'newdata' for",xynames[kp1],"is out of the range of the observed data."))
        }
        
      } # end if(type[[k]] == "par")
    } # for(k in 1:nxvar)
    
    return(NULL)
    
  }