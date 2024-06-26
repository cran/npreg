check_type <- 
  function(mf, type = NULL, xrange = NULL){
    # check and/or guess predictor type
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: 2024-03-27
    
    # predictor types
    # par = parametric (unpenalized)
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
    
    # get types
    if(is.null(type)) {
      type <- vector("list", nxvar)
    } else {
      if(nxvar == 1L){
        type <- list(type[[1]])
        names(type) <- xnames
        if(!any(type[[1]] == alltypes)) stop("Input 'type' is not correctly specified.")
      } else {
        if(!is.list(type)) stop("Input 'type' must be a list with named elements.")
      }
    }
    xrng <- xlev <- vector("list", nxvar)
    rktype <- rep(NA, nxvar)
    names(rktype) <- names(xrng) <- names(xlev) <- xynames[-1]
    nm <- match(xynames[-1], names(type))
    for(k in 1:length(nm)){
      kp1 <- k + 1L
      if(is.na(nm[k])){
        # best guess of rktype
        xc <- class(mf[,kp1])[1]
        if((xc == "factor") | (xc == "character")){
          nlev <- nlevels(mf[,kp1])
          rktype[k] <- ifelse(nlev < 20, "nom", "ran")
          xrng[[k]] <- c(1L, nlev)
          xlev[[k]] <- levels(mf[,kp1])
        } else if(xc == "ordered"){
          rktype[k] <- "ord"
          xrng[[k]] <- c(1L, nlevels(mf[,kp1]))
          xlev[[k]] <- levels(mf[,kp1])
        } else if((xc == "integer") | (xc == "numeric")){
          rktype[k] <- "cub"
          xrng[[k]] <- range(mf[,kp1])
        } else if(xc == "matrix"){
          if(!any(class((mf[,kp1])[1]) == c("integer","numeric"))){
            stop(paste("Input",xynames[kp1],"needs to be a matrix of integers or numerics."))
          }
          ncx <- ncol(mf[,kp1])
          if(ncx > 5L) stop(paste("Input",xynames[kp1],"has too many columns (need 5 or less)."))
          rktype[k] <- ifelse(ncx >= 4L, "tps.qui", "tps.cub")
          xrng[[k]] <- apply(mf[,kp1], 2, range)
        } else {
          stop(paste("Input",xynames[kp1],"needs to be of class factor, ordered, integer, numeric, or matrix."))
        } # end if(xc=="factor")
      } else {
        # evaluate input rktype
        if(!any(type[[nm[k]]] == alltypes)){
          stop(paste("Input of 'type' for",xynames[kp1],"is not correctly specified."))
        }
        rktype[k] <- type[[nm[k]]]
        # compatibility between variable and rktype
        xc <- class(mf[,kp1])[1]
        if(!any(xc == c("character", "factor", "ordered", "integer", "numeric", "matrix"))) stop(paste("Input",xynames[kp1],"needs to be of class factor, ordered, integer, numeric, or matrix."))
        if(rktype[k] == "par"){
          # parametric effect
          if(xc == "matrix"){
            if(ncol(mf[,kp1]) > 1) stop(paste("Input 'type' for",xynames[kp1],"is 'par' but",xynames[kp1],"is a matrix with more than 1 column.\n  Parametric effects should be of class 'factor', 'integer', or 'numeric'."))
            mf[,kp1] <- as.vector(mf[,kp1])
            xc <- class(mf[,kp1])[1]
          }
          if(any(xc == c("character", "factor", "ordered"))){
            if(xc == "character") mf[,kp1] <- as.factor(mf[,kp1])
            xrng[[k]] <- c(1L, nlevels(mf[,kp1]))
            xlev[[k]] <- levels(mf[,kp1])
            attr(xlev[[k]], "ordered") <- ifelse(is.ordered(mf[,kp1]), TRUE, FALSE)
          } else {
            xrng[[k]] <- range(mf[,kp1])
          }
        } else if(rktype[k] == "ran"){
          # random intercept
          if(xc != "factor") {
            warning(paste("Input 'type' for",xynames[kp1],"is 'ran' but",xynames[kp1],"is not a factor.\n  Using as.factor() to coerce",xynames[kp1],"into a factor."))
            mf[,kp1] <- as.factor(mf[,kp1])
          }
          xrng[[k]] <- c(1L, nlevels(mf[,kp1]))
          xlev[[k]] <- levels(mf[,kp1])
        } else if(rktype[k] == "nom"){
          # nominal spline requested
          if(xc != "factor") {
            warning(paste("Input 'type' for",xynames[kp1],"is 'nom' but",xynames[kp1],"is not a factor.\n  Using as.factor() to coerce",xynames[kp1],"into a factor."))
            mf[,kp1] <- as.factor(mf[,kp1])
          }
          xrng[[k]] <- c(1L, nlevels(mf[,kp1]))
          xlev[[k]] <- levels(mf[,kp1])
        } else if(rktype[k] == "ord"){
          # ordinal spline requested
          if(xc != "ordered"){
            warning(paste("Input 'type' for",xynames[kp1],"is 'ord' but",xynames[kp1],"is not an ordered factor.\n  Using as.ordered() to coerce",xynames[kp1],"into an ordered factor."))
            mf[,kp1] <- as.ordered(mf[,kp1])
          }
          xrng[[k]] <- c(1L, nlevels(mf[,kp1]))
          xlev[[k]] <- levels(mf[,kp1])
        } else if(any(rktype[k] == ss.types)){
          # polynomial smoothing spline requested
          if(xc == "matrix"){
            if(ncol(mf[,kp1]) > 1) stop(paste("Input 'type' for",xynames[kp1],"is '",rktype[k],"' but",xynames[kp1],"is a matrix with more than 1 column.\n  Use a thin-plate spline ('tps') for multidimensional predictors."))
            mf[,kp1] <- as.vector(mf[,kp1])
            xc <- class(mf[,kp1])[1]
          }
          if(!any(xc == c("integer", "numeric"))) stop(paste("Input 'type' for",xynames[kp1],"is '",rktype[k],"' but",xynames[kp1],"is not of class 'integer' or 'numeric'.\n  Polynomial smoothing splines require numeric predictors."))
          xrng[[k]] <- range(mf[,kp1])
        } else if(any(rktype[k] == sp.types)){
          # spherical spline requested
          if(xc != "matrix") stop(paste("Input 'type' for",xynames[kp1],"is '",rktype[k],"' but",xynames[kp1],"is not of class 'matrix'.\n  Spherical splines require a 2-dimensional predictor (i.e., n x 2 matrix) containing lat/long."))
          if(ncol(mf[,kp1]) != 2L) stop(paste("Input 'type' for",xynames[kp1],"is '",rktype[k],"' but",xynames[kp1],"is not a matrix with 2 columns.\n  Spherical splines require a 2-dimensional predictor (i.e., n x 2 matrix) containing lat/long."))
          xrng[[k]] <- apply(mf[,kp1], 2, range)
        } else {
          # thin-plate spline requested
          if(xc == "matrix" &&  ncol(mf[,kp1]) > 5L) stop(paste("Input",xynames[kp1],"has too many columns (need 5 or less)."))
          if(!any(class((mf[,kp1])[1])[1] == c("integer", "numeric"))) stop(paste("Input",xynames[kp1],"should be a matrix with numeric (or integer) values."))
          if(xc == "matrix"){
            if(ncol(mf[,kp1]) > 1L && rktype[k] == "tps.lin") stop(paste("Input 'type' for",xynames[kp1],"must be 'tps.cub' or 'tps.qui' because",xynames[kp1],"is multidimensional."))
            if(ncol(mf[,kp1]) > 3L && rktype[k] == "tps.cub") stop(paste("Input 'type' for",xynames[kp1],"must be 'tps.qui' because",xynames[kp1],"has d > 3 dimensions."))
            xrng[[k]] <- apply(mf[,kp1], 2, range)
          } else {
            xrng[[k]] <- range(mf[,kp1])
          }
        } # end if(rktype[k] == "par")
      } # end if(is.na(nm[k]))
    } # for(k in 1:length(nm))
    
    # given xrange?
    if(!is.null(xrange)){
      xrange <- as.list(xrange)
      xrange.names <- names(xrange)
      if(is.null(xrange.names)) stop("Input 'xrange' must be a named list")
      for(k in 1:nxvar){
        kid <- match(xnames[k], xrange.names)
        if(!is.na(kid)){
          obs.range <- xrng[[k]]
          if(obs.range[1] < xrange[[kid]][1]){
            warning("\nxrange$", xrange.names[kid], "[1] = ", xrange[[kid]][1], " > min(data$", xrange.names[kid],") = ", obs.range[1], "\n Redefining to minimum observed data value")
            xrange[[kid]][1] <- obs.range[1]
          }
          if(obs.range[2] > xrange[[kid]][2]){
            warning("\nxrange$", xrange.names[kid], "[2] = ", xrange[[kid]][2], " < max(data$", xrange.names[kid],") = ", obs.range[2], "\n Redefining to maximum observed data value")
            xrange[[kid]][2] <- obs.range[2]
          }
          xrng[[k]] <- xrange[[kid]][1:2]
        } # end if(!is.na(kid))
      } # end for(k in 1:nxvar)
    } # end if(!is.null(xrange))
    
    # return results
    return(list(type = rktype, xrng = xrng, xlev = xlev))
    
  }