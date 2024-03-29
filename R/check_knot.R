check_knot <-
  function(mf, type, xrng, xlev, tprk, knots = NULL){
    # check and/or sample spline knots
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: 2023-04-06
    
    ### initializations
    mt <- attr(mf, "terms")                 # mt contains model info and terms 
    et <- attr(mt,"factors")                # et is effects table
    mfdim <- dim(et)                        # dim of effects table
    nobs <- dim(mf)[1]                      # total number of data points
    nxvar <- mfdim[1] - 1L                  # number of predictors
    nterm <- mfdim[2]                       # number of model terms
    xynames <- row.names(et)
    xnames <- xynames[2:(nxvar+1L)]
    
    
    ### initialize types
    ss.types <- c("lin", "cub", "qui", "per.lin", "per.cub", "per.qui", "per")
    sp.types <- c("sph.2", "sph.3", "sph.4", "sph")
    tp.types <- c("tps.lin", "tps.cub", "tps.qui", "tps")
    
    
    ### no knots provided
    if(is.null(knots)){
      knots <- vector("list", nxvar)
      names(knots) <- xnames
      if(!tprk){
        nk <- ifelse(nxvar == 1L, 10, 5)
        for(j in 1:nxvar) knots[[j]] <- knot1samp(mf[,j+1], nk)
      } else {
        if(nxvar == 1L){
          knots[[1]] <- knot1samp(mf[,2], 10)
        } else {
          if(nxvar < 14L){
            knotid <- bin.sample(mf[,-1], index.return = TRUE)$ix
            if(length(knotid) > 100) knotid <- sample(knotid, size = 100)
          } else {
            knotid <- sample(1:nobs, size = min(100, nobs))
          }
          for(j in 1:nxvar) knots[[j]] <- mf[knotid,j+1]
        }
      }
      return(knots)
    } else {
      if(is.data.frame(knots)) knots <- as.list(knots)
    } # end if(is.null(knots))
    
    
    ### knots for tprk (smoothing spline anova)
    if(tprk){
      # check class of knots input
      if(!is.list(knots) && !is.vector(knots)) {
        stop("When 'tprk' is TRUE, the input 'knots' must be either: 
             1) a scalar giving the total number of knots to sample,
             2) a vector of integers indexing which rows of data are the knots, or
             3) a list with named elements giving the knot values for each predictor**.
             ** Requires the same number of knots for each predictor.")
      }
      
      # input a list of knots
      if(is.list(knots)){
        
        # check if input marginal knots
        dfun <- function(x) ifelse(is.matrix(x), nrow(x), length(x))
        nknots <- sapply(knots, dfun)
        if(min(nknots) != max(nknots)){
          knots <- as.list(expand.grid(knots))
        }
        
        # make named knots list
        given.knots <- knots
        knots <- vector("list", nxvar)
        names(knots) <- xnames
        nknots <- ifelse(is.matrix(given.knots[[1]]), 
                         nrow(given.knots[[1]]), 
                         length(given.knots[[1]]))
        
        # sweep through all variables
        for(j in 1:nxvar){
          
          # check if knots provided for j-th term
          kmj <- match(xnames[j], names(given.knots))
          if(is.na(kmj)){
            # no knots provided; sample some
            knots[[j]] <- knot1samp(mf[,j+1], n = nknots)
            ncheck <- ifelse(is.matrix(knots[[j]]), nrow(knots[[j]]), length(knots[[j]]))
            if(ncheck < nknots){
              if(is.matrix(knots[[j]])){
                knots[[j]] <- rbind(knots[[j]], knots[[j]][sample.int(ncheck, size = nknots - ncheck, replace = TRUE),])
              } else {
                knots[[j]] <- c(knots[[j]], knots[[j]][sample.int(ncheck, size = nknots - ncheck, replace = TRUE)])
              }
            }
            
          } else {
            
            # provided a vector/matrix (giving knots values)
            knotsj <- given.knots[[kmj]]
            ncheck <- ifelse(is.matrix(knotsj), nrow(knotsj), length(knotsj))
            if(ncheck != nknots) stop("When 'tprk' is TRUE and a list of knots values is provided,\n the same number of knots is required for each predictor.")
            kclass <- class(knotsj)[1]
            if(type[j] == "ran"){
              if(kclass != "factor"){
                knotsj <- factor(knotsj, levels = xlev[[j]])
                if(any(is.na(knotsj))) stop(paste("Input knots for",xnames[j],"contain factor levels not present in the data."))
              }
              if(!identical(xlev[[j]], levels(knotsj))) stop(paste("Input 'knots' for",xnames[j],"do not match the levels of the corresponding variable."))
            } else if(type[j] == "nom"){
              if(kclass != "factor"){
                #warning(paste0("Input 'type' for ",xnames[j]," is 'nom' but input knots for ",xnames[j]," are not a factor.\n  Using factor() to coerce knots$",xnames[j]," into a factor."))
                knotsj <- factor(knotsj, levels = xlev[[j]])
                if(any(is.na(knotsj))) stop(paste("Input knots for",xnames[j],"contain factor levels not present in the data."))
              }
              if(!identical(xlev[[j]], levels(knotsj))) stop(paste("Input 'knots' for",xnames[j],"do not match the levels of the corresponding variable."))
            } else if(type[j] == "ord"){
              if(kclass != "ordered"){
                #warning(paste0("Input 'type' for ",xnames[j]," is 'ord' but input knots for ",xnames[j]," are not an ordered factor.\n  Using ordered() to coerce knots$",xnames[j]," into an ordered factor."))
                knotsj <- factor(knotsj, levels = xlev[[j]], ordered = TRUE)
                if(any(is.na(knotsj))) stop(paste("Input knots for",xnames[j],"contain factor levels not present in the data."))
              }
              if(!identical(xlev[[j]], levels(knotsj))) stop(paste("Input 'knots' for",xnames[j],"do not match the levels of the corresponding variable."))
            } else if(any(type[j] == ss.types)){
              if(kclass == "matrix"){
                if(ncol(knotsj) > 1) stop(paste("Input 'type' for",xnames[j],"is '",type[j],"' but input knots for",xnames[j],"is a matrix with more than 1 column."))
                knotsj <- as.vector(knotsj)
                kclass <- class(knotsj)[1]
              }
              if(!any(kclass == c("integer", "numeric"))) stop(paste("Input 'type' for",xnames[j],"is '",type[j],"' but input knots for",xnames[j],"are not of class 'integer' or 'numeric'.\n  Polynomial smoothing splines require numeric predictors and knots."))
              if(min(knotsj) < xrng[[j]][1] | max(knotsj) > xrng[[j]][2]) warning(paste("Input 'knots' for",xnames[j],"are outside of the range of the data."))
            } else if(any(type[j] == sp.types)){
              if(kclass != "matrix") stop(paste("Input 'type' for",xnames[j],"is '",type[j],"' but input knots for",xnames[j],"are not of class 'matrix'.\n  Spherical splines require a 2-dimensional predictor and knots (i.e., n x 2 matrix) containing lat/long."))
              if(ncol(knotsj) != 2L) stop(paste("Input 'type' for",xnames[j],"is '",type[j],"' but input knots for",xnames[j],"are not a matrix with 2 columns.\n  Spherical splines require a 2-dimensional predictor and knots (i.e., n x 2 matrix) containing lat/long."))
              for(k in 1:2) if(min(knotsj[,k]) < xrng[[j]][1,k] | max(knotsj[,k]) > xrng[[j]][2,k]) warning(paste("Input 'knots' for",xnames[j],"are outside of the range of the data for dimension",k))
            } else if(any(type[j] == tp.types)){
              if(class(mf[,j+1])[1] == "matrix"){
                if(ncol(mf[,j+1]) != ncol(knotsj)) stop(paste("Input knots 'type' for",xnames[j],"do not have the same number of columns as input variable."))
                for(k in 1:ncol(knotsj)) if(min(knotsj[,k]) < xrng[[j]][1,k] | max(knotsj[,k]) > xrng[[j]][2,k]) warning(paste("Input 'knots' for",xnames[j],"are outside of the range of the data for dimension",k))
              } else {
                if(min(knotsj) < xrng[[j]][1] | max(knotsj) > xrng[[j]][2]) warning(paste("Input 'knots' for",xnames[j],"are outside of the range of the data."))
              }
            }
            knots[[j]] <- knotsj
            
          } # end if(is.na(kmj))
          
        } # end for(j in 1:nxvar)
        
        return(knots)
        
      } # end if(is.list(knots))
      
      # input a scalar (giving # of knots)
      if(length(knots) == 1L){
        
        # make named knots list
        nknots <- as.integer(knots)
        knots <- vector("list", nxvar)
        names(knots) <- xnames
        
        # feasibility checks
        if(nknots < 3L) stop("Need to use at least three knots.")
        if(nknots > nobs) stop("Input too many knots! Need number of knots less than n.")
        
        # using all data points as knots?
        if(nknots == nobs){
          knotid <- 1:nobs
          for(j in 1:nxvar) knots[[j]] <- mf[knotid,j+1]
          return(knots)
        }
        
        # single predictor variable
        if(nxvar == 1L){
          knots[[1]] <- knot1samp(mf[,2L], nknots)
          return(knots)
        }
        
        # multiple predictor variables
        if(nxvar < 14L){
          knotid <- bin.sample(mf[,-1,drop=FALSE], index.return = TRUE)$ix
          lk <- length(knotid)
          if(lk > nknots){
            knotid <- sample(knotid, nknots)
          } else if(lk < nknots) {
            knotid <- c(knotid, sample(seq(1,nobs)[-knotid], nknots - lk))
          }
        } else {
          knotid <- sample(1:nobs, size = nknots)
        }
        for(j in 1:nxvar) knots[[j]] <- mf[knotid,j+1]
        return(knots)
        
      } # end if(length(knots) == 1L)
      
      # input a vector (giving row indices of knots)
      knotid <- unique(as.integer(knots))
      if(any(knotid < 1) | any(knotid > nobs)) stop("Input 'knots' is a vector, but values are not in the range 1,...,n.")
      if(length(knotid) < 3L) stop("Need to use at least three knots.")
      knots <- vector("list", nxvar)
      names(knots) <- xnames
      for(j in 1:nxvar) knots[[j]] <- mf[knotid,j+1]
      return(knots)
      
    } # end if(tprk)
    
    
    ### knots for !tprk (generalized additive model)
    if(!tprk){
      
      # check class of knots input
      if(!is.vector(knots)) {
        stop("When 'tprk' is FALSE, the input 'knots' must be either: 
             1) a scalar giving the common number of knots for each predictor,
             2) a list with named elements giving the number of knots for each predictor, or
             3) a list with named elements giving the knot values for each predictor.")
      }
      
      # input a list of knots
      if(is.list(knots)){
        
        # make named knots list
        given.knots <- knots
        knots <- vector("list", nxvar)
        names(knots) <- xnames
        
        # sweep through all variables
        for(j in 1:nxvar){
          
          # check if knots provided for j-th term
          kmj <- match(xnames[j], names(given.knots))
          if(is.na(kmj)){
            # no knots provided; sample some
            knots[[j]] <- knot1samp(mf[,j+1])
            
          } else {
            
            knotsj <- given.knots[[kmj]]
            lenknot <- length(knotsj)
            if(lenknot == 1L){
              # provided a scalar (giving # of knots)
              knotsj <- as.integer(knotsj)
              if(knotsj < 1L) stop(paste("Input too few knots for predictor",xnames[j]))
              knots[[j]] <- knot1samp(mf[,j+1], knotsj)
              
            } else {
              # provided a vector/matrix (giving knots values)
              kclass <- class(knotsj)[1]
              if(type[j] == "ran"){
                if(kclass != "factor"){
                  knotsj <- factor(knotsj, levels = xlev[[j]])
                  if(any(is.na(knotsj))) stop(paste("Input knots for",xnames[j],"contain factor levels not present in the data."))
                }
                if(!identical(xlev[[j]], levels(knotsj))) stop(paste("Input 'knots' for",xnames[j],"do not match the levels of the corresponding variable."))
              } else if(type[j] == "nom"){
                if(kclass != "factor"){
                  #warning(paste0("Input 'type' for ",xnames[j]," is 'nom' but input knots for ",xnames[j]," are not a factor.\n  Using factor() to coerce knots$",xnames[j]," into a factor."))
                  knotsj <- factor(knotsj, levels = xlev[[j]])
                  if(any(is.na(knotsj))) stop(paste("Input knots for",xnames[j],"contain factor levels not present in the data."))
                }
                if(!identical(xlev[[j]], levels(knotsj))) stop(paste("Input 'knots' for",xnames[j],"do not match the levels of the corresponding variable."))
              } else if(type[j] == "ord"){
                if(kclass != "ordered"){
                  #warning(paste0("Input 'type' for ",xnames[j]," is 'ord' but input knots for ",xnames[j]," are not an ordered factor.\n  Using ordered() to coerce knots$",xnames[j]," into an ordered factor."))
                  knotsj <- factor(knotsj, levels = xlev[[j]], ordered = TRUE)
                  if(any(is.na(knotsj))) stop(paste("Input knots for",xnames[j],"contain factor levels not present in the data."))
                }
                if(!identical(xlev[[j]], levels(knotsj))) stop(paste("Input 'knots' for",xnames[j],"do not match the levels of the corresponding variable."))
              } else if(any(type[j] == ss.types)){
                if(kclass == "matrix"){
                  if(ncol(knotsj) > 1) stop(paste("Input 'type' for",xnames[j],"is '",type[j],"' but input knots for",xnames[j],"is a matrix with more than 1 column."))
                  knotsj <- as.vector(knotsj)
                  kclass <- class(knotsj)[1]
                }
                if(!any(kclass == c("integer", "numeric"))) stop(paste("Input 'type' for",xnames[j],"is '",type[j],"' but input knots for",xnames[j],"are not of class 'integer' or 'numeric'.\n  Polynomial smoothing splines require numeric predictors and knots."))
                knotsj <- sort(unique(knotsj))
                if(min(knotsj) < xrng[[j]][1] | max(knotsj) > xrng[[j]][2]) warning(paste("Input 'knots' for",xnames[j],"are outside of the range of the data."))
              } else if(any(type[j] == sp.types)){
                if(kclass != "matrix") stop(paste("Input 'type' for",xnames[j],"is '",type[j],"' but input knots for",xnames[j],"are not of class 'matrix'.\n  Spherical splines require a 2-dimensional predictor and knots (i.e., n x 2 matrix) containing lat/long."))
                if(ncol(knotsj) != 2L) stop(paste("Input 'type' for",xnames[j],"is '",type[j],"' but input knots for",xnames[j],"are not a matrix with 2 columns.\n  Spherical splines require a 2-dimensional predictor and knots (i.e., n x 2 matrix) containing lat/long."))
                for(k in 1:2) if(min(knotsj[,k]) < xrng[[j]][1,k] | max(knotsj[,k]) > xrng[[j]][2,k]) warning(paste("Input 'knots' for",xnames[j],"are outside of the range of the data for dimension",k))
              } else if(any(type[j] == tp.types)){
                if(class(mf[,j+1])[1] == "matrix"){
                  if(ncol(mf[,j+1]) != ncol(knotsj)) stop(paste("Input knots 'type' for",xnames[j],"do not have the same number of columns as input variable."))
                  for(k in 1:ncol(knotsj)) if(min(knotsj[,k]) < xrng[[j]][1,k] | max(knotsj[,k]) > xrng[[j]][2,k]) warning(paste("Input 'knots' for",xnames[j],"are outside of the range of the data for dimension",k))
                } else {
                  if(min(knotsj) < xrng[[j]][1] | max(knotsj) > xrng[[j]][2]) warning(paste("Input 'knots' for",xnames[j],"are outside of the range of the data."))
                }
              }
              knots[[j]] <- knotsj
              
            } # end if(lenknot == 1L)
            
          } # end if(is.na(kmj))
          
        } # end for(j in 1:nxvar)
        
        return(knots)
        
      } # end if(is.list(knots))
      
      # input a scalar (giving # of knots)
      if(length(knots) == 1L){
        
        # make named knots list
        nknots <- as.integer(knots)
        knots <- vector("list", nxvar)
        names(knots) <- xnames
        
        # feasibility checks
        if(nknots < 3L) stop("Need to use at least three knots.")
        if(nknots > nobs) stop("Input too many knots! Need number of knots less than n.")
        
        # using all data points as knots?
        if(nknots == nobs){
          knotid <- 1:nobs
          for(j in 1:nxvar) knots[[j]] <- mf[knotid,j+1]
          return(knots)
        }
        
        # select nknots for each predictor
        for(j in 1:nxvar) knots[[j]] <- knot1samp(mf[,j+1], nknots)
        return(knots)
        
      } else {
        stop("When 'tprk' is FALSE, the input 'knots' must be either: 
             1) a scalar giving the common number of knots for each predictor,
             2) a list with named elements giving the number of knots for each predictor, or
             3) a list with named elements giving the knot values for each predictor.")
      } # end if(length(knots) == 1L)
      
    } # end if(!tprk)
    
  }