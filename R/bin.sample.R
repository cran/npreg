bin.sample <- 
  function(x, nbin = 5, size = 1, equidistant = FALSE, 
           index.return = FALSE, breaks.return = FALSE){
    # bin sample vector, matrix, or data frame
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: 2021-04-11
    
    # check 'x'
    xc <- class(x)[1]
    if(any(xc == c("integer", "numeric"))){
      x <- matrix(x)
      xc <- "matrix"
    } else if(any(xc == c("factor", "ordered"))){
      x <- data.frame(x = x)
      xc <- "data.frame"
    }
    if(!any(xc == c("matrix", "data.frame"))) stop("Input 'x' must be an object of class 'matrix' or 'data.frame'.")
    nobs <- nrow(x)
    nvar <- ncol(x)
    
    # check number of bins
    nbin <- as.integer(nbin)
    if(length(nbin) != nvar) nbin <- rep(nbin[1], nvar)
    if(any(nbin < 1L)) stop("Input 'nbin' must contain non-negative integers.")
    
    # check size
    size <- as.integer(size[1])
    if(size < 1L) stop("Input 'size' must be a positive integer.")
    
    # check equidistant
    equidistant <- as.logical(equidistant)
    if(length(equidistant) != nvar) equidistant <- rep(equidistant[1], nvar)
    
    # check breaks.return and index.return
    breaks.return <- as.logical(breaks.return[1])
    index.return <- as.logical(index.return[1])
    
    # convert data.frame
    factorID <- rep(FALSE, nvar)
    if(xc == "data.frame"){
      xclass <- sapply(x, function(x) class(x)[1])
      if(any(xclass == "factor") | any(xclass == "ordered")){
        facid <- c(which(xclass == "factor"), which(xclass == "ordered"))
        factorID[facid] <- TRUE
        for(j in facid){
          nbin[j] <- nlevels(x[,j])
          x[,j] <- as.integer(x[,j])
        }
      }
      x <- as.matrix(x)
    }
    
    # initializations (if returning breaks)
    if(breaks.return) {
      binmat <- matrix(NA, nobs, nvar)
      binmid <- vector("list", nvar)
    }
    
    # multidimensional binning
    binvec <- rep(1L, nobs)
    alpha <- 1L
    for(j in 1:nvar){
      if(factorID[j]){
        newbrk <- seq(0.5, nbin[j] + 0.5, length.out = (nbin[j] + 1))
      } else if(equidistant[j]){
        newbrk <- seq(min(x[,j]), max(x[,j]), length.out = (nbin[j] + 1L))
      } else {
        newbrk <- quantile(x[,j], probs = seq(0, 1, length.out = (nbin[j] + 1L)))
      }
      newbin <- .bincode(x[,j], breaks = newbrk, include.lowest = TRUE)
      binvec <- binvec + alpha * (newbin - 1L)
      alpha <- alpha * nbin[j]
      if(breaks.return){
        binmat[,j] <- newbin
        factor.correct <- ifelse(factorID[j], 0.5, 0)
        binmid[[j]] <- newbrk[1:nbin[j]] + factor.correct
      }
    } # end for(j in 1:nvar)
    
    # bin-sample
    binsam <- function(x, nsize) {
      lenx <- length(x)
      if(lenx == 1L){
        return(x)
      } else {
        return(sample(x, min(nsize, lenx)))
      }
    }
    
    # return breaks
    if(breaks.return){
      ix <- unname(unlist(tapply(1:nobs, binvec, binsam, nsize = 1)))
      bx <- matrix(0, nrow = length(ix), ncol = nvar)
      for(j in 1:nvar) bx[,j] <- binmid[[j]][binmat[ix,j]]
      if(nvar == 1L) bx <- as.vector(bx)
    } 
    
    # return samples
    ix <- as.vector(unname(unlist(tapply(1:nobs, binvec, binsam, nsize = size))))
    if(breaks.return | index.return){
      if(breaks.return & index.return){
        return(list(x = x[ix,], ix = ix, bx = bx))
      } else if(breaks.return) {
        return(list(x = x[ix,], bx = bx))
      } else {
        return(list(x = x[ix,], ix = ix))
      }
    } else {
      return(x[ix,])
    }
    
  }