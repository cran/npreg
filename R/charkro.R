charkro <- 
  function(x, y, sep = ":"){
    
    if(is.null(x) | is.null (y)) return(NULL)
    x <- as.character(x)
    y <- as.character(y)
    nx <- length(x)
    ny <- length(y)
    z <- vector("character", nx * ny)
    for(i in 1:nx){
      zind <- (1 + (i - 1) * ny):(ny + (i - 1) * ny)
      z[zind] <- paste(x[i], y, sep = sep)
    }
    z
    
  }