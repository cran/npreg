plotci <- 
  function(x, y, se, level = 0.95, crit.val = NULL, 
           add = FALSE, col.ci = "gray80", alpha = 0.5, 
           bars = FALSE, bw = 0.2, linkinv = NULL, ...){
    
    # check x and y inputs
    if(missing(y)){
      if(is.list(x) | is.data.frame(x)){
        y <- as.numeric(x$y)
        xx <- as.numeric(x$x)
        n <- length(xx)
        if(length(y) != n) stop("If 'x' is a list, the 'x' and 'y' elements must have length n.")
        if(!is.null(x$se)) {
          se <- x$se
          if(length(se) != n) stop("If 'x' is a list, the 'se' element must have length n.")
        }
        x <- xx
      } else if(is.matrix(x)){
        y <- as.numeric(x[,2])
        if(ncol(x) > 2L) se <- as.numeric(x[,3])
        x <- as.numeric(x[,1])
        n <- length(x)
      } else {
        n <- length(x)
        y <- as.numeric(x)
        x <- 1:n
      }
      if(any(is.na(y)) | any(is.infinite(y)) | any(is.nan(y))) 
        stop("'x' and 'y' cannot contain missing (NA), infinite (Inf), or undefined (NaN) values")
    } else {
      x <- as.numeric(x)
      y <- as.numeric(y)
      n <- length(x)
      if(length(y) != n) stop("'x' and 'y' must have the same length")
    }
    
    # check se input
    missing.se <- TRUE
    if(!missing(se)){
      missing.se <- FALSE
      se <- as.numeric(se)
      if(length(se) != n) stop("Input 'se' must have the same length as 'x' and 'y'.")
      if(any(se < 0)) stop("Input 'se' must contain non-negative standard error values.")
    } 
    
    # check level & crit.val inputs
    if(is.null(crit.val)){
      level <- as.numeric(level)
      if(level <= 0 | level >=1 ) stop("Input 'level' must satisfy:  0 < level < 1")
      crit.val <- qnorm(1 - (1 - level)/2)
    } else {
      crit.val <- as.numeric(crit.val[1])
      if(crit.val <= 0) stop("Input 'crit.val' must be a positive scalar.")
    }
    
    # check linkinv input
    link <- FALSE
    if(!is.null(linkinv)){
      link <- TRUE
      if(!is.function(linkinv)) stop("Input 'linkinv' must be a function.")
      ylink <- linkinv(y)
      if(length(ylink) != n) stop("Input 'linkinv' function produced 'y' of incorrect length.\nThe function results must satisfy:  length(y) = length(linkinv(y))")
    } 
    
    # sort data
    sx <- sort(x, index = TRUE)
    y <- y[sx$ix]
    if(link) ylink <- ylink[sx$ix]
    if(!missing.se) se <- se[sx$ix]
    x <- x[sx$ix]
    
    # get arguments
    args <- list(...)
    
    # add x and y
    args$x <- x
    args$y <- if(link) ylink else y
    
    # missing type?
    if(is.null(args$type)) args$type <- "l"
    
    # missing lwd?
    if(is.null(args$lwd)) args$lwd <- 2
    
    # missing labels?
    if(is.null(args$xlab)) args$xlab <- "x"
    if(is.null(args$ylab)) args$ylab <- "y"
    
    # missing y limits?
    if(is.null(args$ylim) && !missing.se) {
      fitrng <- range(c(y - crit.val * se, y + crit.val * se))
      rngdif <- fitrng[2] - fitrng[1]
      args$ylim <- c(fitrng[1] - 0.05 * rngdif, fitrng[2] + 0.05 * rngdif)
    }
    
    # plot x vs y
    if(add){
      do.call(lines.default, args)
    } else {
      do.call(plot.default, args)
    }
    
    # add se
    if(!missing.se){
      col <- rgb(t(col2rgb(col.ci) / 255), alpha = alpha)
      if(!bars){
        xpoly <- c(x, rev(x))
        ypoly <- c(y - crit.val * se, rev(y + crit.val * se))
        if(link) ypoly <- linkinv(ypoly)
        polygon(x = xpoly, y = ypoly, col = col, border = col)
        do.call(lines, args)
      } else {
        for(i in 1:n){
          iargs <- args
          iargs$col <- col
          iargs$x <- rep(iargs$x[i], 2)
          iargs$y <- c(args$y[i] - crit.val * se[i],
                       args$y[i] + crit.val * se[i])
          if(link) iargs$y <- linkinv(iargs$y)
          do.call(lines, iargs)
          if(bw > 0){
            yval <- iargs$y
            xrng <- range(x)
            xrng <- xrng[2] - xrng[1]
            sw <- bw * xrng / n
            iargs$x <- c(iargs$x[1] - sw, iargs$x[1] + sw)
            iargs$y <- rep(yval[1], 2)
            do.call(lines, iargs)
            iargs$y <- rep(yval[2], 2)
            do.call(lines, iargs)
          }
        }
        do.call(lines, args)
      }
    }
    
  } # end 