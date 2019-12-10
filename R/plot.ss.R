plot.ss <-
  function(x, n = 201, ...){
    xmin <- x$fit$min
    xseq <- seq(xmin, xmin + x$fit$range, length.out = n)
    df <- x$df
    x <- predict(object = x, x = xseq)
    if(!any(names(list(...)) == "main")){
      main <- paste0("Smoothing Spline (df = ", round(df,2),")")
      plotci(x, main = main, ...)
    } else {
      plotci(x, ...)
    }
  }