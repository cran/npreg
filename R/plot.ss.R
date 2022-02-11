plot.ss <-
  function(x, n = 201, ci = TRUE, xseq = NULL, ...){
    if(is.null(xseq)) {
      xseq <- seq(x$fit$min, x$fit$min + x$fit$range, length.out = n)
    }
    df <- x$df
    x <- predict(object = x, x = xseq)
    if(!ci) x$se <- 0
    if(!any(names(list(...)) == "main")){
      main <- paste0("Smoothing Spline (df = ", round(df,2),")")
      plotci(x, main = main, ...)
    } else {
      plotci(x, ...)
    }
  }