plot.ss <-
  function(x, se = TRUE, ...){
    if(se) x$se <- x$sigma * sqrt(x$lev)
    if(!any(names(list(...)) == "main")){
      main <- paste0("Smoothing Spline (df = ", round(x$df,2),")")
      plotci(x, main = main, ...)
    } else {
      plotci(x, ...)
    }
  }