plot.gsm <-
  function(x, terms = x$terms, se = TRUE, n = 201, intercept = FALSE,
           ask = prod(par("mfcol")) < length(terms) && dev.interactive(),
           zero.line = TRUE, zero.lty = 3, zero.col = "black", ncolor = 21, 
           colors = NULL, rev = FALSE, zlim = NULL, lty.col = NULL, 
           legend.xy = "top", main = NULL, xlab = NULL, ylab = NULL, ...){
    # Plot Effects of Generalized Smooth Model Fits
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Update: 2024-03-28
    
    if(!inherits(x, "gsm")) stop("Input 'x' must be of class 'gsm'.")
    class(x) <- "sm"
    plot(x, terms = terms, se = se, n = n, intercept = intercept, ask = ask, 
         zero.line = zero.line, zero.lty = zero.lty, zero.col = zero.col,
         ncolor = ncolor, colors = colors, rev = rev, zlim = zlim, 
         lty.col = lty.col, legend.xy = legend.xy, main = main,
         xlab = xlab, ylab = ylab, ...)
     
  }