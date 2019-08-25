penalty_poly <-
  function(x, m = 2, xmin = min(x), xmax = max(x), 
           periodic = FALSE, rescale = FALSE){
    # Polynomial Smoothing Spline Penalty
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Update: 2019-04-04
    
    # initializations
    m <- as.integer(m)
    if(m < 1L | m > 3L) stop("Input 'm' must be 1 (linear), 2 (cubic), or 3 (quintic).")
    x <- (as.numeric(x) - xmin) / (xmax - xmin)
    
    # get kernel function name
    kerns <- c("lin", "cub", "qui")
    fname <- paste0(kerns[m], "kern", 0)
    
    # evaluate kernel function
    Q <- outer(X = x, Y = x, FUN = fname, periodic = periodic)
    
    # rescale?
    if(rescale) Q <- Q / mean(diag(Q))
    
    # return result
    return(Q)
    
  }