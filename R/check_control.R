check_control <-
  function(lower = 0, upper = 1, tol = 1e-08, iterlim = 5000L,
           epsilon = 1e-08, maxit = 25, epsilon.out = 1e-6, 
           maxit.out = 10, print.level = 0L){
    if (!is.numeric(lower))
      stop("value of 'control$lower' must be a numeric")
    if (!is.numeric(upper) || upper <= lower)
      stop("value of 'control$upper' must be > value of 'control$lower'")
    if (!is.numeric(tol) || tol <= 0) 
      stop("value of 'control$tol' must be > 0")
    if (!is.numeric(iterlim) || iterlim <= 0) 
      stop("value of 'control$iterlim' must be > 0")
    if (!is.numeric(epsilon) || epsilon <= 0) 
      stop("value of 'control$epsilon' must be > 0")
    if (!is.numeric(maxit) || maxit <= 0) 
      stop("value of 'control$maxit' must be > 0")
    if (!is.numeric(epsilon.out) || epsilon.out <= 0) 
      stop("value of 'control$epsilon' must be > 0")
    if (!is.numeric(maxit.out) || maxit.out <= 0) 
      stop("value of 'control$maxit' must be > 0")
    print.level <- as.integer(print.level)
    if (!any(print.level == c(0L, 1L, 2L)))
      stop("Input 'control$print.level' must be either:\n0 (no printing), 1 (minimal printing), or 2 (full printing).")
    list(lower = lower, upper = upper, tol = tol, iterlim = iterlim, 
         print.level = print.level, epsilon = epsilon, maxit = maxit, 
         epsilon.out = epsilon.out, maxit.out = maxit.out)
  }