check_control <-
  function(epsilon = 1e-08, maxit = 25, trace = FALSE, 
           lower = 0, upper = 1, tol = 1e-08){
    if (!is.numeric(epsilon) || epsilon <= 0) 
      stop("value of 'control$epsilon' must be > 0")
    if (!is.numeric(maxit) || maxit <= 0) 
      stop("value of 'control$maxit' must be > 0")
    if (!is.logical(trace))
      stop("value of 'control$trace' must be logical (T/F)")
    if (!is.numeric(lower) || lower < 0)
      stop("value of 'control$lower' must be positive")
    if (!is.numeric(upper) || upper <= lower)
      stop("value of 'control$upper' must be > value of 'control$lower'")
    if (!is.numeric(tol) || tol <= 0) 
      stop("value of 'control$tol' must be > 0")
    list(epsilon = epsilon, maxit = maxit, trace = trace,
         lower = lower, upper = upper, tol = tol)
  }