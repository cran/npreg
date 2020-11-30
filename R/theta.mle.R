theta.mle <-
  function(y, mu, theta, wt = 1, 
           maxit = 100, tol = .Machine$double.eps^0.5){
    # MLE of theta for NegBinomial
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2020-08-23
    
    iter <- 0
    delta <- tol + 1
    if(missing(theta)) theta <- 1 / mean(wt * (y / mu - 1)^2)
    while(iter < maxit & abs(delta) > tol){
      g <- theta.grad(theta, y, mu, wt)
      i <- theta.info(theta, y, mu, wt)
      delta <- g / i
      theta <- abs(theta + delta)
      iter <- iter + 1
    }
    attr(theta, "SE") <- ifelse(i > 0, 1 / sqrt(i), NA)
    attr(theta, "iter") <- iter
    theta
    
  } # end theta.mle

theta.grad <-
  function(theta, y, mu, wt){
    sum(wt * (1 + log(theta) - digamma(theta) - (theta + y)/(theta + mu) - 
                log(theta + mu) + digamma(y + theta)) )
  }

theta.info <-
  function(theta, y, mu, wt){
    sum(wt * (-1/theta + trigamma(theta) + (mu - y)/(theta + mu)^2 +
                1/(theta + mu) - trigamma(y + theta)) )
  }