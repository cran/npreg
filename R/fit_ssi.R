fit_ssi <- 
  function(x, y = NULL, w = NULL, df, spar = NULL, lambda = NULL, 
           method = c("GCV", "OCV", "GACV", "ACV", "REML", "ML", "AIC", "BIC"),
           m = 2L, periodic = FALSE, all.knots = FALSE, nknots = .nknots.smspl, 
           knots = NULL, keep.data = TRUE, df.offset = 0, penalty = 1, 
           control.spar = list(), tol = 1e-6 * IQR(x), bernoulli = TRUE,
           xmin = NULL, xmax = NULL, homosced = FALSE, iter.max = 1L){
    
    # fit smoothing spline model to data
    m0 <- ss(x = x, y = y, w = w, df = df, spar = spar, lambda = lambda,
             method = method, m = m, periodic = periodic, all.knots = all.knots,
             nknots = .nknots.smspl, knots = knots, keep.data = keep.data,
             df.offset = df.offset, penalty = penalty, control.spar = control.spar,
             tol = tol, bernoulli = bernoulli, xmin = xmin, xmax = xmax)
    
    # get standardized residuals
    res <- rstandard(m0)
    
    # fit smoothing spline model to absolute residuals
    m1 <- ss(x = x, y = abs(res), method = method, m = m, periodic = periodic, 
             all.knots = all.knots, nknots = .nknots.smspl, knots = knots, 
             keep.data = keep.data, df.offset = df.offset, penalty = penalty, 
             control.spar = control.spar, tol = tol, bernoulli = bernoulli, 
             xmin = xmin, xmax = xmax)
    
    # updates weights
    w <- weights(m0) / predict(m1, x = x)$y^2
    
    # fit smoothing spline model to data with updated weights
    m0 <- ss(x = x, y = y, w = w, df = df, spar = spar, lambda = lambda,
             method = method, m = m, periodic = periodic, all.knots = all.knots,
             nknots = .nknots.smspl, knots = knots, keep.data = keep.data,
             df.offset = df.offset, penalty = penalty, control.spar = control.spar,
             tol = tol, bernoulli = bernoulli, xmin = xmin, xmax = xmax)
    
    # iterative updates?
    if(iter.max == 1L){
      m0 <- c(m0, list(homosced = homosced, iter.max = iter.max, iter = 1L))
      class(m0) <- "ss"
      return(m0)
    }
    
    # set convergence tolerance
    cnvg.tol <- m0$fit$control.spar$tol
    
    # initializations
    iter <- 1L
    delta <- cnvg.tol + 1
    yhat.old <- predict(m0, x = x)$y
    while(iter < iter.max && delta > cnvg.tol){
      
      # get standardized residuals
      res <- rstandard(m0)
      
      # fit smoothing spline model to absolute residuals
      m1 <- ss(x = x, y = abs(res), method = method, m = m, periodic = periodic, 
               all.knots = all.knots, nknots = .nknots.smspl, knots = knots, 
               keep.data = keep.data, df.offset = df.offset, penalty = penalty, 
               control.spar = control.spar, tol = tol, bernoulli = bernoulli, 
               xmin = xmin, xmax = xmax)
      
      # updates weights
      w <- weights(m0) / predict(m1, x = x)$y^2
      
      # fit smoothing spline model to data with updated weights
      m0 <- ss(x = x, y = y, w = w, df = df, spar = spar, lambda = lambda,
               method = method, m = m, periodic = periodic, all.knots = all.knots,
               nknots = .nknots.smspl, knots = knots, keep.data = keep.data,
               df.offset = df.offset, penalty = penalty, control.spar = control.spar,
               tol = tol, bernoulli = bernoulli, xmin = xmin, xmax = xmax)
      
      # check for convergence
      yhat <- predict(m0, x = x)$y
      delta <- mean((yhat.old - yhat)^2) / (mean(yhat.old^2) + cnvg.tol)
      iter <- iter + 1L
      yhat.old <- yhat
      
    } # end while(iter < iter.max && delta > cnvg.tol)
    
    # return results
    m0 <- c(m0, list(homosced = homosced, iter.max = iter.max, iter = iter, delta = delta))
    class(m0) <- "ss"
    return(m0)
    
  }