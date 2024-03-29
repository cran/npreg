predict.ss <-
  function(object, x, deriv = 0, se.fit = TRUE, ...){
    # predict method for class "ss"
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2020-10-22
    
    X <- basis.poly(x = x, knots = object$fit$knot, 
                    m = object$fit$m, d = deriv,
                    xmin = object$fit$min,
                    xmax = object$fit$min + object$fit$range,
                    periodic = object$fit$periodic,
                    intercept = TRUE, bernoulli = object$fit$bernoulli)
    pred <- data.frame(x = x, y = X %*% object$fit$coef)
    if(se.fit){
      pred$se <- sqrt(rowSums((X %*% object$fit$cov.sqrt)^2))
    }
    pred
    
  } # end predict.ss