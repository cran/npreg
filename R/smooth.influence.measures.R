smooth.influence.measures <- 
  function(model, infl = smooth.influence(model)){
    # regression diagnostic measures for smooth models
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2022-05-03
    
    # is influential?
    is.influential <- function(infmat, n, df) {
      k <- ncol(infmat) - 4L
      if (n <= k) stop("too few cases i with s_ii > 0), n < k")
      absmat <- abs(infmat)
      r <- cbind(absmat[, 1L:k] > 1,
                 absmat[, k + 1] > 3 * sqrt(df/(n - df)),
                 abs(1 - infmat[, k + 2]) > (3 * df)/(n - df),
                 pf(infmat[, k + 3], df, n - df) > 0.5,
                 infmat[, k + 4] > (3 * df)/n)
      attributes(r) <- attributes(infmat)
      r
    }
    
    # make influence matrix
    infmat <- dfbetas(model, infl = infl)
    colnames(infmat) <- paste0("dfb.", c("1_", abbreviate(names(coef(model)[-1]))))
    infmat <- cbind(infmat, 
                    dffit = dffits(model, infl = infl),
                    cov.r = cov.ratio(model, infl = infl),
                    cook.d = cooks.distance(model),
                    hat = infl$hat)
    
    # which cases are influential?
    is.inf <- is.influential(infmat, sum(infl$hat > 0), model$df)
    
    # return results
    ans <- list(infmat = infmat, is.inf = is.inf, call = model$call)
    class(ans) <- "infl"
    ans
    
  } # end smooth.influence.measures