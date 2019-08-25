fit_gsm <-
  function(y, depe, lambda = NULL, tprk = TRUE, 
           control = check_control(), method = "GCV",
           family = check_family(gaussian)){
    # fit semiparametric model
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2019-08-20
    
    
    #########***#########   INITIALIZATIONS   #########***#########
    
    # info
    nobs <- nrow(depe$K)
    nsdim <- ncol(depe$K)
    if(!tprk){
      Nknots <- sapply(depe$J, ncol)
      depe$J <- do.call(cbind, depe$J)
    }
    nknots <- ncol(depe$J)
    nullindx <- 1:nsdim
    
    # reparameterize contrast space
    eps <- .Machine$double.eps
    if(tprk){
      
      Qeig <- eigen(depe$Q, symmetric = TRUE)
      Qrnk <- sum(Qeig$values > nknots * eps * Qeig$values[1])
      if(Qrnk == 1L){
        Qprj <- matrix(Qeig$vectors[,1] / sqrt(Qeig$values[1]), nrow = nknots, ncol = 1)
      } else{
        Qprj <- Qeig$vectors[,1:Qrnk] %*% diag(1 / sqrt(Qeig$values[1:Qrnk]))
      }
      Rmat <- depe$J %*% Qprj
      
    } else {
      
      cknots <- c(0, cumsum(Nknots))
      Rmat <- Qprj <- vector("list", length(depe$Q))
      for(k in 1:length(depe$Q)){
        indx <- seq(cknots[k] + 1, cknots[k+1])
        Qeig <- eigen(depe$Q[[k]], symmetric = TRUE)
        Qrnk <- sum(Qeig$values > nknots * eps * Qeig$values[1])
        if(Qrnk == 1L){
          Qprj[[k]] <- matrix(Qeig$vectors[,1] / sqrt(Qeig$values[1]), nrow = nknots, ncol = 1)
        } else{
          Qprj[[k]] <- Qeig$vectors[,1:Qrnk] %*% diag(1 / sqrt(Qeig$values[1:Qrnk]))
        }
        Rmat[[k]] <- depe$J[,indx] %*% Qprj[[k]]
      }
      Rmat <- do.call("cbind", Rmat)
      Qrnk <- ncol(Rmat)
      
    } # end if(tprk)
    
    # reverse transformation
    Tmat <- matrix(0, nsdim + nknots, nsdim + Qrnk)
    Tmat[nullindx,nullindx] <- diag(nsdim)
    #Tmat[nullindx,-nullindx] <- (-1) * solve(crossprod(X.w[,nullindx])) %*% crossprod(X.w[,nullindx], R.w)
    if(tprk){
      Tmat[-nullindx,-nullindx] <- Qprj
    } else {
      row.offset <- col.offset <- nsdim
      for(k in 1:length(Qprj)){
        nrowk <- nrow(Qprj[[k]])
        ncolk <- ncol(Qprj[[k]])
        Tmat[row.offset + 1:nrowk, col.offset + 1:ncolk] <- Qprj[[k]]
        row.offset <- row.offset + nrowk
        col.offset <- col.offset + ncolk
      }
    }
    #Tmat <- Tmat %*% Xsvd$v %*% diag(1 / Xsvd$d)
    
    #########***#########   ESTIMATE COEFS   #########***#########
    
    # get initial beta0
    beta0 <- family$linkfun(mean(y))
    
    # estimate coefficients (and possibly tune lambda)
    if(is.null(lambda)){
      
      # optimize lambda
      interval <- c(control$lower, control$upper)
      opt <- optimize(f = tune.gsm, interval = interval,
                      y = y, Kmat = depe$K, Rmat = Rmat, 
                      weights = depe$weights, beta0 = beta0,
                      tprk = tprk, control = control, 
                      family = family, method = method, 
                      tol = control$tol)
      
      # collect results
      spar <- opt$minimum
      lambda <- 256^(3*(spar-1))
      crit <- as.numeric(opt$objective)
      df <- attr(opt$objective, "df")
      iter <- attr(opt$objective, "iter")
      
      # estimate coefficients
      opt <-tune.gsm(spar, y = y, Kmat = depe$K, Rmat = Rmat, 
                     weights = depe$weights, beta0 = beta0,
                     tprk = tprk, control = control, 
                     family = family, method = method, 
                     return.coef = TRUE)
      beta <- opt$coef
      
    } else {
      
      # fit with given lambda
      spar <- 1 + log(lambda, base = 256) / 3
      opt <- tune.gsm(spar, y = y, Kmat = depe$K, Rmat = Rmat, 
                      weights = depe$weights, beta0 = beta0,
                      tprk = tprk, control = control, 
                      family = family, method = method, 
                      return.coef = TRUE)
      
      # collect results
      crit <- as.numeric(opt$objective)
      df <- attr(opt$objective, "df")
      iter <- attr(opt$objective, "iter")
      beta <- opt$coef
      
    } # end if(is.null(lambda))
    
    # evaluate log-likelihood
    eta <- cbind(depe$K, Rmat) %*% beta
    mu <- family$linkinv(eta)
    devres <- family$dev.resids(y, mu, depe$weights)
    dev <- sum(devres)
    famLL <- sum( family$logLik(y, nobs, mu, depe$weights, dev) )
    
    # null deviance
    wtdmu <- sum(depe$weights * y) / sum(depe$weights)
    nulldev <- sum(family$dev.resids(y, wtdmu, depe$weights))
    
    # get standard errors of linear predictors
    se.lp <- sqrt(rowSums((cbind(depe$K, Rmat) %*% opt$cov.sqrt)^2))
    
    # retransform coefficients and cov.sqrt
    beta <- Tmat %*% beta
    opt$cov.sqrt <- Tmat %*% opt$cov.sqrt
    
    # name coefficients
    names(beta) <- c(colnames(depe$K), colnames(depe$J))
    
    # smoothness penalty
    if(tprk){
      penalty <- as.numeric(crossprod(beta[-nullindx], depe$Q %*% beta[-nullindx]))
    } else {
      penalty <- 0
      scoefs <- beta[-nullindx]
      for(k in 1:length(depe$Q)){
        indx <- seq(cknots[k] + 1, cknots[k+1])
        penalty <- penalty + as.numeric(crossprod(scoefs[indx], depe$Q[[k]] %*% scoefs[indx]))
      }
    }
    
    # R-squared
    r.squared <- as.numeric(cor(y, mu)^2)
    
    # collect results
    res <- list(iter = iter,
                residuals = opt$resid,
                null.deviance = nulldev,
                linear.predictors = as.numeric(eta), 
                se.lp = se.lp,
                deviance = dev,
                cv.crit = crit,
                df = df, 
                nsdf = nsdim,
                r.squared = r.squared,
                dispersion = opt$disp,
                logLik = famLL, 
                aic = 2*(df - famLL),
                bic = log(nobs)*df - 2*famLL,
                spar = 1 + log(lambda, base = 256)/3,
                lambda = lambda, 
                penalty = penalty,
                coefficients = beta,
                cov.sqrt = opt$cov.sqrt)
    return(res)
    
  } # end fit_gsm

tune.gsm <-
  function(spar, y, Kmat, Rmat, weights, beta0,
           tprk = TRUE, control = check_control(),
           family = check_family(gaussian), 
           method = "GCV", return.coef = FALSE){
    
    # initialize parameters for IRLS
    lambda <- 256^(3*(spar-1))
    nobs <- length(y)
    nlam <- nobs * lambda
    if(family$family != "gaussian"){
      mustart <- NULL
      eval(family$initialize)
      mu <- mustart
    } else {
      mu <- rep(0, nobs)
    }
    eta <- family$linkfun(mu)
    Xmat <- cbind(Kmat, Rmat)
    Pmat <- diag(rep(c(0, nlam), times = c(ncol(Kmat), ncol(Rmat))))
    
    # initialize log-likelihood
    devres <- family$dev.resids(y, mu, weights)
    dev <- sum(devres)
    suppressWarnings({famLL.old <- sum( family$logLik(y, nobs, mu, weights, dev) )})
    if(is.nan(famLL.old) | is.na(famLL.old) | is.infinite(famLL.old)) famLL.old <- (-1)*.Machine$double.xmax
    
    # irls update...
    iter <- 0L
    vtol <- control$epsilon + 1
    oldcoef <- rep(c(beta0, 0), times = c(1, ncol(Xmat) - 1L))
    while(iter < control$maxit & vtol > control$epsilon){
      
      # update design matrices
      mueta <- family$mu.eta(eta)
      varx <- family$variance(mu)
      vsqrt <- as.numeric(sqrt(mueta^2/varx))
      Kw <- as.matrix(Kmat * vsqrt)
      Rw <- Rmat * vsqrt
      
      # form pseudo-response
      etaold <- eta
      yp <- (eta + (y - mu) / mueta) * vsqrt
      
      # form crossproduct matrices
      KtK <- crossprod(Kw)
      KtR <- crossprod(Kw, Rw)
      RtR <- crossprod(Rw)
      XtX <- rbind(cbind(KtK, KtR), cbind(t(KtR), RtR))
      Kty <- crossprod(Kw, yp)
      Rty <- crossprod(Rw, yp)
      Xty <- rbind(Kty, Rty)
      
      # update coefficients and eta
      XtXeig <- eigen(XtX + Pmat, symmetric = TRUE)
      nze <- sum(XtXeig$values > XtXeig$values[1] * .Machine$double.eps)
      XtXisqrt <- XtXeig$vectors[,1:nze] %*% diag(1/sqrt(XtXeig$values[1:nze]))
      XtXi <- tcrossprod(XtXisqrt)
      coef <- XtXi %*% Xty
      eta <- Xmat %*% coef
      
      # check if eta is valid
      while(!family$valideta(eta)){
        coef <- (coef + oldcoef) / 2
        eta <- Xmat %*% coef
      }
      
      # update mu and check for convergence
      mu <- family$linkinv(eta)
      devres <- family$dev.resids(y, mu, weights)
      dev <- sum(devres)
      famLL <- sum( family$logLik(y, nobs, mu, weights, dev) )
      vtol <- abs( (famLL - famLL.old) / (famLL.old + 1e-4) )
      iter <- iter + 1L
      famLL.old <- famLL
      oldcoef <- coef
      
    } # end while(iter < control$maxit & vtol > control$epsilon)
    
    # get degrees of freedom
    lev <- rowSums((Xmat %*% XtXisqrt)^2)
    df <- sum(vsqrt^2 * lev)
    
    # define loss value
    if(method == "GCV"){
      #wlev <- df / nobs
      #etacv <- (eta - wlev * yp / vsqrt) / (1 - wlev)
      etacv <- (eta - mean(lev) * yp * vsqrt) / (1 - df/nobs)
      mucv <- family$linkinv(etacv)
      theta <- family$theta(mucv)
      cumfun <- family$cumulant(mucv)
      val <- -mean(y * theta - cumfun)
    } else if(method == "OCV"){
      wlev <- vsqrt^2 * lev
      etacv <- (eta - wlev * yp / vsqrt) / (1 - wlev)
      mucv <- family$linkinv(etacv)
      theta <- family$theta(mucv)
      cumfun <- family$cumulant(mucv)
      val <- -mean(y * theta - cumfun)
    } else if(method == "GACV"){
      theta <- family$theta(mu)
      cumfun <- family$cumulant(mu)
      LLpart <- mean(y * theta - cumfun)
      CVpart <- (sum(lev) / (nobs - df)) * mean(y * (y - mu))
      val <- CVpart - LLpart
    } else if(method == "ACV"){
      theta <- family$theta(mu)
      cumfun <- family$cumulant(mu)
      LLpart <- mean(y * theta - cumfun)
      wlev <- vsqrt^2 * lev
      CVpart <- mean( (wlev / (1 - wlev)) * y * (y - mu) / vsqrt^2 )
      val <- CVpart - LLpart
    } else if(method == "AIC"){
      val <- 2 * (df - famLL)
    } else if(method == "BIC"){
      val <- log(nobs) * df - 2 * famLL
    }
    
    # attach attributes
    attr(val, "iter") <- iter
    attr(val, "df") <- df
    
    # return results
    if(return.coef) {
      if(family$family == "gaussian"){
        disp <- dev / (nobs - df)
      } else if(family$family == "binomial"){
        disp <- 1 / mean(weights)
      } else if(family$family == "poisson"){
        disp <- 1
      } else if(family$family == "Gamma"){
        disp <- dev / (nobs - df)
      } else if(family$family == "inverse.gaussian"){
        disp <- dev / (nobs - df)
      }
      return(list(objective = val, coef = as.numeric(coef),
                  cov.sqrt = sqrt(disp) * XtXisqrt, 
                  disp = disp, resid = as.numeric(yp / vsqrt - mu)))
    } else {
      return(val)
    }
    
  }

