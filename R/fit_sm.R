fit_sm <-
  function(y, depe, df = NULL, lambda = NULL, 
           tprk = TRUE, control = list(), method = "GCV"){
    # fit semiparametric model
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2022-03-22
    
    
    #########***#########   REPARAMETERIZATION   #########***#########
    
    # info
    n <- nrow(depe$K)
    nsdim <- ncol(depe$K)
    if(!tprk){
      Nknots <- sapply(depe$J, ncol)
      depe$J <- do.call(cbind, depe$J)
    }
    nknots <- ncol(depe$J)
    nullindx <- 1:nsdim
    
    # weighted matrices
    wsqrt <- sqrt(depe$weights)
    y.w <- y * wsqrt
    yss <- sum(y.w^2)
    depe$K <- depe$K * wsqrt
    depe$J <- depe$J * wsqrt
    
    # reparameterize contrast space
    if(tprk){
      
      Qisqrt <- msqrt(depe$Q, inverse = TRUE, checkx = FALSE)
      R.w <- depe$J %*% Qisqrt
      Qrnk <- ncol(Qisqrt)
      
    } else {
      
      cknots <- c(0, cumsum(Nknots))
      R.w <- Qisqrt <- vector("list", length(depe$Q))
      for(k in 1:length(depe$Q)){
        indx <- seq(cknots[k] + 1, cknots[k+1])
        Qisqrt[[k]] <- msqrt(depe$Q[[k]], inverse = TRUE, checkx = FALSE)
        R.w[[k]] <- depe$J[,indx] %*% Qisqrt[[k]]
      }
      R.w <- do.call("cbind", R.w)
      Qrnk <- ncol(R.w)
      
    } # end if(tprk)
    
    # define X
    XsvdN <- svd(depe$K)
    X.w <- cbind(depe$K, R.w - XsvdN$u %*% crossprod(XsvdN$u, R.w))
    
    # sse for null space
    beta0 <- crossprod(XsvdN$u, y.w)
    fit0 <- as.numeric(XsvdN$u %*% beta0)
    sse0 <- yss - 2 * sum(y.w * fit0) + sum(beta0^2)
    lev0 <- rowSums(XsvdN$u^2)
    
    # SVD of weighted X
    XsvdC <- svd(X.w[,-nullindx,drop=FALSE])
    bvec <- crossprod(XsvdC$u, y.w)
    dvec <- 1/XsvdC$d^2
    
    # SVD of ML or REML
    if(method == "REML") {
      avec <- XsvdC$d^2
      nval <- sum(log(XsvdN$d^2))
      mliw <- mean(log(1/depe$weights[depe$weights > 0]))
      const <- mliw + (nval / n) + (n - nsdim) * (1 + log(2 * pi) - log(n - nsdim) ) / n
    } else if (method == "ML"){
      avec <- svd(R.w, nu = 0, nv = 0)$d^2
      mliw <- mean(log(1/depe$weights[depe$weights > 0]))
      const <- mliw + 1 + log(2*pi) - log(n)
    }
    
    # reverse transformation
    Tmat <- matrix(0, nsdim + nknots, nsdim + Qrnk)
    Tmat[nullindx,nullindx] <- diag(nsdim)
    Tmat[nullindx,-nullindx] <- (-1) * solve(crossprod(X.w[,nullindx])) %*% crossprod(X.w[,nullindx], R.w)
    if(tprk){
      Tmat[-nullindx,-nullindx] <- Qisqrt
    } else {
      row.offset <- col.offset <- nsdim
      for(k in 1:length(Qisqrt)){
        nrowk <- nrow(Qisqrt[[k]])
        ncolk <- ncol(Qisqrt[[k]])
        Tmat[row.offset + 1:nrowk, col.offset + 1:ncolk] <- Qisqrt[[k]]
        row.offset <- row.offset + nrowk
        col.offset <- col.offset + ncolk
      }
    }
    Tmat[,nullindx] <- Tmat[,nullindx] %*% XsvdN$v %*% diag(1 / XsvdN$d, nrow = nsdim, ncol = nsdim)
    #Tmat[,-nullindx] <- Tmat[,-nullindx] %*% XsvdC$v %*% diag(1 / XsvdC$d)
    Tmat <- cbind(Tmat[,nullindx], Tmat[,-nullindx] %*% XsvdC$v %*% diag(1 / XsvdC$d))
    
    # get coefficient names
    coefnames <- c(colnames(depe$K), colnames(depe$J))
    
    # remove junk
    Q <- depe$Q
    rm(depe, Qisqrt, X.w, R.w)
    
    
    #########***#########   ESTIMATE COEFS   #########***#########
    
    df.offset <- 0
    penalty <- 1
    if(is.null(df) && is.null(lambda)){
      interval <- c(control$lower, control$upper)
      if(method == "GCV"){
        
        opt <- optimize(f = tune.gcv.ss, interval = interval,
                        bvec = bvec, dvec = dvec, n = n, yss = sse0,
                        nsdim = nsdim, df.offset = df.offset, penalty = penalty,
                        tol = control$tol)
        spar <- opt$minimum
        lambda <- 256^(3*(spar-1))
        cv.crit <- opt$objective
        shrink <- 1 / (1 + n * lambda * dvec)
        beta <- bvec * shrink
        df <- sum(shrink) + nsdim
        sse <- cv.crit * n * (1 - (df.offset + penalty * df) / n)^2
        sigma <- sqrt(sse / (n - df))
        fit <- as.numeric(fit0 + XsvdC$u %*% beta) / wsqrt
        se.fit <- sigma * sqrt(lev0 + rowSums((XsvdC$u %*% diag(sqrt(shrink)))^2)) / wsqrt
        beta <- c(beta0, beta)
        shrink <- c(rep(1, nsdim), shrink)
        n2LL <- NULL
        
      } else if(method == "OCV") {
        
        opt <- optimize(f = tune.ocv.ss, interval = interval,
                        bvec = bvec, dvec = dvec, n = n, 
                        u = XsvdC$u, y = y.w, fit0 = fit0,
                        lev0 = lev0, tol = control$tol)
        spar <- opt$minimum
        lambda <- 256^(3*(spar-1))
        cv.crit <- opt$objective
        shrink <- 1 / (1 + n * lambda * dvec)
        beta <- bvec * shrink
        fit <- as.numeric(fit0 + XsvdC$u %*% beta) / wsqrt
        df <- sum(shrink) + nsdim
        sse <- sse0 - 2 * sum(bvec * beta) + sum(beta^2)
        sigma <- sqrt(sse / (n - df))
        se.fit <- sigma * sqrt(lev0 + rowSums((XsvdC$u %*% diag(sqrt(shrink)))^2)) / wsqrt
        beta <- c(beta0, beta)
        shrink <- c(rep(1, nsdim), shrink)
        n2LL <- NULL
        
      } else if(method == "GACV") {
        
        opt <- optimize(f = tune.gacv.ss, interval = interval,
                        bvec = bvec, dvec = dvec, n = n, yss = sse0,
                        nsdim = nsdim, const = yss/2, tol = control$tol)
        spar <- opt$minimum
        lambda <- 256^(3*(spar-1))
        cv.crit <- opt$objective
        shrink <- 1 / (1 + n * lambda * dvec)
        beta <- bvec * shrink
        fit <- as.numeric(fit0 + XsvdC$u %*% beta) / wsqrt
        df <- sum(shrink) + nsdim
        sse <- sse0 - 2 * sum(bvec * beta) + sum(beta^2)
        sigma <- sqrt(sse / (n - df))
        se.fit <- sigma * sqrt(lev0 + rowSums((XsvdC$u %*% diag(sqrt(shrink)))^2)) / wsqrt
        n2LL <- NULL
        beta <- c(beta0, beta)
        shrink <- c(rep(1, nsdim), shrink)
        
      } else if(method == "ACV") {
        
        opt <- optimize(f = tune.acv.ss, interval = interval,
                        bvec = bvec, dvec = dvec, n = n, yss = sse0,
                        u = XsvdC$u, y = y.w, fit0 = fit0,
                        lev0 = lev0, const = yss/2, tol = control$tol)
        spar <- opt$minimum
        lambda <- 256^(3*(spar-1))
        cv.crit <- opt$objective
        shrink <- 1 / (1 + n * lambda * dvec)
        beta <- bvec * shrink
        fit <- as.numeric(fit0 + XsvdC$u %*% beta) / wsqrt
        df <- sum(shrink) + nsdim
        sse <- sse0 - 2 * sum(bvec * beta) + sum(beta^2)
        sigma <- sqrt(sse / (n - df))
        se.fit <- sigma * sqrt(lev0 + rowSums((XsvdC$u %*% diag(sqrt(shrink)))^2)) / wsqrt
        n2LL <- NULL
        beta <- c(beta0, beta)
        shrink <- c(rep(1, nsdim), shrink)
        
      } else if(any(method == c("REML", "ML"))){
        
        opt <- optimize(f = tune.mle.ss, interval = interval,
                        avec = avec, bvec = bvec, dvec = dvec, 
                        n = n, yss = yss - sum(y.w * fit0), 
                        m = ifelse(method == "REML", nsdim, 0),
                        const = const, tol = control$tol)
        spar <- opt$minimum
        lambda <- 256^(3*(spar-1))
        n2LL <- n * opt$objective
        shrink <- 1 / (1 + n * lambda * dvec)
        beta <- bvec * shrink
        sse <- sse0 - 2 * sum(bvec * beta) + sum(beta^2)
        sseML <- yss - sum(y.w * fit0) - sum(bvec * beta)
        sigma <- sqrt(sseML / (n - ifelse(method == "REML", nsdim, 0)))
        df <- sum(shrink) + nsdim
        cv.crit <- (sse / n) / (1 - (df.offset + penalty * df) / n)^2
        fit <- as.numeric(fit0 + XsvdC$u %*% beta) / wsqrt
        se.fit <- sigma * sqrt(lev0 + rowSums((XsvdC$u %*% diag(sqrt(shrink)))^2)) / wsqrt
        beta <- c(beta0, beta)
        shrink <- c(rep(1, nsdim), shrink)
        
      } else if(any(method == c("AIC", "BIC"))){
        
        const <- ifelse(method == "AIC", 2, log(n))
        opt <- optimize(f = tune.aic.ss, interval = interval,
                        bvec = bvec, dvec = dvec, n = n, yss = sse0,
                        nsdim = nsdim, const = const,
                        tol = control$tol)
        spar <- opt$minimum
        lambda <- 256^(3*(spar-1))
        ic <- opt$objective + n * (1 + log(2*pi))
        shrink <- 1 / (1 + n * lambda * dvec)
        beta <- bvec * shrink
        df <- sum(shrink) + nsdim
        sse <- sse0 - 2 * sum(bvec * beta) + sum(beta^2)
        cv.crit <- (sse / n) / (1 - (df.offset + penalty * df) / n)^2
        fit <- as.numeric(fit0 + XsvdC$u %*% beta) / wsqrt
        n2LL <- ic - const * df
        sigma <- sqrt(sse / (n - df))
        se.fit <- sigma * sqrt(lev0 + rowSums((XsvdC$u %*% diag(sqrt(shrink)))^2)) / wsqrt
        beta <- c(beta0, beta)
        shrink <- c(rep(1, nsdim), shrink)
        
      } # end if(method == "GCV")
      
    } else {
      
      # convert df to lambda
      if(is.null(lambda)){
        interval <- 256^(3 * c(control$lower - 1, control$upper - 1))
        getlam <- optimize(df2lambda, interval = interval, df = df, n = n, 
                           nsdim = nsdim, dvec = dvec, tol = .Machine$double.eps)
        lambda <- getlam$minimum
        crit <- getlam$objective
      } # end if(is.null(lambda))
      
      # fitting with lambda
      shrink <- 1 / (1 + n * lambda * dvec)
      beta <- bvec * shrink
      df <- sum(shrink) + nsdim
      sse <- sse0 - 2 * sum(bvec * beta) + sum(beta^2)
      fit <- as.numeric(fit0 + XsvdC$u %*% beta) / wsqrt
      lev <- lev0 + rowSums((XsvdC$u %*% diag(sqrt(shrink)))^2)
      n2LL <- sigma <- NULL
      if(method == "GCV"){
        cv.crit <- (sse / n) / (1 - (df.offset + penalty * df) / n)^2
      } else if(method == "OCV"){
        cv.crit <- mean(((y - fit) / (1 - lev))^2)
      } else if(method == "GACV"){
        sse1 <- sse0 - sum(bvec * beta)
        sse2 <- sse0 + sum(bvec^2 * (shrink^2 - 2*shrink))
        cv.crit <- ( (1/2) * sse2  + (df / (n - df)) * sse1 - yss/2 ) / n
      } else if(method == "ACV"){
        sse1 <- sum((lev / (1 - lev)) * y * (y - fit))
        sse2 <- sse0 + sum(bvec^2 * (shrink^2 - 2*shrink))
        cv.crit <- ( (1/2) * sse2  + sse1 - yss/2 ) / n
      } else if(any(method == c("REML", "ML"))){
        cv.crit <- (sse / n) / (1 - (df.offset + penalty * df) / n)^2
        sseML <- yss - sum(y.w * fit0) - sum(bvec * beta)
        r <- length(avec)
        part1 <- sum(log(n * lambda + avec))
        part2 <- r * log(n * lambda)
        part3 <- (n - ifelse(method == "REML", nsdim, 0)) * log(sseML)
        n2LL <- part1 - part2 + part3 + n * const
        sigma <- sqrt(sseML / (n - ifelse(method == "REML", nsdim, 0)))
      } else if(any(method == c("AIC", "BIC"))){
        cv.crit <- (sse / n) / (1 - (df.offset + penalty * df) / n)^2
        const <- ifelse(method == "AIC", 2, log(n))
        n2LL <- n * log(sse / n) + n * (1 + log(2*pi))
        ic <- n2LL + const * df
      }
      if(is.null(sigma)) sigma <- sqrt(sse / (n - df))
      se.fit <- sigma * sqrt(lev) / wsqrt
      beta <- c(beta0, beta)
      shrink <- c(rep(1, nsdim), shrink)
      
    } # end if(is.null(lambda))
    
    # retransform coefficients
    beta <- as.numeric(Tmat %*% beta)
    names(beta) <- coefnames
    
    # coef cov sqrt
    cov.sqrt <- sigma * Tmat %*% diag(sqrt(shrink))
    
    # smoothness penalty
    if(tprk){
      penalty <- as.numeric(crossprod(beta[-nullindx], Q %*% beta[-nullindx]))
    } else {
      penalty <- 0
      scoefs <- beta[-nullindx]
      for(k in 1:length(Q)){
        indx <- seq(cknots[k] + 1, cknots[k+1])
        penalty <- penalty + as.numeric(crossprod(scoefs[indx], Q[[k]] %*% scoefs[indx]))
      }
    }
    
    # R-squared
    r.squared <- as.numeric(cor(y, fit)^2)
    
    
    #########***#########   RETURN RESULTS   #########***#########
    
    res <- list(fitted.values = fit, se.fit = se.fit,
                sse = sse, cv.crit = cv.crit, 
                df = df, nsdf = nsdim, 
                r.squared = r.squared, sigma = sigma, 
                logLik = if(!is.null(n2LL)) (-1/2) * n2LL,
                aic = if(method == "AIC") ic, bic = if(method == "BIC") ic,
                spar = 1 + log(lambda, base = 256) / 3, lambda = lambda, 
                penalty = penalty, coefficients = beta, cov.sqrt = cov.sqrt)
    return(res)
    
  } # end fit_sm