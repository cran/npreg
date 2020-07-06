ss <-
  function(x, y = NULL, w = NULL, df, spar = NULL, lambda = NULL, 
           method = c("GCV", "OCV", "GACV", "ACV", "REML", "ML", "AIC", "BIC"),
           m = 2L, periodic = FALSE, all.knots = FALSE, nknots = .nknots.smspl, 
           knots = NULL, keep.data = TRUE, df.offset = 0, penalty = 1, 
           control.spar = list(), tol = 1e-6 * IQR(x)){
    # smoothing spline in R
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2020-06-27
    
    
    #########***#########   INITIAL CHECKS   #########***#########

    # check x and y
    if(is.null(y) | missing(y)){
      if(is.list(x)){
        if(length(x) != 2L) stop("If 'x' is a list, it must contain 2 elements (x and y).")
        y <- as.numeric(x[[2]])
        x <- as.numeric(x[[1]])
        n <- length(x)
        if(length(y) != n) stop("If 'x' is a list, the 2 elements must have length n.")
      } else if(is.matrix(x)){
        if(ncol(x) != 2L) stop("If 'x' is a matrix, it must contain 2 columns (x and y).")
        y <- as.numeric(x[,2])
        x <- as.numeric(x[,1])
        n <- length(x)
      } else {
        n <- length(x)
        y <- as.numeric(x)
        x <- 1:n
      }
      if(any(is.na(y)) | any(is.infinite(y)) | any(is.nan(y))) 
        stop("'x' and 'y' cannot contain missing (NA), infinite (Inf), or undefined (NaN) values")
    } else {
      x <- as.numeric(x)
      y <- as.numeric(y)
      n <- length(x)
      if(length(y) != n) stop("'x' and 'y' must have the same length")
      if(any(is.na(x)) | any(is.infinite(x)) | any(is.nan(x)))
        stop("'x' and 'y' cannot contain missing (NA), infinite (Inf), or undefined (NaN) values")
      if(any(is.na(y)) | any(is.infinite(y)) | any(is.nan(y))) 
        stop("'x' and 'y' cannot contain missing (NA), infinite (Inf), or undefined (NaN) values")
    }
    xmin <- min(x)
    xmax <- max(x)
    xrange <- xmax - xmin
    
    # check method
    methods <- c("GCV", "OCV", "GACV", "ACV", "REML", "ML", "AIC", "BIC")
    method <- as.character(method[1])
    method <- pmatch(toupper(method), methods)
    if(is.na(method)) stop("Invalid 'method' input.")
    method <- methods[method]
    
    # check w
    if(is.null(w)){
      w <- 1
      no.wghts <- TRUE
      if(any(method == c("REML", "ML"))) mliw <- 0
    } else {
      no.wghts <- FALSE
      w <- as.numeric(w)
      if(length(w) != n) stop("'w' must have the same length as inputs 'x' and 'y'")
      if(any(w < 0)) stop("'w' must contain non-negative weights")
      if(all(w == 0)) stop("some 'w' must be non-zero")
      w <- w * sum(w > 0) / sum(w)
      if(any(method == c("REML", "ML"))) mliw <- mean(log(1/w[w > 0]))
    }
    
    # keeping data?
    if(keep.data){
      if(no.wghts){
        data.orig <- data.frame(x = x, y = y)
      } else {
        data.orig <- data.frame(x = x, y = y, w = w)
      }
    } 
    
    # check m (penalty order)
    m <- as.integer(m[1])
    if(m < 1L | m > 3L) stop("'m' must be 1 (linear), 2 (cubic), or 3 (quintic)")
    
    # check df
    if(!missing(df)){
      df <- as.numeric(df[1])
      if(df < m | df > n) stop("'df' must satisfy:  m < df < n")
    }
    
    # check spar
    if(!is.null(spar)){
      spar <- as.numeric(spar[1])
      lambda <- 256^(3*(spar-1))
    }
    
    # check lambda
    if(!is.null(lambda)){
      lambda <- as.numeric(lambda[1])
      if(any(lambda < 0)) stop("'lambda' must satisfy:  lambda >= 0")
      spar <- 1 + log(lambda, base = 256) / 3
    }
    
    # check nknots
    if(!is.function(nknots)){
      nknots <- as.integer(nknots[1])
      if(nknots < 1L) stop("'nknots' must be a positive integer\n(or a function returning a positive integer)")
    }
    
    # check knots
    if(!is.null(knots)){
      knots <- sort(unique(as.numeric(knots)))
      nknots <- length(knots)
      if(nknots > n) stop("'knots' must satisfy:  length(knots) <= length(x)")
    }
    
    # check df.offset
    df.offset <- as.numeric(df.offset[1])
    if(df.offset < 0) stop("'df.offset' must be a non-negative scalar.")
    
    # check penalty
    penalty <- as.numeric(penalty[1])
    
    # check control parameters
    control.spar <- as.list(control.spar)
    if(is.null(control.spar$lower)) {
      control.spar$lower <- -1.5
    } else {
      control.spar$lower <- as.numeric(control.spar$lower[1])
    }
    if(is.null(control.spar$upper)){
      control.spar$upper <- 1.5
    } else {
      control.spar$upper <- as.numeric(control.spar$upper[1])
    }
    if(control.spar$upper <= control.spar$lower) stop("'control.spar$lower' and 'control.spar$upper' must satisfy:  low < high")
    if(is.null(control.spar$tol)) {
      control.spar$tol <- 1e-8
    } else {
      control.spar$tol <- as.numeric(control.spar$tol[1])
      if(control.spar$tol <= 0) stop("'control.spar$tol' must be positive")
    }
      
    # check tol
    tol <- as.numeric(tol[1])
    if(tol <= 0) stop("'tol' must be positive scalar")
    
    
    #########***#########   ROUNDING   #########***#########
    
    # find unique x
    xbar <- mean(x)
    rx <- round((x - xbar) / tol)
    xnd <- !duplicated(rx)
    nux <- sum(xnd)
    
    # sort by unique x
    ux.sort <- sort(x[xnd], index.return = TRUE)
    rux <- (rx[xnd])[ux.sort$ix]

    # get sufficient statistics
    if(nux < n){
      spux <- split(data.frame(w = w, wy = w * y, wy2 = w * y^2),
                    f = factor(rx, levels = rux))
      data <- as.data.frame(cbind(t(sapply(spux, colSums)), x = ux.sort$x), row.names = 1:nux)
      rm(spux)
    } else {
      data <- data.frame(w = w, wy = w * y, wy2 = w * y^2, x = x)[ux.sort$ix,]
    }
    yss <- sum(data$wy2)
    
    # remove junk
    rm(w, x, y, xbar, rx, xnd, ux.sort, rux)
    
    
    #########***#########   DETERMINE KNOTS   #########***#########
    
    if(is.null(knots)){
      if(all.knots){
        # use all unique x as knots
        knots <- data$x
        nknots <- nux
      } else {
        # use 'nknots' quantiles of x as knots
        if(is.function(nknots)){
          nknotfun <- nknots
          nknots <- as.integer(nknotfun(nux)[1])
        }
        if(nknots > nux) stop("Too many knots! Need length(unique(x)) >= length(unique(knots))")
        #knots <- quantile(data$x, probs = seq(0, 1, length.out = nknots + 1L)[1:nknots])
        knots <- data$x[seq(1, nrow(data), length.out = nknots)]
      } # end if(all.knots)
    } else {
      # use provided 'knots'
      if(nknots > nux) stop("Too many knots! Need length(unique(knots)) <= length(unique(x))")
      if(min(knots) < xmin | max(knots) > xmax) warning("Input 'knots' are outside of range of input 'x'")
    }
    
    
    #########***#########   BASIS AND PENALTY   #########***#########
    
    # make basis function matrix
    X <- cbind(1, basis_poly(x = data$x, knots = knots, 
                             m = m, xmin = xmin, xmax = xmax, 
                             periodic = periodic))
    nsdim <- ifelse(periodic, 1, m)
    Q <- penalty_poly(x = knots, m = m, xmin = xmin, 
                      xmax = xmax, periodic = periodic)
    
    # EVD of Q
    eps <- .Machine$double.eps
    Qeig <- eigen(Q, symmetric = TRUE)
    Qrnk <- sum(Qeig$values > nknots * eps * Qeig$values[1])
    if(Qrnk == 1L){
      Qprj <- matrix(Qeig$vectors[,1] / sqrt(Qeig$values[1]), nrow = nknots, ncol = 1)
    } else{
      Qprj <- Qeig$vectors[,1:Qrnk] %*% diag(1 / sqrt(Qeig$values[1:Qrnk]))
    }
    
    # reparameterize X
    nullindx <- 1:nsdim
    X.w <- sqrt(data$w) * X
    R.w <- X.w[,-nullindx,drop=FALSE] %*% Qprj
    XsvdN <- svd(X.w[,nullindx,drop=FALSE])
    X.w <- cbind(X.w[,nullindx], R.w - XsvdN$u %*% crossprod(XsvdN$u, R.w))
    
    # sse for null space
    y.w <- data$wy / sqrt(data$w)
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
      const <- mliw + (nval / n) + (n - m) * (1 + log(2 * pi) - log(n - m) ) / n
    } else if (method == "ML"){
      avec <- svd(R.w, nu = 0, nv = 0)$d^2
      const <- mliw + 1 + log(2*pi) - log(n)
    }
    
    # reverse transformation
    Tmat <- matrix(0, nsdim + nknots, nsdim + Qrnk)
    Tmat[nullindx,nullindx] <- diag(nsdim)
    Tmat[nullindx,-nullindx] <- (-1) * solve(crossprod(X.w[,nullindx])) %*% crossprod(X.w[,nullindx], R.w)
    Tmat[-nullindx,-nullindx] <- Qprj
    Tmat[,nullindx] <- Tmat[,nullindx,drop=FALSE] %*% XsvdN$v %*% diag(1 / XsvdN$d, nrow = nsdim, ncol = nsdim)
    Tmat[,-nullindx] <- Tmat[,-nullindx] %*% XsvdC$v %*% diag(1 / XsvdC$d)
    
    # remove junk
    rm(X, Qeig, Qprj, Qrnk, X.w, R.w)
    
    
    #########***#########   ESTIMATE COEFS   #########***#########
    
    # if lambda is provided, use that
    # else if spar is provided use that
    # else if df is provided use that
    # else tune via GCV, OCV, REML, or ML
    
    # tune lambda
    crit <- NA
    tunelambda <- (missing(df) && is.null(spar) && is.null(lambda))
    if(tunelambda){
      
      interval <- c(control.spar$lower, control.spar$upper)
      if(method == "GCV"){
        
        opt <- optimize(f = tune.gcv.ss, interval = interval,
                        bvec = bvec, dvec = dvec, n = n, yss = sse0,
                        nsdim = nsdim, df.offset = df.offset, penalty = penalty,
                        tol = control.spar$tol)
        spar <- opt$minimum
        lambda <- 256^(3*(spar-1))
        cv.crit <- opt$objective
        shrink <- 1 / (1 + n * lambda * dvec)
        beta <- bvec * shrink
        df <- sum(shrink) + nsdim
        sse <- cv.crit * n * (1 - (df.offset + penalty * df) / n)^2
        fit <- as.numeric(fit0 + XsvdC$u %*% beta) / sqrt(data$w)
        lev <- lev0 + rowSums((XsvdC$u %*% diag(sqrt(shrink)))^2)
        n2LL <- NULL
        sigma <- sqrt(sse / (n - df))
        beta <- c(beta0, beta)
        shrink <- c(rep(1, nsdim), shrink)
        
      } else if(method == "OCV") {
        
        opt <- optimize(f = tune.ocv.ss, interval = interval,
                        bvec = bvec, dvec = dvec, n = n, 
                        u = XsvdC$u, y = y.w, fit0 = fit0,
                        lev0 = lev0, tol = control.spar$tol)
        spar <- opt$minimum
        lambda <- 256^(3*(spar-1))
        cv.crit <- opt$objective
        shrink <- 1 / (1 + n * lambda * dvec)
        beta <- bvec * shrink
        fit <- as.numeric(fit0 + XsvdC$u %*% beta) / sqrt(data$w)
        df <- sum(shrink) + nsdim
        sse <- sse0 - 2 * sum(bvec * beta) + sum(beta^2)
        lev <- lev0 + rowSums((XsvdC$u %*% diag(sqrt(shrink)))^2)
        n2LL <- NULL
        sigma <- sqrt(sse / (n - df))
        beta <- c(beta0, beta)
        shrink <- c(rep(1, nsdim), shrink)
        
      } else if(method == "GACV") {
        
        opt <- optimize(f = tune.gacv.ss, interval = interval,
                        bvec = bvec, dvec = dvec, n = n, yss = sse0,
                        nsdim = nsdim, const = yss/2, tol = control.spar$tol)
        spar <- opt$minimum
        lambda <- 256^(3*(spar-1))
        cv.crit <- opt$objective
        shrink <- 1 / (1 + n * lambda * dvec)
        beta <- bvec * shrink
        fit <- as.numeric(fit0 + XsvdC$u %*% beta) / sqrt(data$w)
        df <- sum(shrink) + nsdim
        sse <- sse0 - 2 * sum(bvec * beta) + sum(beta^2)
        lev <- lev0 + rowSums((XsvdC$u %*% diag(sqrt(shrink)))^2)
        n2LL <- NULL
        sigma <- sqrt(sse / (n - df))
        beta <- c(beta0, beta)
        shrink <- c(rep(1, nsdim), shrink)
        
      } else if(method == "ACV") {
        
        opt <- optimize(f = tune.acv.ss, interval = interval,
                        bvec = bvec, dvec = dvec, n = n, yss = sse0,
                        u = XsvdC$u, y = y.w, fit0 = fit0,
                        lev0 = lev0, const = yss/2, tol = control.spar$tol)
        spar <- opt$minimum
        lambda <- 256^(3*(spar-1))
        cv.crit <- opt$objective
        shrink <- 1 / (1 + n * lambda * dvec)
        beta <- bvec * shrink
        fit <- as.numeric(fit0 + XsvdC$u %*% beta) / sqrt(data$w)
        df <- sum(shrink) + nsdim
        sse <- sse0 - 2 * sum(bvec * beta) + sum(beta^2)
        lev <- lev0 + rowSums((XsvdC$u %*% diag(sqrt(shrink)))^2)
        n2LL <- NULL
        sigma <- sqrt(sse / (n - df))
        beta <- c(beta0, beta)
        shrink <- c(rep(1, nsdim), shrink)
        
      } else if(any(method == c("REML", "ML"))){
        
        opt <- optimize(f = tune.mle.ss, interval = interval,
                        avec = avec, bvec = bvec, dvec = dvec, 
                        n = n, yss = yss - sum(y.w * fit0), 
                        m = ifelse(method == "REML", nsdim, 0),
                        const = const, tol = control.spar$tol)
        spar <- opt$minimum
        lambda <- 256^(3*(spar-1))
        n2LL <- n * opt$objective
        shrink <- 1 / (1 + n * lambda * dvec)
        beta <- bvec * shrink
        sse <- sse0 - 2 * sum(bvec * beta) + sum(beta^2)
        df <- sum(shrink) + nsdim
        cv.crit <- (sse / n) / (1 - (df.offset + penalty * df) / n)^2
        fit <- as.numeric(fit0 + XsvdC$u %*% beta) / sqrt(data$w)
        lev <- lev0 + rowSums((XsvdC$u %*% diag(sqrt(shrink)))^2)
        sseML <- yss - sum(y.w * fit0) - sum(bvec * beta)
        sigma <- sqrt(sseML / (n - ifelse(method == "REML", nsdim, 0)))
        beta <- c(beta0, beta)
        shrink <- c(rep(1, nsdim), shrink)
        
      } else if(any(method == c("AIC", "BIC"))){
        
        const <- ifelse(method == "AIC", 2, log(n))
        opt <- optimize(f = tune.aic.ss, interval = interval,
                        bvec = bvec, dvec = dvec, n = n, yss = sse0,
                        nsdim = nsdim, const = const,
                        tol = control.spar$tol)
        spar <- opt$minimum
        lambda <- 256^(3*(spar-1))
        ic <- opt$objective + n * (1 + log(2*pi))
        shrink <- 1 / (1 + n * lambda * dvec)
        beta <- bvec * shrink
        df <- sum(shrink) + nsdim
        sse <- sse0 - 2 * sum(bvec * beta) + sum(beta^2)
        cv.crit <- (sse / n) / (1 - (df.offset + penalty * df) / n)^2
        fit <- as.numeric(fit0 + XsvdC$u %*% beta) / sqrt(data$w)
        lev <- lev0 + rowSums((XsvdC$u %*% diag(sqrt(shrink)))^2)
        n2LL <- ic - const * df
        sigma <- sqrt(sse / (n - df))
        beta <- c(beta0, beta)
        shrink <- c(rep(1, nsdim), shrink)
        
      } # end if(method == "GCV")
      
    } else {
      
      # convert df to lambda
      if(is.null(lambda)){
        interval <- 256^(3 * c(control.spar$lower - 1, control.spar$upper - 1))
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
      fit <- as.numeric(fit0 + XsvdC$u %*% beta) / sqrt(data$w)
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
        part3 <- (n - ifelse(method == "REML", m, 0)) * log(sseML)
        n2LL <- part1 - part2 + part3 + n * const
        sigma <- sqrt(sseML / (n - ifelse(method == "REML", nsdim, 0)))
      } else if(any(method == c("AIC", "BIC"))){
        cv.crit <- (sse / n) / (1 - (df.offset + penalty * df) / n)^2
        const <- ifelse(method == "AIC", 2, log(n))
        n2LL <- n * log(sse / n) + n * (1 + log(2*pi))
        ic <- n2LL + const * df
      }
      if(is.null(sigma)) sigma <- sqrt(sse / (n - df))
      beta <- c(beta0, beta)
      shrink <- c(rep(1, nsdim), shrink)
      
    } # end if(tunelambda)
    
    
    #########***#########   RETURN RESULTS   #########***#########
    
    # evaluate penalty and retransform coefficients
    if(!periodic & m > 1){
      nlab <- "x"
      if(m == 3L) nlab <- c(nlab, "x^2")
    } else {
      nlab <- NULL
    }
    coefnames <- c("(Intercept)", nlab, paste0("s(x).", 1:nknots))
    beta <- as.numeric(Tmat %*% beta)
    names(beta) <- coefnames
    penalty <- as.numeric(crossprod(beta[-nullindx], Q %*% beta[-nullindx]))
    
    # coef cov sqrt
    cov.sqrt <- sigma * Tmat %*% diag(sqrt(shrink))
    
    # collect results
    fitinfo <- list(n = n, knot = knots, nk = nsdim + nknots, coef = beta,
                    min = xmin, range = xmax - xmin, m = m, periodic = periodic,
                    cov.sqrt = cov.sqrt, weighted = !no.wghts)
    ss <- list(x = data$x, y = fit, w = data$w, yin = data$wy / ifelse(data$w > 0, data$w, 1),
               tol = tol, data = if(keep.data) data.orig, lev = lev, 
               cv.crit = cv.crit, pen.crit = sse, crit = crit, df = df, 
               spar = spar, lambda = lambda, fit = fitinfo, call = match.call(),
               sigma = sigma, logLik = if(!is.null(n2LL)) (-1/2) * n2LL,
               aic = if(method == "AIC") ic, bic = if(method == "BIC") ic,
               penalty = penalty, method = method)
    class(ss) <- "ss"
    return(ss)
    
  } # end ss.R

# print function
print.ss <-
  function(x, ...){
    cat("\nCall:\n")
    print(x$call)
    cat("\nSmoothing Parameter  spar =", x$spar, "  lambda =", x$lambda)
    cat("\nEquivalent Degrees of Freedom (Df)", x$df)
    if(x$fit$weighted){
      cat("\nPenalized Criterion (weighted RSS)", x$pen.crit)
    } else {
      cat("\nPenalized Criterion (RSS)", x$pen.crit)
    }
    if(x$method == "GCV"){
      cat("\nGeneralized Cross-Validation (GCV)", x$cv.crit,"\n\n")
    } else if(x$method == "OCV"){
      cat("\nOrdinary Cross-Validation (OCV)", x$cv.crit,"\n\n")
    } else if(x$method == "GACV"){
      cat("\nGeneralized Approximate Cross-Validation (GACV)", x$cv.crit,"\n\n")
    } else if(x$method == "ACV"){
      cat("\nApproximate Cross-Validation (ACV)", x$cv.crit,"\n\n")
    } else if(x$method == "REML"){
      cat("\nLog Likelihood (REML)", x$logLik,"\n\n")
    } else if(x$method == "ML"){
      cat("\nLog Likelihood (ML)", x$logLik,"\n\n")
    } else if(x$method == "AIC"){
      cat("\nAkaike's Information Criterion (AIC)", x$aic,"\n\n")
    } else if(x$method == "BIC"){
      cat("\nBayesian Information Criterion (BIC)", x$bic,"\n\n")
    }
    
  } # end print.ss

# gcv tuning function
tune.gcv.ss <-
  function(spar, bvec, dvec, n, yss, nsdim = 1, df.offset = 0, penalty = 1){
    lambda <- 256^(3*(spar-1))
    shrink <- 1 / (1 + n * lambda * dvec)
    beta <- bvec * shrink
    df <- sum(shrink) + nsdim
    sse <- yss - 2 * sum(bvec * beta) + sum(beta^2)
    (sse / n) / (1 - (df.offset + penalty * df) / n)^2
  } # end tune.gcv.ss.R

# ocv tuning function
tune.ocv.ss <-
  function(spar, bvec, dvec, n, u, y, fit0, lev0){
    lambda <- 256^(3*(spar-1))
    shrink <- 1 / (1 + n * lambda * dvec)
    beta <- bvec * shrink
    fit <- fit0 + u %*% beta
    lev <- lev0 + rowSums((u %*% diag(sqrt(shrink)))^2)
    mean(((y - fit) / (1 - lev))^2)
  } # end tune.ocv.ss.R

# gacv tuning function
tune.gacv.ss <-
  function(spar, bvec, dvec, n, yss, const = 0, nsdim = 1){
    lambda <- 256^(3*(spar-1))
    shrink <- 1 / (1 + n * lambda * dvec)
    beta <- bvec * shrink
    df <- sum(shrink) + nsdim
    sse1 <- yss - sum(bvec * beta)
    sse2 <- yss + sum(bvec^2 * (shrink^2 - 2*shrink))
    ( (1/2) * sse2  + (df / (n - df)) * sse1 - const ) / n
  } # end tune.gacv.ss.R

# acv tuning function
tune.acv.ss <-
  function(spar, bvec, dvec, n, yss, u, y, fit0, lev0, const = 0){
    lambda <- 256^(3*(spar-1))
    shrink <- 1 / (1 + n * lambda * dvec)
    beta <- bvec * shrink
    fit <- fit0 + u %*% beta
    lev <- lev0 + rowSums((u %*% diag(sqrt(shrink)))^2)
    sse1 <- sum((lev / (1 - lev)) * y * (y - fit))
    sse2 <- yss + sum(bvec^2 * (shrink^2 - 2*shrink))
    ( (1/2) * sse2  + sse1 - const ) / n
  } # end tune.acv.ss.R

# reml/ml tuning function
tune.mle.ss <-
  function(spar, avec, bvec, dvec, n, yss, m = 0, const = 0){
    lambda <- 256^(3*(spar-1))
    nlam <- n * lambda
    shrink <- 1 / (1 + nlam * dvec)
    beta <- bvec * shrink
    sse <- yss - sum(bvec * beta)
    r <- length(avec)
    part1 <- sum(log(nlam + avec)) / n
    part2 <- (r / n) * log(nlam)
    part3 <- ((n - m) / n) * log(sse)
    part1 - part2 + part3 + const
  } # end tune.mle.ss.R

# aic/bic tuning function
tune.aic.ss <-
  function(spar, bvec, dvec, n, yss, nsdim = 1, const = 2){
    lambda <- 256^(3*(spar-1))
    shrink <- 1 / (1 + n * lambda * dvec)
    beta <- bvec * shrink
    df <- sum(shrink) + nsdim
    sse <- yss - 2 * sum(bvec * beta) + sum(beta^2)
    n * log(sse / n) + const * df
  } # end tune.aic.ss.R

df2lambda <- 
  function(lambda, df, n, nsdim, dvec){
    (df - nsdim - sum(1 / (1 + n * lambda * dvec)))^2
  } # end df2lambda.R
