boot.gsm <-
  function(object, statistic, ..., R = 9999, level = 0.95, bca = TRUE, 
           method = c("cases", "resid", "param"), fix.lambda = TRUE, 
           fix.thetas = TRUE, cov.mat = FALSE, boot.dist = FALSE, 
           verbose = TRUE, parallel = FALSE, cl = NULL){
    # Bootstrap a fit Generalized Smooth Model
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Update: 2022-03-21
    
    
    
    #########***#########   INITIAL CHECKS   #########***#########
    
    # check 'object'
    if(class(object) != "gsm") stop("Input 'object' must be of class 'gsm'")
    if(is.null(object$data)) stop("Input 'object' has no data. Need to refit model with data.")
    n <- nrow(object$data)
    
    # check 'statistic'
    if(missing(statistic)){
      nostat <- TRUE
      statistic <- function(object, ...) predict(object)
    } else {
      nostat <- FALSE
      statistic <- as.function(statistic)
    }
    t0 <- statistic(object, ...)
    if(!is.vector(t0))
      stop("Output of the 'statistic' function must be a vector")
    nstat <- length(t0)
    varnames <- names(t0)
    
    # check 'R'
    R <- as.integer(R[1])
    if(R < 1L) stop("Input 'R' must be a positive integer")
    nboot <- R + 1L
    
    # check 'level'
    level <- as.numeric(level[1])
    if(level <= 0 | level >= 1) 
      stop("Input 'level' must be between 0 and 1")
    nlevel <- length(level)
    alpha <- 1 - level
    
    # check 'bca'
    bca <- as.logical(bca[1])
    if(!any(bca == c(TRUE, FALSE))) stop("Input 'bca' must be TRUE/FALSE")
    
    # check 'method'
    method <- as.character(method[1])
    method <- pmatch(tolower(method), c("cases", "resid", "param"))
    if(is.na(method)) stop("Input method must be cases, resid, or param")
    method <- c("cases", "resid", "param")[method]
    
    # check 'fix.lambda'
    if(fix.lambda){
      lambda <- object$lambda
    } else {
      lambda <- NULL
    }
    if(method == "resid" && !fix.lambda) {
      fix.lambda <- TRUE
      lambda <- object$lambda
      warning("fix.lambda = TRUE is required when method = 'resid'.\nBootstrapping was conducted with lambda fixed.")
    }
    
    # check 'fix.thetas'
    if(fix.thetas){
      thetas <- object$specs$thetas
    } else {
      thetas <- NULL
    }
    if(method == "resid" && !fix.thetas) {
      fix.thetas <- TRUE
      thetas <- object$specs$thetas
      warning("fix.thetas = TRUE is required when method = 'resid'.\nBootstrapping was conducted with thetas fixed.")
    }
    
    # check 'parallel' and 'cl'
    make.cl <- FALSE
    if(parallel){
      if(is.null(cl)){
        make.cl <- TRUE
        ncores <- parallel::detectCores()
        cl <- parallel::makeCluster(ncores)
      } else {
        if(!any(class(cl) == "cluster")) stop("Input 'cl' must be an object of class 'cluster'.")
      }
    }
    
    
    
    ######***######   BOOTSTRAP   ######***######
    
    # initialize bootstrap distribution
    bootdist <- matrix(0, nboot, nstat)
    colnames(bootdist) <- varnames
    bootdist[1,] <- t0
    
    # setup progress bar
    if(verbose && !parallel){
      pbar <- txtProgressBar(min = 1, max = nboot + ifelse(bca, n, 0), style = 3)
    }
    
    # create bootstrap distribution
    if(method == "cases"){
      
      if(parallel){
        
        parfun <- function(x){
          index <- sample.int(n, replace = TRUE)
          bootdata <- object$data[index,]
          temp.object <- gsm(formula = object$formula, family = object$family,
                             data = bootdata, types = as.list(object$types), 
                             tprk = object$specs$tprk, knots = object$specs$knots, 
                             skip.iter = object$specs$skip.iter, lambda = lambda, 
                             control = object$specs$control, method = object$method, 
                             xrange = object$specs$xrng, thetas = thetas, mf = bootdata)
          if(nostat) temp.object$data <- object$data
          statistic(temp.object, ...)
        }
        bootdist <- t(cbind(t0, parallel::parSapply(cl = cl, X = 2:nboot, FUN = parfun)))
        
      } else {
        
        for(i in 2:nboot){
          index <- sample.int(n, replace = TRUE)
          bootdata <- object$data[index,]
          temp.object <- gsm(formula = object$formula, family = object$family,
                             data = bootdata, types = as.list(object$types), 
                             tprk = object$specs$tprk, knots = object$specs$knots, 
                             skip.iter = object$specs$skip.iter, lambda = lambda, 
                             control = object$specs$control, method = object$method, 
                             xrange = object$specs$xrng, thetas = thetas, mf = bootdata)
          if(nostat) temp.object$data <- object$data
          bootdist[i,] <- statistic(temp.object, ...)
          if (verbose) setTxtProgressBar(pbar, i)
          
        } # end for(i in 2:nboot)
        
      } # end if(parallel)
      
    } else if(method == "resid"){
      
      xinv <- tcrossprod(object$cov.sqrt) %*% t(model.matrix(object)) / object$dispersion
      temp.object <- object
      
      if(parallel){
        
        parfun <- function(x){
          index <- sample.int(n, replace = TRUE)
          pseudoy <- object$linear.predictors + object$residuals[index]
          temp.object$coefficients <- xinv %*% pseudoy
          names(temp.object$coefficients) <- names(object$coefficients)
          statistic(temp.object, ...)
        }
        bootdist <- t(cbind(t0, parallel::parSapply(cl = cl, X = 2:nboot, FUN = parfun)))
        
      } else {
        
        for(i in 2:nboot){
          index <- sample.int(n, replace = TRUE)
          pseudoy <- object$linear.predictors + object$residuals[index]
          temp.object$coefficients <- xinv %*% pseudoy
          names(temp.object$coefficients) <- names(object$coefficients)
          bootdist[i,] <- statistic(temp.object, ...)
          if (verbose) setTxtProgressBar(pbar, i)
        } # end for(i in 2:nboot)
        
      } # end if(parallel)
      
    } else {
      
      bootdata <- object$data
      
      if(parallel){
        
        parfun <- function(x){
          
          bootdata[,1] <- rexpfam(object)
          temp.object <- gsm(formula = object$formula, family = object$family,
                             data = bootdata, types = as.list(object$types), 
                             tprk = object$specs$tprk, knots = object$specs$knots, 
                             skip.iter = object$specs$skip.iter, lambda = lambda, 
                             control = object$specs$control, method = object$method, 
                             xrange = object$specs$xrng, thetas = thetas, mf = bootdata)
          if(nostat) temp.object$data <- object$data
          statistic(temp.object, ...)
        } # end parfun
        bootdist <- t(cbind(t0, parallel::parSapply(cl = cl, X = 2:nboot, FUN = parfun)))
        
      } else {
        
        for(i in 2:nboot){
          bootdata[,1] <- rexpfam(object)
          temp.object <- gsm(formula = object$formula, family = object$family,
                             data = bootdata, types = as.list(object$types), 
                             tprk = object$specs$tprk, knots = object$specs$knots, 
                             skip.iter = object$specs$skip.iter, lambda = lambda, 
                             control = object$specs$control, method = object$method, 
                             xrange = object$specs$xrng, thetas = thetas, mf = bootdata)
          if(nostat) temp.object$data <- object$data
          bootdist[i,] <- statistic(temp.object, ...)
          if (verbose) setTxtProgressBar(pbar, i)
        } # end for(i in 2:nboot)
        
      } # end if(parallel)
      
    } # end if(method == "cases")
    
    
    ######***######   CONFIDENCE INTERVALS   ######***######
    
    # bca interval?
    if(bca){
      
      # initialize jackknife distribution
      jackstat <- matrix(0, n, nstat)
      colnames(jackstat) <- varnames
      
      # create jackknife distribution
      if(parallel){
        
        jackfun <- function(i){
          index <- (1:n)[-i]
          jackdata <- object$data[index,]
          temp.object <- gsm(formula = object$formula, family = object$family,
                             data = jackdata, types = as.list(object$types), 
                             tprk = object$specs$tprk, knots = object$specs$knots, 
                             skip.iter = object$specs$skip.iter, lambda = lambda, 
                             control = object$specs$control, method = object$method, 
                             xrange = object$specs$xrng, thetas = thetas, mf = jackdata)
          if(nostat) temp.object$data <- object$data
          statistic(temp.object, ...)
        } # end jackfun
        jackstat <- t(parallel::parSapply(cl = cl, X = 1:n, FUN = jackfun))
        
      } else {
        
        for(i in 1:n){
          index <- (1:n)[-i]
          jackdata <- object$data[index,]
          temp.object <- gsm(formula = object$formula, family = object$family,
                             data = jackdata, types = as.list(object$types), 
                             tprk = object$specs$tprk, knots = object$specs$knots, 
                             skip.iter = object$specs$skip.iter, lambda = lambda, 
                             control = object$specs$control, method = object$method, 
                             xrange = object$specs$xrng, thetas = thetas, mf = jackdata)
          if(nostat) temp.object$data <- object$data
          jackstat[i,] <- statistic(temp.object, ...)
          if (verbose) setTxtProgressBar(pbar, nboot + i)
        } # end for(i in 1:n)
        
      } # end if(parallel)
      
      # get quantiles
      z1 <- qnorm(alpha/2)
      z2 <- qnorm(1 - alpha/2)
      
      # bias correction and acceleration
      confint <- matrix(0, nstat, 2)
      meanjack <- colMeans(jackstat)
      z0 <- acc <- rep(0, nstat)
      for(j in 1:nstat){
        z0[j] <- qnorm(mean(bootdist[,j] < t0[j]))
        acc[j] <- sum((meanjack[j] - jackstat[,j])^3) / (6 * sum((meanjack[j] - jackstat[,j])^2)^(3/2))
        alphas1 <- pnorm(z0[j] + (z0[j] + z1) / (1 - acc[j]*(z0[j] + z1)))
        alphas2 <- pnorm(z0[j] + (z0[j] + z2) / (1 - acc[j]*(z0[j] + z2)))
        probs <- c(alphas1, alphas2)
        confint[j,] <- quantile(bootdist[,j], probs = probs)
      }
      colnames(confint) <- c("lwr", "upr")
      attr(confint, "level") <- level
      
    } else {
      
      z0 <- acc <- rep(0, nstat)
      confint <- t(apply(bootdist, 2, quantile, probs = c(alpha/2, 1 - alpha/2)))
      colnames(confint) <- c("lwr", "upr")
      attr(confint, "level") <- level
      
    } # end if(bca)
    
    #########***#########   RETURN RESULTS   #########***#########
    
    # close progress bar or cluster
    if(parallel && make.cl) parallel::stopCluster(cl)
    if(verbose && !parallel) close(pbar)
    
    # create list
    res <- list(t0 = t0,
                se = apply(X = bootdist, MARGIN = 2, FUN = sd),
                bias = apply(X = bootdist, MARGIN = 2, FUN = mean) - t0,
                cov = if(cov.mat) cov(bootdist) else NULL,
                ci = confint, 
                boot.dist = if(boot.dist) bootdist else NULL,
                object = object, R = R, level = level, bca = bca, 
                method = method, fix.lambda = fix.lambda, 
                fix.thetas = fix.thetas, bias.correct = z0, acceleration = acc)
    
    # class and return
    class(res) <- "boot.gsm"
    return(res)
    
  } # end boot.gsm

print.boot.gsm <-
  function(x, header = FALSE, ...){
    nstat <- length(x$t0)
    nplab <- ifelse(x$method == "param", "Parametric", "Nonparametric")
    cilab <- ifelse(x$bca, "BCa", "Percentile")
    cat(c("\n", paste0(nplab, " Bootstrap of Generalized Smooth Model\n")), 
        sep = "")
    cat("using R = ", x$R, " bootstrap replicates\n", sep = "")
    cat("\n    method: ", x$method, sep = "")
    cat("\n     level: ", x$level, sep = "")
    cat("\n   confint: ", cilab, sep = "")
    cat("\nfix.lambda: ", x$fix.lambda, sep = "")
    cat("\nfix.thetas: ", x$fix.lambda, "\n", sep = "")
    if(header){
      cat(c("\n", paste0(round(x$lev * 100, 2), "% Confidence Interval\n")))
      print(head(x$ci, ...))
    }
    cat("\n")
  } # end print.boot.gsm

# simulate from exp fam dist
rexpfam <- function(object){
  n <- nrow(object$data)
  mu <- fitted(object)
  wts <- object$weights
  if(is.null(wts)) wts <- rep(1, n)
  if(object$family$family == "gaussian"){
    return(rnorm(n = n, mean = mu, sd = sqrt(object$dispersion / wts)))
  } else if(object$family$family == "binomial"){
    return(rbinom(n = n, size = wts, prob = mu) / wts)
  } else if(object$family$family == "poisson"){
    if(any(wts != 1)) warning("ignoring prior weights")
    return(rpois(n = n, lambda = mu))
  } else if(object$family$family == "NegBin"){
    if(any(wts != 1)) warning("ignoring prior weights")
    return(rnbinom(n = n, size = object$family$theta, mu = mu))
  } else if(object$family$family == "Gamma"){
    if(any(wts != 1)) message("using weights as shape parameters")
    shape <- wts / object$dispersion
    return(rgamma(n = n, shape = shape, rate = shape / mu))
  } else if(object$family$family == "inverse.gaussian"){
    if (!requireNamespace("SuppDists", quietly = TRUE)) 
      stop("need the 'SuppDists' package to simulate from 'inverse.gaussian' family")
    if (any(wts != 1)) message("using weights as inverse variances")
    return(SuppDists::rinvGauss(n = n, nu = mu, lambda = wts / object$dispersion))
  }
}
