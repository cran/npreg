boot.ss <-
  function(object, statistic, ..., 
           R = 9999, level = 0.95, bca = TRUE, 
           method = c("cases", "resid", "param"), 
           fix.lambda = TRUE, cov.mat = FALSE, 
           boot.dist = FALSE, verbose = TRUE,
           parallel = FALSE, cl = NULL){
    # Bootstrap a fit Smoothing Spline
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Update: 2021-10-29
    
    
    
    #########***#########   INITIAL CHECKS   #########***#########
    
    # check 'object'
    if(class(object) != "ss") stop("Input 'object' must be of class 'ss'")
    if(is.null(object$data)) stop("Input 'object' has no data: object$data is NULL\n Need to refit model with 'keep.data = TRUE'")
    n <- nrow(object$data)
    xmax <- object$fit$min + object$fit$range
    
    # check 'statistic'
    if(missing(statistic)){
      statistic <- function(object, ...){
        xseq <- seq(object$fit$min, object$fit$min + object$fit$range, length.out = 201)
        predict(object, x = xseq)$y
      }
      x0 <- seq(object$fit$min, object$fit$min + object$fit$range, length.out = 201)
    } else {
      statistic <- as.function(statistic)
      x0 <- NA
    }
    t0 <- statistic(object, ...)
    if(!is.vector(t0))
      stop("Output of the 'statistic' function must be a vector")
    nstat <- length(t0)
    varnames <- names(t0)
    if(is.na(x0[1])) x0 <- 1:nstat
    
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
          temp.object <- suppressWarnings(
            {
              ss(x = object$data$x[index], 
                 y = object$data$y[index],
                 w = object$data$w[index], 
                 lambda = lambda,
                 method = object$method,
                 m = object$fit$m,
                 periodic = object$fit$periodic,
                 knots = object$fit$knot,
                 df.offset = object$fit$df.offset,
                 penalty = object$fit$penalty,
                 control.spar = object$fit$control.spar,
                 tol = object$tol,
                 bernoulli = object$fit$bernoulli,
                 xmin = object$fit$min,
                 xmax = xmax)
            }
          )
          statistic(temp.object, ...)
        } # end parfun
        bootdist <- t(cbind(t0, parallel::parSapply(cl = cl, X = 2:nboot, FUN = parfun)))
        
      } else {
        
        for(i in 2:nboot){
          
          index <- sample.int(n, replace = TRUE)
          temp.object <- suppressWarnings(
            {
              ss(x = object$data$x[index], 
                 y = object$data$y[index],
                 w = object$data$w[index], 
                 lambda = lambda,
                 method = object$method,
                 m = object$fit$m,
                 periodic = object$fit$periodic,
                 knots = object$fit$knot,
                 df.offset = object$fit$df.offset,
                 penalty = object$fit$penalty,
                 control.spar = object$fit$control.spar,
                 tol = object$tol,
                 bernoulli = object$fit$bernoulli,
                 xmin = object$fit$min,
                 xmax = xmax)
            }
          )
          bootdist[i,] <- statistic(temp.object, ...)
          if (verbose) setTxtProgressBar(pbar, i)
          
        } # end for(i in 2:nboot)
        
      } # end if(parallel)
      
    } else if(method == "resid") {
      
      yhat <- predict(object, x = object$data$x)$y
      ehat <- object$data$y - yhat
      
      if(parallel){
        
        parfun <- function(x){
          index <- sample.int(n, replace = TRUE)
          temp.object <- suppressWarnings(
            {
              ss(x = object$data$x, 
                 y = yhat + ehat[index],
                 w = object$data$w, 
                 lambda = lambda,
                 method = object$method,
                 m = object$fit$m,
                 periodic = object$fit$periodic,
                 knots = object$fit$knot,
                 df.offset = object$fit$df.offset,
                 penalty = object$fit$penalty,
                 control.spar = object$fit$control.spar,
                 tol = object$tol,
                 bernoulli = object$fit$bernoulli,
                 xmin = object$fit$min,
                 xmax = xmax)
            }
          )
          statistic(temp.object, ...)
        } # end parfun
        bootdist <- t(cbind(t0, parallel::parSapply(cl = cl, X = 2:nboot, FUN = parfun)))
        
      } else {
        
        for(i in 2:nboot){
          
          index <- sample.int(n, replace = TRUE)
          temp.object <- suppressWarnings(
            {
              ss(x = object$data$x, 
                 y = yhat + ehat[index],
                 w = object$data$w, 
                 lambda = lambda,
                 method = object$method,
                 m = object$fit$m,
                 periodic = object$fit$periodic,
                 knots = object$fit$knot,
                 df.offset = object$fit$df.offset,
                 penalty = object$fit$penalty,
                 control.spar = object$fit$control.spar,
                 tol = object$tol,
                 bernoulli = object$fit$bernoulli,
                 xmin = object$fit$min,
                 xmax = xmax)
            }
          )
          bootdist[i,] <- statistic(temp.object, ...)
          if (verbose) setTxtProgressBar(pbar, i)
        
      } # end if(parallel)
        
      } # end for(i in 2:nboot)
      
    } else {
      
      yhat <- predict(object, x = object$data$x)$y
      
      if(parallel){
        
        parfun <- function(x){
          temp.object <- suppressWarnings(
            {
              ss(x = object$data$x, 
                 y = yhat + rnorm(n, sd = object$sigma),
                 w = object$data$w, 
                 lambda = lambda,
                 method = object$method,
                 m = object$fit$m,
                 periodic = object$fit$periodic,
                 knots = object$fit$knot,
                 df.offset = object$fit$df.offset,
                 penalty = object$fit$penalty,
                 control.spar = object$fit$control.spar,
                 tol = object$tol,
                 bernoulli = object$fit$bernoulli,
                 xmin = object$fit$min,
                 xmax = xmax)
            }
          )
          statistic(temp.object, ...)
        } # end parfun
        bootdist <- t(cbind(t0, parallel::parSapply(cl = cl, X = 2:nboot, FUN = parfun)))
        
      } else {
        
        for(i in 2:nboot){
          
          temp.object <- suppressWarnings(
            {
              ss(x = object$data$x, 
                 y = yhat + rnorm(n, sd = object$sigma),
                 w = object$data$w, 
                 lambda = lambda,
                 method = object$method,
                 m = object$fit$m,
                 periodic = object$fit$periodic,
                 knots = object$fit$knot,
                 df.offset = object$fit$df.offset,
                 penalty = object$fit$penalty,
                 control.spar = object$fit$control.spar,
                 tol = object$tol,
                 bernoulli = object$fit$bernoulli,
                 xmin = object$fit$min,
                 xmax = xmax)
            }
          )
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
          temp.object <- suppressWarnings(
            {
              ss(x = object$data$x[index], 
                 y = object$data$y[index],
                 w = object$data$w[index], 
                 lambda = lambda,
                 method = object$method,
                 m = object$fit$m,
                 periodic = object$fit$periodic,
                 knots = object$fit$knot,
                 df.offset = object$fit$df.offset,
                 penalty = object$fit$penalty,
                 control.spar = object$fit$control.spar,
                 tol = object$tol,
                 bernoulli = object$fit$bernoulli,
                 xmin = object$fit$min,
                 xmax = xmax)
            }
          )
          statistic(temp.object, ...)
        } # end jackfun
        jackstat <- t(parallel::parSapply(cl = cl, X = 1:n, FUN = jackfun))
        
      } else {
        
        for(i in 1:n){
          index <- (1:n)[-i]
          temp.object <- suppressWarnings(
            {
              ss(x = object$data$x[index], 
                 y = object$data$y[index],
                 w = object$data$w[index], 
                 lambda = lambda,
                 method = object$method,
                 m = object$fit$m,
                 periodic = object$fit$periodic,
                 knots = object$fit$knot,
                 df.offset = object$fit$df.offset,
                 penalty = object$fit$penalty,
                 control.spar = object$fit$control.spar,
                 tol = object$tol,
                 bernoulli = object$fit$bernoulli,
                 xmin = object$fit$min,
                 xmax = xmax)
            }
          )
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
      
      z0 <- acc <-rep(0, nstat)
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
                x0 = x0, bias.correct = z0, acceleration = acc)
    
    # class and return
    class(res) <- "boot.ss"
    return(res)
    
  } # end boot.ss

print.boot.ss <-
  function(x, header = FALSE, ...){
    nstat <- length(x$t0)
    nplab <- ifelse(x$method == "param", "Parametric", "Nonparametric")
    cilab <- ifelse(x$bca, "BCa", "Percentile")
    cat(c("\n", paste0(nplab, " Bootstrap of Smoothing Spline\n")), 
        sep = "")
    cat("using R = ", x$R, " bootstrap replicates\n", sep = "")
    cat("\n    method: ", x$method, sep = "")
    cat("\n     level: ", x$level, sep = "")
    cat("\n   confint: ", cilab, sep = "")
    cat("\nfix.lambda: ", x$fix.lambda, "\n", sep = "")
    if(header){
      cat(c("\n", paste0(round(x$lev * 100, 2), "% Confidence Interval\n")))
      print(head(x$ci, ...))
    }
    cat("\n")
  } # end print.boot.ss

plot.boot.ss <-
  function(x, n = 201, ci = TRUE, xseq = NULL, ...){
    if(is.null(xseq)) xseq <- x$x0
    plot(x$object, n = n, ci = FALSE, xseq = xseq, ...)
    if(ci) plotci(x$x0, x$t0, ci = x$ci, add = TRUE)
  }

