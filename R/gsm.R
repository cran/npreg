gsm <- 
  function(formula, family = gaussian, data, weights, types = NULL, tprk = TRUE, 
           knots = NULL, skip.iter = TRUE, spar = NULL, lambda = NULL, control = list(),
           method = c("GCV", "OCV", "GACV", "ACV", "PQL", "AIC", "BIC"),
           xrange = NULL, thetas = NULL, mf = NULL){
    # generalized smooth model
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2022-03-30
    
    
    #########***#########   CHECKS   #########***#########
    
    # check 'family'
    if(missing(family)) family <- gaussian
    family <- check_family(family)
    
    # check 'spar'
    if(!is.null(spar)){
      spar <- as.numeric(spar[1])
      lambda <- 256^(3*(spar-1))
    }
    
    # check 'lambda'
    if(!is.null(lambda)){
      lambda <- as.numeric(lambda[1])
      if(lambda < 0) stop("Input 'lambda' must be a non-negative scalar.")
    }
    
    # check 'control'
    control <- do.call("check_control", control)
    
    # check 'method'
    method <- as.character(method[1])
    method <- pmatch(toupper(method), c("GCV", "OCV", "GACV", "ACV", "PQL", "AIC", "BIC"))
    if(is.na(method)) stop("Invalid 'method' input.")
    method <- c("GCV", "OCV", "GACV", "ACV", "PQL","AIC", "BIC")[method]
    
    
    #########***#########   FORMULA   #########***#########
    
    # create model frame
    if(is.null(mf)){
      mf <- match.call()
      m <- match(c("formula","data","weights"),names(mf),0L)
      mf <- mf[c(1L, m)]
      mf[[1L]] <- as.name("model.frame")
      mf <- eval(mf, parent.frame())
    }
    mt <- attr(mf, "terms")                 # mt contains model info and terms 
    et <- attr(mt,"factors")                # et is effects table
    mfdim <- dim(et)                        # dim of effects table
    nobs <- dim(mf)[1]                      # total number of data points
    nxvar <- mfdim[1] - 1L                  # number of predictors
    nterm <- mfdim[2]                       # number of model terms
    xynames <- row.names(et)
    xnames <- xynames[2:(nxvar+1L)]
    tnames <- colnames(et)
    if(any(colSums(et>0L)>3L)){stop("Four-way (and higher-order) interactions are not supported.")}
    for(k in 1:nxvar){
      if(class(mf[,k+1])[1] == "character") mf[,k+1] <- factor(mf[,k+1])
    }
    
    # get response
    yvar <- model.response(mf)
    if(family$family == "binomial"){
      if(any(class(yvar)[1] == c("factor", "ordered"))){
        lev1 <- levels(yvar)[1]
        yvar <- ifelse(yvar == lev1, 0, 1)
      }
      yvar <- as.matrix(yvar + 0.0)
    } else if(family$family == "poisson"){
      yvar <- as.matrix(as.integer(yvar))
    } else {
      yvar <- as.matrix(yvar + 0.0)
    }
    if(ncol(yvar) > 1){ stop("Response must be unidimensional (vector).") }
    
    
    #########***#########   WEIGHTS   #########***#########
    xlist <- as.list(mf[,-1,drop=FALSE])
    if(any(names(xlist) == "(weights)")) {
      no.weights <- FALSE
      xlist <- xlist[-which(names(xlist) == "(weights)")]
      weights <- as.numeric(mf[,"(weights)"])
      if(any(weights < 0)) stop("Input 'weights' must be non-negative.")
      if(family$family != "binomial") weights <- weights * sum(weights > 0) / sum(weights)
    } else {
      no.weights <- TRUE
      weights <- rep(1, nobs)
    }
    
    
    #########***#########   GET TYPES   #########***#########
    chty <- check_type(mf, types, xrange)
    
    
    #########***#########   GET KNOTS   #########***#########
    knot <- check_knot(mf = mf, type = chty$type, xrng = chty$xrng, 
                       xlev = chty$xlev, tprk = tprk, knots = knots)
    
    
    #########***#########   GET RKHS   #########***#########
    rkhs <- build_rkhs(x = xlist, type = chty$type, knots = knot, 
                       xrng = chty$xrng, xlev = chty$xlev)
    
    
    #########***#########   DESIGN & PENALTY   #########***#########
    depe <- build_depe(Etab = et, rkhs = rkhs, tprk = tprk, thetas = thetas)
    depe$weights <- weights
    
    
    #########***#########   FIT MODEL   #########***#########
    
    # fit model
    fit <- fit_gsm(y = yvar, depe = depe, lambda = lambda,
                   tprk = tprk, control = control, method = method,
                   family = family)
    
    # use Algorithm 3.2 from Gu and Wahba (1991)
    nullindx <- 1:fit$nsdf
    allzero <- (sum(abs(fit$coefficients[-nullindx])) == 0)
    if((length(depe$thetas) > 1L) && is.null(thetas) && !allzero){
      
      # rebuild design and penalty
      if(tprk){
        depe <- build_depe2(Etab = et, rkhs = rkhs, 
                            thetas = depe$thetas,
                            Jcoef = fit$coef[-nullindx])
      } else {
        depe <- build_depe2(Etab = et, thetas = depe$thetas,
                            Jcoef = fit$coef[-nullindx],
                            depe = depe, tprk = FALSE)
      }
      depe$weights <- weights
      
      # refit model
      fit <- fit_gsm(y = yvar, depe = depe, lambda = lambda,
                     tprk = tprk, control = control, method = method,
                     family = family)
      
      # deep tuning?
      if(!skip.iter){
        opt <- nlm(f = tune.deep.gsm, p = log(depe$thetas),
                   spar = fit$spar, y = yvar, Etab = et, rkhs = rkhs, 
                   weights = depe$weights, tprk = tprk, method = method, 
                   family = fit$family, control = control, gradtol = control$tol, 
                   iterlim = control$iterlim, print.level = control$print.level)
        thetas <- exp(opt$estimate)
        names(thetas) <- names(depe$thetas)
        depe <- build_depe(Etab = et, rkhs = rkhs, tprk = tprk, thetas = thetas)
        depe$weights <- weights
        fit <- fit_gsm(y = yvar, depe = depe, lambda = fit$lambda,
                       tprk = tprk, control = control, method = method,
                       family = family)
        fit$iter <- c(irpls = fit$iter, nlm = opt$iter)
        if(family$family == "NegBin" && !family$fixed.theta){
          iter.out <- 0L
          vtol.out <- control$epsilon.out + 1
          obj.old <- as.numeric(opt$minimum)
          while(iter.out < control$maxit.out & vtol.out > control$epsilon.out){
            mu <- family$linkinv(cbind(depe$K, depe$J) %*% fit$coefficients)
            theta.hat <- theta.mle(y = yvar, mu = mu, wt = depe$weights)
            family <- NegBin(theta = theta.hat, link = family$link)
            opt <- nlm(f = tune.deep.gsm, p = log(depe$thetas),
                       spar = fit$spar, y = yvar, Etab = et, rkhs = rkhs,
                       weights = depe$weights, tprk = tprk, method = method,
                       family = family, control = control, 
                       gradtol = control$tol, iterlim = control$iterlim, 
                       print.level = control$print.level)
            thetas <- exp(opt$estimate)
            names(thetas) <- names(depe$thetas)
            vtol.out <- abs( (opt$minimum - obj.old) / (obj.old + 1e-4) )
            iter.out <- iter.out + 1L
            obj.old <- opt$minimum
            depe <- build_depe(Etab = et, rkhs = rkhs, tprk = tprk, thetas = thetas)
            depe$weights <- weights
            fit <- fit_gsm(y = yvar, depe = depe, lambda = fit$lambda,
                           tprk = tprk, control = control, method = method,
                           family = family)
            fit$iter <- c(inner = fit$iter, outter = attr(theta.hat, "iter"), nlm = opt$iter)
          } # end while(iter.out < control$maxit.out & vtol.out > control$epsilon.out)
          family$fixed.theta <- FALSE
          attr(fit$family$theta, "SE") <- attr(theta.hat, "SE")
          attr(fit$family$theta, "iter") <- attr(theta.hat, "iter")
        } # end if(family$family == "NegBin" && !family$fixed.theta)
      } # end if(!skip.iter)
      
    } # end if((length(depe$thetas) > 1L) && !allzero)
    
    
    #########***#########   RETURN RESULTS   #########***#########
    
    # check for NegBin with theta = NULL
    if(family$family == "NegBin" && !family$fixed.theta){
      family <- fit$family
      family$fixed.theta <- FALSE
      fit$family <- NULL
    } else {
      fit$family <- NULL
    }
    
    # fit specs
    specs <- list(knots = knot, thetas = depe$thetas, 
                  xrng = chty$xrng, xlev = chty$xlev,
                  tprk = tprk, skip.iter = skip.iter,
                  control = control)
    
    # collect results
    fit$specs <- specs
    fit$data <- mf
    fit$types <- chty$type
    fit$terms <- colnames(et)
    fit$method <- method
    fit$formula <- formula
    fit$weights <- if(!no.weights) weights
    fit$call <- match.call()
    fit$family <- family
    fit <- fit[c(4:length(fit), 1:3)]
    
    # transform GCV for Gaussian
    if(family$family == "gaussian" && (method %in% c("GCV", "OCV"))){
      fit$cv.crit <- fit$cv.crit * 2 + mean(yvar^2)
    }
    
    # class and return
    class(fit) <- "gsm"
    return(fit)
    
  } # end gsm

# print function
print.gsm <-
  function(x, ...){
    cat("\nCall:\n")
    print(x$call)
    #cat("\nFormula: ", as.character(as.expression(x$formula)))
    cat("\nSmoothing Parameter:  spar =", x$spar, "  lambda =", x$lambda)
    cat("\nEquivalent Degrees of Freedom (Df): ", x$df)
    weighted <- ifelse(any(colnames(x$data) == "(weights)"), TRUE, FALSE)
    if(weighted){
      cat("\nPenalized Criterion (weighted deviance): ", x$deviance)
    } else {
      cat("\nPenalized Criterion (deviance): ", x$deviance)
    }
    if(x$method == "GCV"){
      cat("\nGeneralized Cross-Validation (GCV): ", x$cv.crit,"\n\n")
    } else if(x$method == "OCV"){
      cat("\nOrdinary Cross-Validation (OCV): ", x$cv.crit,"\n\n")
    } else if(x$method == "GACV"){
      cat("\nGeneralized Approximate Cross-Validation (GACV): ", x$cv.crit,"\n\n")
    } else if(x$method == "ACV"){
      cat("\nApproximate Cross-Validation (ACV): ", x$cv.crit,"\n\n")
    } else if(x$method == "PQL"){
      cat("\nPenalized Quasi-Likelihood (PQL): ", x$cv.crit,"\n\n")
    } else if(x$method == "AIC"){
      cat("\nAkaike Information Criterion (AIC): ", x$aic,"\n\n")
    } else if(x$method == "BIC"){
      cat("\nBayesian Information Criterion (BIC): ", x$bic,"\n\n")
    }
  } # end print.gsm
