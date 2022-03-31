sm <- 
  function(formula, data, weights, types = NULL, tprk = TRUE, knots = NULL, 
           skip.iter = TRUE, df, spar = NULL, lambda = NULL, control = list(),
           method = c("GCV", "OCV", "GACV", "ACV", "REML", "ML", "AIC", "BIC"),
           xrange = NULL, thetas = NULL, mf = NULL){
    # smooth model
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2022-03-30
    
    
    #########***#########   CHECKS   #########***#########
    
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
    if(is.null(control$lower)){
      control$lower <- -1.5
    } else {
      control$lower <- as.numeric(control$lower[1])
    }
    if(is.null(control$upper)){
      control$upper <- 1.5
    } else {
      control$upper <- as.numeric(control$upper[1])
    }
    if(control$upper <= control$lower) stop("Search bounds must satisfy:  control$lower < control$upper")
    if(is.null(control$tol)){
      control$tol <- 1e-8
    } else {
      control$tol <- as.numeric(control$tol[1])
      if(control$tol <= 0) stop("Input 'control$tol' must be a positive scalar.")
    }
    if(is.null(control$iterlim)){
      control$iterlim <- 5000L
    } else {
      control$iterlim <- as.integer(control$iterlim[1])
      if(control$iterlim < 1L) stop("Input 'control$iterlim' must be a positive integer.")
    }
    if(is.null(control$print.level)){
      control$print.level <- 0L
    } else {
      control$print.level <- as.integer(control$print.level[1])
      if(!any(control$print.level == c(0L, 1L, 2L))) stop("Input 'control$print.level' must be either:\n0 (no printing), 1 (minimal printing), or 2 (full printing).")
    }
    
    # check 'method'
    methods <- c("GCV", "OCV", "GACV", "ACV", "REML", "ML", "AIC", "BIC")
    method <- as.character(method[1])
    method <- pmatch(toupper(method), methods)
    if(is.na(method)) stop("Invalid 'method' input.")
    method <- methods[method]
    
    
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
    yvar <- as.matrix(model.response(mf, "numeric")+0.0)   # response variable
    if(ncol(yvar)>1){stop("Response must be unidimensional (vector).")}
    for(k in 1:nxvar){
      if(class(mf[,k+1])[1] == "character") mf[,k+1] <- factor(mf[,k+1])
    }
    
    
    #########***#########   WEIGHTS   #########***#########
    xlist <- as.list(mf[,-1,drop=FALSE])
    if(any(names(xlist) == "(weights)")) {
      no.weights <- FALSE
      xlist <- xlist[-which(names(xlist) == "(weights)")]
      weights <- as.numeric(mf[,"(weights)"])
      if(any(weights < 0)) stop("Input 'weights' must be non-negative.")
      weights <- weights * sum(weights > 0) / sum(weights)
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
    
    # check df
    if(!missing(df)){
      df <- as.numeric(df[1])
      if(df < ncol(depe$K) | df > nobs) stop("'df' must satisfy:  m < df < n")
      if(!skip.iter) warning("Input 'df' will be inexact when 'skip.iter = FALSE'")
    } else {
      df <- NULL
    }
    
    # fit model
    fit <- fit_sm(y = yvar, depe = depe, df = df, lambda = lambda,
                  tprk = tprk, control = control, method = method)
    
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
      fit <- fit_sm(y = yvar, depe = depe, df = df, lambda = lambda,
                    tprk = tprk, control = control, method = method)
      
      # deep tuning?
      if(!skip.iter){
        opt <- nlm(f = tune.deep.sm, p = log(depe$thetas),
                   spar = fit$spar, y = yvar, Etab = et, rkhs = rkhs, 
                   weights = depe$weights, tprk = tprk, method = method, 
                   gradtol = control$tol, iterlim = control$iterlim,
                   print.level = control$print.level)
        thetas <- exp(opt$estimate)
        names(thetas) <- names(depe$thetas)
        depe <- build_depe(Etab = et, rkhs = rkhs, tprk = tprk, thetas = thetas)
        depe$weights <- weights
        fit <- fit_sm(y = yvar, depe = depe, df = df, lambda = fit$lambda,
                      tprk = tprk, control = control, method = method)
        fit$iter <- opt$iter
      }
      
    } # end if(update && (length(depe$thetas) > 1L) && !allzero)
    
    
    #########***#########   RETURN RESULTS   #########***#########
    
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
    
    # class and return
    class(fit) <- "sm"
    return(fit)
    
  } # end sm

# print function
print.sm <-
  function(x, ...){
    cat("\nCall:\n")
    print(x$call)
    #cat("\nFormula: ", as.character(as.expression(x$formula)))
    cat("\nSmoothing Parameter:  spar =", x$spar, "  lambda =", x$lambda)
    cat("\nEquivalent Degrees of Freedom (Df): ", x$df)
    weighted <- ifelse(any(colnames(x$data) == "(weights)"), TRUE, FALSE)
    if(weighted){
      cat("\nPenalized Criterion (weighted RSS): ", x$sse)
    } else {
      cat("\nPenalized Criterion (RSS): ", x$sse)
    }
    if(x$method == "GCV"){
      cat("\nGeneralized Cross-Validation (GCV): ", x$cv.crit,"\n\n")
    } else if(x$method == "OCV"){
      cat("\nOrdinary Cross-Validation (OCV): ", x$cv.crit,"\n\n")
    } else if(x$method == "GACV"){
      cat("\nGeneralized Approximate Cross-Validation (GACV)", x$cv.crit,"\n\n")
    } else if(x$method == "ACV"){
      cat("\nApproximate Cross-Validation (ACV)", x$cv.crit,"\n\n")
    } else if(x$method == "REML"){
      cat("\nLog Likelihood (REML): ", x$logLik,"\n\n")
    } else if(x$method == "ML"){
      cat("\nLog Likelihood (ML): ", x$logLik,"\n\n")
    } else if(x$method == "AIC"){
      cat("\nAkaike's Information Criterion (AIC)", x$aic,"\n\n")
    } else if(x$method == "BIC"){
      cat("\nBayesian Information Criterion (BIC)", x$bic,"\n\n")
    }
  } # end print.sm
