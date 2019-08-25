sm <- 
  function(formula, data, weights, types = NULL, tprk = TRUE, knots = NULL, 
           update = TRUE, spar = NULL, lambda = NULL, control = list(),
           method = c("GCV", "OCV", "GACV", "ACV", "REML", "ML", "AIC", "BIC")){
    # semiparametric model
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2019-08-20
    
    
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
    
    # check 'method'
    methods <- c("GCV", "OCV", "GACV", "ACV", "REML", "ML", "AIC", "BIC")
    method <- as.character(method[1])
    method <- pmatch(toupper(method), methods)
    if(is.na(method)) stop("Invalid 'method' input.")
    method <- methods[method]
    
    #########***#########   FORMULA   #########***#########
    
    # create model frame
    mf <- match.call()
    m <- match(c("formula","data","weights"),names(mf),0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
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
    
    
    #########***#########   GET TYPES   #########***#########
    chty <- check_type(mf, types)
    
    
    #########***#########   GET KNOTS   #########***#########
    knot <- check_knot(mf = mf, type = chty$type, xrng = chty$xrng, 
                       xlev = chty$xlev, tprk = tprk, knots = knots)
    
    
    #########***#########   GET RKHS   #########***#########
    xlist <- as.list(mf[,-1,drop=FALSE])
    if(any(names(xlist) == "(weights)")) {
      xlist <- xlist[-which(names(xlist) == "(weights)")]
    }
    rkhs <- build_rkhs(x = xlist, type = chty$type, knots = knot, xrng = chty$xrng)
    
    
    #########***#########   DESIGN & PENALTY   #########***#########
    depe <- build_depe(Etab = et, rkhs = rkhs, tprk = tprk)
    if(missing(weights)){
      no.weights <- TRUE
      depe$weights <- rep(1, nobs)
    } else {
      no.weights <- FALSE
      depe$weights <- as.numeric(mf[,"(weights)"])
      if(any(depe$weights < 0)) stop("Input 'weights' must be non-negative.")
      depe$weights <- depe$weights * sum(depe$weights > 0) / sum(depe$weights)
    }
    
    
    #########***#########   FIT MODEL   #########***#########
    
    # fit model
    fit <- fit_sm(y = yvar, depe = depe, lambda = lambda,
                  tprk = tprk, control = control, method = method)
    
    # use Algorithm 3.2 from Gu and Wahba (1991)
    nullindx <- 1:fit$nsdf
    allzero <- (sum(abs(fit$coefficients[-nullindx])) == 0)
    if(update && (length(depe$thetas) > 1L) && !allzero){
      
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
      if(missing(weights)){
        depe$weights <- rep(1, nobs)
      } else {
        depe$weights <- as.numeric(mf[,"(weights)"])
      }
      
      # refit model
      fit <- fit_sm(y = yvar, depe = depe, lambda = lambda,
                    tprk = tprk, control = control, method = method)
      
    } # end if(update && (length(depe$thetas) > 1L) && !allzero)
    
    
    #########***#########   RETURN RESULTS   #########***#########
    
    # fit specs
    specs <- list(knots = knot, thetas = depe$thetas, 
                  xrng = chty$xrng, xlev = chty$xlev,
                  tprk = tprk)
    
    # collect results
    fit$specs <- specs
    fit$data <- mf
    fit$types <- chty$type
    fit$terms <- colnames(et)
    fit$method <- if(is.null(lambda)) method
    fit$formula <- formula
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
    if(is.null(x$method) || x$method == "GCV"){
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
