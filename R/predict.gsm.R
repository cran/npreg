predict.gsm <- 
  function(object, newdata = NULL, se.fit = FALSE, 
           type = c("link", "response", "terms"), 
           terms = NULL, na.action = na.pass,
           intercept = NULL, combine = TRUE, design = FALSE, 
           check.newdata = TRUE, ...){
    # predict method for class "gsm"
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2019-04-27
    
    
    #########***#########   CHECKS   #########***#########
    
    ### check if newdata provided
    provided.newdata <- TRUE
    if(is.null(newdata)) {
      newdata <- object$data[,-1,drop=FALSE]
      provided.newdata <- FALSE
    }
    
    ### check se.fit
    se.fit <- as.logical(se.fit[1])
    if(!any(se.fit == c(TRUE, FALSE))) stop("Input 'se.fit' must be a logical (TRUE/FALSE).")
    
    ### check type
    topts <- c("link", "response", "terms")
    type <- pmatch(as.character(type[1]), topts)
    if(is.na(type)) stop("Input 'type' must be 'response' or 'terms'.")
    type <- topts[type]
    
    ### check terms
    if(is.null(terms)){
      all.terms <- TRUE
      terms <- object$terms
    } else {
      all.terms <- FALSE
      tidx <- match(terms, object$terms)
      if(any(is.na(tidx))) stop("Input 'terms' contains terms that are not included in the model.")
      terms <- object$terms[sort(tidx)]
      all.terms <- ifelse(identical(terms, object$terms), TRUE, FALSE)
    }
    nterms <- length(terms)
    
    ### check na.action
    if(!is.function(na.action)) stop("Input 'na.action' must be a function determining what should be done with missing values in newdata.")
    
    ### check intercept
    if(missing(intercept) || is.null(intercept)){
      intercept <- ifelse(type != "terms" && all.terms, TRUE, FALSE)
    } else {
      intercept <- as.logical(intercept[1])
      if(!any(intercept == c(TRUE, FALSE))) stop("Input 'intercept' must be NULL or a logical (TRUE/FALSE).")
    }
    
    ### check combine
    combine <- as.logical(combine[1])
    if(!any(combine == c(TRUE, FALSE))) stop("Input 'combine' must be a logical (TRUE/FALSE).")
    
    ### check design
    design <- as.logical(design[1])
    if(!any(design == c(TRUE, FALSE))) stop("Input 'design' must be a logical (TRUE/FALSE).")
    
    ### check check.newdata
    check.newdata <- as.logical(check.newdata[1])
    if(!any(check.newdata == c(TRUE, FALSE))) stop("Input 'check.newdata' must be a logical (TRUE/FALSE).")
    
    
    #########***#########   RETURN TERMS   #########***#########
    
    ### check type is terms
    if(type == "terms"){
      
      # initializations
      fit <- vector("list", nterms)
      names(fit) <- terms
      if(se.fit) se <- vector("list", nterms)
      if(design) {
        X <- vector("list", nterms)
        names(X) <- terms
      }
      if(!provided.newdata) newdata <- NULL
      
      # loop through terms
      for(t in 1:nterms){
        tfit <- predict(object, newdata = newdata, se.fit = se.fit, 
                        type = "link", terms = terms[t], na.action = na.action,
                        intercept = intercept, combine = combine, 
                        design = design, check.newdata = check.newdata)
        if(se.fit | design){
          fit[[t]] <- tfit$fit
          if(se.fit) se[[t]] <- tfit$se.fit
          if(design) X[[t]] <- tfit$X
        } else {
          fit[[t]] <- tfit
        } 
      } # end for(t in 1:nterms)
      
      # collect results
      if(combine){
        fit <- do.call("cbind", fit)
        colnames(fit) <- terms
        if(se.fit) {
          se <- do.call("cbind", se)
          colnames(se) <- terms
        }
      } else {
        fit <- list(p = sapply(fit, function(x) x$p),
                    s = sapply(fit, function(x) x$s))
        colnames(fit$p) <- colnames(fit$s) <- terms
        if(se.fit){
          se <- list(p = sapply(se, function(x) x$p),
                     s = sapply(se, function(x) x$s))
          colnames(se$p) <- colnames(se$s) <- terms
        }
      }
      
      # return terms results
      attr(fit, "constant") <- object$coefficients[1]
      if(se.fit | design) {
        res <- list(linear.predictors = fit)
        if(se.fit) res$se.fit <- se
        if(design) res$X <- X
        return(res)
      } else {
        return(fit)
      }
      
    } # end if(return.type == "terms")
    
    
    #########***#########   MODEL FRAME   #########***#########
    
    ### check formula
    spform <- strsplit(as.character(object$formula), split = " ~ ")
    if(length(spform) == 3L && spform[[3]] == "."){
      object$formula <- as.formula(paste0(spform[[2]], " ~ ", paste0(object$terms, collapse = " + ")))
    }
    
    ### get effect table
    etab <- attr(terms(object$formula), "factors")
    
    ### zero unwanted terms (if applicable)
    if(!identical(terms, object$terms)){
      rmterms <- !(object$terms %in% terms)
      etab[,rmterms] <- rep(0, nrow(etab))
      etab <- etab[, colSums(etab) > 0, drop=FALSE]
      etab <- etab[c(1, which(rowSums(etab) > 0)), ,drop=FALSE]
    }
    varnames <- names(which(rowSums(etab) > 0))
    
    ### does newdata have needed variables?
    nidx <- match(varnames, names(newdata))
    if(provided.newdata && any(is.na(nidx))) stop("Input 'newdata' is missing variables needed\n to form some of the requested 'terms'.")
    newdata <- newdata[nidx]
    
    ### formula and model frame
    form <- as.formula(paste("y  ~", paste(terms, collapse = " + ")))
    mf <- model.frame(form, data = cbind(y = 1, newdata), na.action = na.action)
    
    
    #########***#########   MODEL MATRIX   #########***#########
    
    ### get fitting info for newdata
    dfun <- function(x) ncol(as.matrix(x))
    tidx <- match(varnames, names(object$types))
    types <- object$types[tidx]         # kernel types
    knots <- object$specs$knots[tidx]   # spline knots
    xdim <- sapply(knots, dfun)         # get dimension
    xrng <- object$specs$xrng[tidx]     # data ranges
    xlev <- object$specs$xlev[tidx]     # data levels
    
    ### check new data (if provided)
    if(provided.newdata && check.newdata){
      chty <- check_type2(mf = mf, type = as.list(types),
                          xdim = xdim, xrng = xrng, xlev = xlev)
    } # end if(provided.newdata)
    
    ### get reproducing kernel functions
    rkhs <- pred_rkhs(x = as.list(mf[,-1,drop=FALSE]), type = types, 
                      knots = knots, xrng = xrng, xlev = xlev)
    
    ### get theta names
    newfn <- colnames(etab)
    thetanames <- gsub(".c", "", names(object$specs$thetas), fixed = TRUE)
    thetanames <- gsub(".n", "", thetanames, fixed = TRUE)
    thetas <- NULL
    for(k in 1:length(newfn)){
      ix <- which(thetanames == newfn[k])
      thetas <- c(thetas, object$specs$thetas[ix])
    }
    
    ### build design matrix
    depe <- pred_depe(Etab = etab, rkhs = rkhs, 
                      tprk = object$specs$tprk, thetas = thetas)
    
    ### include intercept?
    if(!intercept) {
      if(ncol(depe$K) > 1L){
        depe$K <- depe$K[,-1,drop=FALSE]
      } else {
        depe$K <- NULL
      }
    }
    
    
    #########***#########   PREDICTIONS   #########***#########
    
    ### fitted values for null space
    Kmatch <- NULL
    Kfit <- 0
    if(!is.null(depe$K)){
      Kold <- names(object$coef[1:object$nsdf])
      Knew <- colnames(depe$K)
      Kmatch <- match(Knew, Kold)
      Kfit <- depe$K %*% object$coef[Kmatch]
    }
    
    ### fitted values for contrast space
    Jmatch <- NULL
    Jfit <- 0
    if(!is.null(depe$J)){
      Jold <- names(object$coef[(1+object$nsdf):length(object$coef)])
      Jnew <- colnames(depe$J)
      Jmatch <- match(Jnew, Jold) + object$nsdf
      Jfit <- depe$J %*% object$coef[Jmatch]
    }
    
    ### combine null and contrast fits?
    if(combine){
      fit <- as.numeric(Kfit + Jfit)
      if(se.fit){
        coefid <- c(Kmatch, Jmatch)
        se <- sqrt(rowSums( (cbind(depe$K, depe$J) %*% object$cov.sqrt[coefid,])^2 ))
      }
      if(design) X <- cbind(depe$K, depe$J)
    } else {
      fit <- data.frame(p = as.numeric(Kfit), s = as.numeric(Jfit))
      if(se.fit){
        se.p <- 0
        if(!is.null(Kmatch)) se.p <- sqrt(rowSums( (depe$K %*% object$cov.sqrt[Kmatch,])^2 ))
        if(!is.null(Jmatch)) se.s <- sqrt(rowSums( (depe$J %*% object$cov.sqrt[Jmatch,])^2 ))
        se <- data.frame(p = se.p, s = se.s)
      }
      if(design) X <- list(p = depe$K, s = depe$J)
    }
    
    ### return results
    if(se.fit | design){
      if(type == "link"){
        res <- list(fit = fit)
        if(se.fit) res$se.fit <- se
      } else {
        res <- list(fit = object$family$linkinv(fit))
        if(se.fit) res$se.fit <- sqrt((se^2) * (linkinvd(object$family$link)(fit))^2)
      }
      if(design) res$X <- X
      return(res)
    } else {
      if(type == "link"){
        return(fit)
      } else {
        return(object$family$linkinv(fit))
      }
    }
    
  } # end predict.smnet