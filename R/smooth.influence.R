smooth.influence <-
  function(model, do.coef = TRUE){
    # regression diagnostic basics for smooth models
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2022-07-19
    
    
    # check model
    clamod <- class(model)
    if(!any(clamod == c("ss", "sm", "gsm")))
      stop("Input 'model' must of be class 'ss', 'sm', or 'gsm'.")
    
    # check 'ss' class for data
    if(clamod == "ss" && is.null(model$data))
      stop("Input 'model$data' is NULL. You must refit model with 'keep.data = TRUE'.")
    
    # number of subjects and coefficients
    n <- nrow(model$data)
    m <- length(coef(model))
    
    # model response
    y <- if(clamod == "ss") model$data$y else model$data[,1]
    
    # weighted design and inverse crossproduct
    wsqrt <- sqrt(weights(model))
    if(clamod == "gsm"){
      mueta <- model$family$mu.eta(model$linear.predictors)
      varx <- model$family$variance(fitted(model))
      wsqrt <- wsqrt * sqrt(mueta^2 / varx)
    }
    M <- wsqrt * model.matrix(model)
    Ginv <- vcov(model) / ifelse(clamod == "gsm", model$dispersion, model$sigma^2)
    MGinv <- M %*% Ginv
    
    # gaussian model?
    isGaussian <- (clamod != "gsm") || (model$family$family == "gaussian")
    
    # weighted response, fitted values, residuals, and crossproducts
    r.w <- weighted.residuals(model)
    GMtMG <- crossprod(MGinv)
    if(isGaussian){
      y.w <- wsqrt * y
      f.w <- wsqrt * fitted(model)
      fMGinv <- crossprod(f.w, MGinv)
    } else {
      r.pearson <- residuals(model, type = "pearson")
    }
    
    # hat values
    hat <- hatvalues(model)
    
    # calculate coefficients
    if(do.coef) 
      coefficients <- if(isGaussian) {
        (r.w / (1 - hat)) * MGinv
      } else {
        (r.pearson / (1 - hat)) * MGinv
      }
    
    # squared norm of i-th column of smoothing matrix
    dfun0 <- function(x) as.numeric(t(x) %*% GMtMG %*% x)
    si.si <- apply(M, 1, dfun0)
    
    # calculate deviance
    if(isGaussian){
      dvec1 <- (1 + hat) * y.w^2 - (2 + hat) * 2 * y.w * f.w + (3 + hat) * f.w^2
      dvec2 <- r.w * rowSums(M * matrix(fMGinv, nrow = n, ncol = m, byrow = TRUE))
      dvec3 <- r.w^2 * (si.si - hat^2) / (1 - hat)
      delta <- dvec1 + 2 * dvec2 - dvec3
      deviance <- deviance(model) - delta / (1 - hat)
    } else {
      deviance <- deviance(model) - model$dispersion * r.pearson^2 * hat / (1 - hat)^2
    }
    
    # calculate df
    df <- model$df + (si.si - hat) / (1 - hat)
    
    # calculate sigma
    sigma <- sqrt(deviance / (n - df - 1))
    
    # return results
    res <- list(hat = hat, coefficients = if(do.coef) coefficients,
                deviance = deviance, df = df, sigma = sigma, wt.res = r.w)
    return(res)
    
  } # end function(model, do.coef = TRUE)

influence.ss <- influence.sm <-
  function(model, do.coef = TRUE, ...){
    smooth.influence(model, do.coef = do.coef, ...)
  }

influence.gsm <-
  function(model, do.coef = TRUE, ...){
    res <- smooth.influence(model, do.coef = do.coef, ...)
    pRes <- na.omit(residuals(model, type = "pearson"))[weights(model) != 0]
    #pRes <- naresid(model$na.action, pRes)
    names(res)[names(res) == "wt.res"] <- "dev.res"
    c(res, list(pear.res = pRes))
  }