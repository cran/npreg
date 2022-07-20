rstandard.ss <- rstandard.sm <- 
  function(model, infl = NULL, sd = model$sigma, 
           type = c("sd.1", "predictive"), ...){
    types <- c("sd.1", "predictive")
    type <- pmatch(as.character(type[1]), types)
    if(is.na(type)) stop("Input 'type' must be 'sd.1' or 'predictive'")
    type <- types[type]
    wt.res <- weighted.residuals(model)
    hat <- hatvalues(model)
    if(type == "sd.1"){
      res <- wt.res / (sd * sqrt(1 - hat))
    } else {
      res <- wt.res / (1 - hat)
    }
    res[is.infinite(res)] <- NaN
    res
  }

rstandard.gsm <- 
  function(model, infl = NULL, type = c("deviance", "pearson"), ...){
    types <- c("deviance", "pearson")
    type <- pmatch(as.character(type[1]), types)
    if(is.na(type)) stop("Input 'type' must be 'deviance' or 'pearson'")
    type <- types[type]
    hat <- hatvalues(model)
    res <- residuals(model, type = type) / sqrt(model$dispersion * (1 - hat))
    res[is.infinite(res)] <- NaN
    res
  }

