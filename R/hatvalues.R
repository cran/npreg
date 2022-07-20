hatvalues.ss <- 
  function(model, ...){
    if(is.null(model$data)) {
      warning("Input 'object' has no data. Returning object$lev as leverages.")
      hat <- model$lev
    } else {
      hat <- weights(model) * (predict(model, x = model$data$x)$se / model$sigma)^2
    }
    hat[hat > 1 - 10 * .Machine$double.eps] <- 1
    return(hat)
  }

hatvalues.sm <- 
  function(model, ...){
    hat <- weights(model) * (model$se.fit / model$sigma)^2 
    hat[hat > 1 - 10 * .Machine$double.eps] <- 1
    return(hat)
  }

hatvalues.gsm <- 
  function(model, ...){
    hat <- weights(model) * model$se.lp^2 / model$dispersion
    hat[hat > 1 - 10 * .Machine$double.eps] <- 1
    return(hat)
  }