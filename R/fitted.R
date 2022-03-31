# gsm objects
fitted.values.gsm <- function(object, ...){
  ginv <- object$family$linkinv
  return(ginv(object$linear.predictors))
}
fitted.gsm <- function(object, ...){
  ginv <- object$family$linkinv
  return(ginv(object$linear.predictors))
}

# sm objects
fitted.values.sm <- function(object, ...){
  return(object$fitted.values)
}
fitted.sm <- function(object, ...){
  return(object$fitted.values)
}


# ss objects
fitted.values.ss <- function(object, ...){
  if(is.null(object$data)){
    warning("Input 'object' has no data. Returning object$y as fitted values.")
    return(object$y)
  } else {
    return(predict.ss(object, x = object$data$x, se.fit = FALSE)$y)
  }
}
fitted.ss <- function(object, ...){
  if(is.null(object$data)){
    warning("Input 'object' has no data. Returning object$y as fitted values.")
    return(object$y)
  } else {
    return(predict.ss(object, x = object$data$x, se.fit = FALSE)$y)
  }
}
