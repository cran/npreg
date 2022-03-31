# gsm
residuals.gsm <- function(object, type = c("deviance", "pearson", "working",
                                           "response", "partial"), ...){
  types <- c("deviance", "pearson", "working", "response", "partial")
  type <- pmatch(as.character(type[1]), types)
  if(is.na(type)) stop("Input 'type' must be one of the following options:\n'deviance', 'pearson', 'working', 'response', or 'partial'")
  type <- types[type]
  r <- object$residuals
  y <- object$data[,1]
  mu <- fitted(object)
  wt <- object$weights
  if(is.null(wt)) wt <- rep(1, nrow(object$data))
  if(type == "deviance"){
    d.res <- sqrt(pmax(object$family$dev.resids(y, mu, wt), 0))
    return( ifelse(y > mu, d.res, -d.res) )
  } else if(type == "pearson"){
    return( (y - mu) * sqrt(wt) / sqrt(object$family$variance(mu)) )
  } else if(type == "working") {
    return(r)
  } else if(type == "response"){
    return(y - mu)
  } else if(type == "partial"){
    return(r + predict(object, type = "terms"))
  }
}
resid.gsm <- function(object, type = c("deviance", "pearson", "working",
                                       "response", "partial"), ...){
  types <- c("deviance", "pearson", "working", "response", "partial")
  type <- pmatch(as.character(type[1]), types)
  if(is.na(type)) stop("Input 'type' must be one of the following options:\n'deviance', 'pearson', 'working', 'response', or 'partial'")
  type <- types[type]
  r <- object$residuals
  y <- object$data[,1]
  mu <- fitted(object)
  wt <- object$weights
  if(is.null(wt)) wt <- rep(1, nrow(object$data))
  if(type == "deviance"){
    d.res <- sqrt(pmax(object$family$dev.resids(y, mu, wt), 0))
    return( ifelse(y > mu, d.res, -d.res) )
  } else if(type == "pearson"){
    return( (y - mu) * sqrt(wt) / sqrt(object$family$variance(mu)) )
  } else if(type == "working") {
    return(r)
  } else if(type == "response"){
    return(y - mu)
  } else if(type == "partial"){
    return(r + predict(object, type = "terms"))
  }
}


# sm
residuals.sm <- function(object, ...){
  return(object$data[, 1] - object$fitted.values)
}
resid.sm <- function(object, ...){
  return(object$data[, 1] - object$fitted.values)
}

# ss
residuals.ss <- function(object, ...){
  if(is.null(object$data)){
    stop("Input 'object' has no data, which is needed to calculate residuals.")
  }
  return(object$data$y - fitted(object))
}
resid.ss <- function(object, ...){
  if(is.null(object$data)){
    stop("Input 'object' has no data, which is needed to calculate residuals.")
  }
  return(object$data$y - fitted(object))
}
