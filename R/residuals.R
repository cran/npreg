# gsm
residuals.gsm <- 
  function(object, type = c("deviance", "pearson", "working",
                            "response", "partial"), ...){
    types <- c("deviance", "pearson", "working", "response", "partial")
    type <- pmatch(as.character(type[1]), types)
    if(is.na(type)) stop("Input 'type' must be one of the following options:\n'deviance', 'pearson', 'working', 'response', or 'partial'")
    type <- types[type]
    r <- object$residuals
    y <- object$data[,1]
    mu <- fitted(object)
    wt <- weights(object)
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
residuals.sm <- 
  function(object, type = c("working", "response", "deviance", 
                            "pearson", "partial"), ...){
    types <- c("working", "response", "deviance", "pearson", "partial")
    type <- pmatch(as.character(type[1]), types)
    if(is.na(type)) stop("Input 'type' must be one of the following options:\n'working', 'response', 'deviance', 'pearson', or 'partial'")
    type <- types[type]
    r <- object$data[, 1] - object$fitted.values
    if(any(type == c("deviance", "pearson")) && !is.null(object$weights)) 
      r <- r * sqrt(weights(object))
    if(type == "partial")
      r <- r + predict(object, type = "terms")
    return(r)
  }

# ss
residuals.ss <- 
  function(object, type = c("working", "response", "deviance", 
                            "pearson", "partial"), ...){
    if(is.null(object$data)){
      stop("Input 'object' has no data, which is needed to calculate residuals.")
    }
    types <- c("working", "response", "deviance", "pearson", "partial")
    type <- pmatch(as.character(type[1]), types)
    if(is.na(type)) stop("Input 'type' must be one of the following options:\n'working', 'response', 'deviance', 'pearson', or 'partial'")
    type <- types[type]
    r <- object$data$y - fitted(object)
    if(any(type == c("deviance", "pearson"))) 
      r <- r * sqrt(weights(object))
    if(type == "partial"){
      pt <- as.matrix(predict(object, x = object$data$x)$y - object$fit$coef[1])
      colnames(pt) <- "x"
      attr(pt, "constant") <- object$fit$coef[1]
      r <- r + pt
    }
    return(r)
  }