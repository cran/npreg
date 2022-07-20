weights.ss <- 
  function(object, ...){
    if(is.null(object$data)) {
      warning("Input 'object' has no data. Returning object$w as weights")
      w <- object$w
    } else {
      w <- object$data$w
      if(is.null(w)) w <- rep(1, nrow(object$data))
    }
    return(w)
  }

weights.sm <- weights.gsm <-
  function(object, ...){
    w <- object$weights
    if(is.null(w)) w <- rep(1, nrow(object$data))
    return(w)
  }