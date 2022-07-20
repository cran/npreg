model.matrix.ss <-
  function(object, ...){
    #if(class(object) != "ss") stop("Input 'object' must be of class 'ss'")
    if(!inherits(object, "ss")) stop("Input 'object' must be of class 'ss'")
    if(is.null(object$data)) stop("object$data is NULL. Need to refit model with 'keep.data = TRUE'")
    return(basis.poly(x = object$data$x, knots = object$fit$knot, 
                      m = object$fit$m, xmin = object$fit$min, 
                      xmax = object$fit$min + object$fit$range, 
                      periodic = object$fit$periodic, intercept = TRUE, 
                      bernoulli = object$fit$bernoulli))
  }

model.matrix.sm <-
  function(object, ...){
    #if(class(object) != "sm") stop("Input 'object' must be of class 'sm'")
    if(!inherits(object, "sm")) stop("Input 'object' must be of class 'sm'")
    return(predict(object, design = TRUE)$X)
  }

model.matrix.gsm <-
  function(object, ...){
    #if(class(object) != "gsm") stop("Input 'object' must be of class 'gsm'")
    if(!inherits(object, "gsm")) stop("Input 'object' must be of class 'gsm'")
    return(predict(object, design = TRUE)$X)
  }