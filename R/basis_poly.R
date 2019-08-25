basis_poly <-
  function(x, knots, m = 2, d = 0, 
           xmin = min(x), xmax = max(x), 
           periodic = FALSE, rescale = FALSE,
           intercept = FALSE){
    # Polynomial Smoothing Spline Basis
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Update: 2019-04-04
    
    ### check x and k
    x <- as.numeric(x)
    knots <- as.numeric(knots)
    
    ### initializations
    m <- as.integer(m[1])
    if(m < 1L | m > 3L) stop("Input 'm' must be 1 (linear), 2 (cubic), or 3 (quintic).")
    x <- (x - xmin) / (xmax - xmin)
    knots <- (knots - xmin) / (xmax - xmin)
    
    ### check d
    d <- as.integer(d[1])
    if(d < 0 | d > 2) stop("Input 'd' must be 0, 1, or 2.")
    
    ### get kernel function name
    kerns <- c("lin", "cub", "qui")
    fname <- paste0(kerns[m], "kern", d)
    
    ### evaluate kernel function (or its derivative)
    Xs <- outer(X = x, Y = knots, FUN = fname, periodic = periodic)
    colnames(Xs) <- paste("knot", 1:length(knots), sep = ".")
    
    ### evaluate null space basis (or its derivative)
    if(periodic || m == 1L){
      if(intercept){
        Xp <- matrix(ifelse(d > 0, 0, 1), nrow = length(x), ncol = 1)
        colnames(Xp) <- "(Intercept)"
      } else {
        Xp <- NULL
      }
    } else {
      if(d == 0){
        Xp <- as.matrix(x - 1/2)
        if(m == 3L) Xp <- cbind(Xp, ((x - 1/2)^2 - (1/12))/2 )
      } else if(d == 1){
        Xp <- matrix(1, nrow = length(x), ncol = 1)
        if(m == 3L) Xp <- cbind(Xp, x - 1/2)
      } else {
        Xp <- matrix(0, nrow = length(x), ncol = 1)
        if(m == 3L) Xp <- cbind(Xp, 1)
      }
      nullnames <- "x"
      if(m == 3L) nullnames <- c(nullnames, "x^2")
      if(intercept) {
        Xp <- cbind(ifelse(d > 0, 0, 1), Xp)
        nullnames <- c("(Intercept)", nullnames)
      }
      colnames(Xp) <- nullnames
    } # end if(periodic || m == 1L)
    
    # rescale?
    if(rescale){
      fname <- paste0(kerns[m], "kern", 0)
      theta <- 0
      for(k in 1:length(knots)) {
        fargs <- list(x = knots[k], y = knots[k], periodic = periodic)
        theta <- theta + do.call(fname, fargs)
      }
      theta <- length(knots) / theta
    } else {
      theta <- 1
    }
    
    # return results
    return(cbind(Xp, theta * Xs))
    
  } # end basis_poly.R

# linear smoothing spline (deriv = 0)
linkern0 <- 
  function(x, y, periodic = FALSE){
    val <- ( (abs(x - y) - 1/2)^2 - (1/12) ) / 2
    if(!periodic) val <- val + (x - 1/2) * (y - 1/2)
    val
  }

# linear smoothing spline (deriv = 1)
linkern1 <- 
  function(x, y, periodic = FALSE){
    dif <- x - y
    sgn <- sign(dif)
    sgn[sgn == 0] <- 1
    val <-  sgn * (abs(dif) - 1/2)
    if(!periodic) val <- val + (y - 1/2)
    val
  }

# cubic smoothing spline (deriv = 0)
cubkern0 <- 
  function(x, y, periodic = FALSE){
    dif <- abs(x - y)
    val <- 0 - ((dif - 1/2)^4 - (((dif - 1/2)^2)/2) + 7/240) / 24
    if(!periodic) {
      vx <- ( (x - 1/2)^2 - 1/12 ) / 2
      vy <- ( (y - 1/2)^2 - 1/12 ) / 2
      val <- val + vx * vy
    }
    val
  }

# cubic smoothing spline (deriv = 1)
cubkern1 <- 
  function(x, y, periodic = FALSE){
    dif <- x - y
    adif <- abs(dif)
    sgn <- sign(dif)
    sgn[sgn == 0] <- 1
    val <- 0 - sgn * (4*((adif - 1/2)^3) - (adif - 1/2))/24
    if(!periodic) {
      vx <- x - 1/2
      vy <- ( (y - 1/2)^2 - 1/12 ) / 2
      val <- val + vx * vy
    }
    val
  }

# cubic smoothing spline (deriv = 2)
cubkern2 <- 
  function(x, y, periodic = FALSE){
    dif <- abs(x - y)
    val <- 0 - ((dif - 1/2)^2 - (1/12))/2 
    if(!periodic) {
      vx <- 1
      vy <- ( (y - 1/2)^2 - 1/12 ) / 2
      val <- val + vx * vy
    }
    val
  }

# quintic smoothing spline (deriv = 0)
quikern0 <- 
  function(x, y, periodic = FALSE){
    dif <- abs(x - y)
    val <- ( ((dif - 1/2)^6)/30 - ((dif - 1/2)^4)/24 + 7*((dif - 1/2)^2)/480 - 31/40320 ) / 24
    if(!periodic) {
      vx <- ( 4 * (x - 1/2)^3 - (x - 1/2) ) / 24
      vy <- ( 4 * (y - 1/2)^3 - (y - 1/2) ) / 24
      val <- val + vx * vy
    }
    val
  }

# quintic smoothing spline (deriv = 1)
quikern1 <- 
  function(x, y, periodic = FALSE){
    dif <- x - y
    adif <- abs(dif)
    sgn <- sign(dif)
    sgn[sgn == 0] <- 1
    val <- sgn * (((adif - 1/2)^5)/5 - ((adif - 1/2)^3)/6 + 7*(adif - 1/2)/240)/24
    if(!periodic) {
      vx <- ( (x - 1/2)^2 - (1 / 12) ) / 2
      vy <- (4 * (y - 1/2)^3 - (y - 1/2) ) / 24
      val <- val + vx * vy
    }
    val
  }

# quintic smoothing spline (deriv = 2)
quikern2 <- 
  function(x, y, periodic = FALSE){
    dif <- abs(x - y)
    val <- ((dif - 1/2)^4 - (((dif - 1/2)^2)/2) + 7/240)/24
    if(!periodic) {
      vx <- x - 1/2
      vy <- (4 * (y - 1/2)^3 - (y - 1/2) ) / 24
      val <- val + vx * vy 
    }
    val
  }
