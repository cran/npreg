check_family <-
  function(family, ...){
    
    # check family
    if (is.character(family)) 
      family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
      family <- family()
    if (is.null(family$family)) {
      print(family)
      stop("'family' not recognized")
    }
    
    # add logLik function to family
    if(family$family == "gaussian"){
      family$logLik <- function(y, n, mu, wt, dev){
        nobs <- length(y)
        0.5 * (sum(log(wt)) - nobs * (log(dev/nobs * 2 * pi) + 1))
      }
      family$theta <- function(mu){
        mu
      }
      family$cumulant <- function(mu){
        mu^2/2
      }
      family$canonical <- ifelse(family$link == "identity", TRUE, FALSE)
    } else if(family$family == "binomial"){
      family$logLik <- function(y, n, mu, wt, dev){
        m <- if (any(n > 1)) 
          n
        else wt
        ifelse(m > 0, (wt/m), 0) * dbinom(round(m * y), round(m), mu, log = TRUE)
      }
      family$theta <- function(mu){
        log(mu / (1 - mu))
      }
      family$cumulant <- function(mu){
        -log(1 - mu)
      }
      family$canonical <- ifelse(family$link == "logit", TRUE, FALSE)
    } else if(family$family == "poisson"){
      family$logLik <- function(y, n, mu, wt, dev){
        dpois(y, mu, log = TRUE) * wt
      }
      family$theta <- function(mu){
        log(mu)
      }
      family$cumulant <- function(mu){
        mu
      }
      family$canonical <- ifelse(family$link == "log", TRUE, FALSE)
    } else if(family$family == "Gamma"){
      family$logLik <- function(y, n, mu, wt, dev){
        disp <- dev / sum(wt)
        dgamma(y, 1/disp, scale = mu * disp, log = TRUE) * wt
      }
      family$theta <- function(mu){
        -1/mu
      }
      family$cumulant <- function(mu){
        log(mu)
      }
      family$canonical <- ifelse(family$link == "inverse", TRUE, FALSE)
    } else if(family$family == "inverse.gaussian"){
      family$logLik <- function(y, n, mu, wt, dev){
        disp <- dev / sum(wt)
        -(1/2) * wt * (log(disp * 2 * pi) + 1 + 3 * log(y))
      }
      family$theta <- function(mu){
        -1/mu^2
      }
      family$cumulant <- function(mu){
        -2/mu
      }
      family$canonical <- ifelse(family$link == "1/mu^2", TRUE, FALSE)
    } else {
      stop("'family' not recognized")
    }
    
    family
    
  }