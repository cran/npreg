linkinvd <- 
  function(link){
    # first derivative of linkinv
    
    # check link
    links <- c("logit", "probit", "cauchit", "cloglog", "identity", "log", "sqrt", "1/mu^2", "inverse")
    link <- as.character(link[1])
    linkmatch <- pmatch(link, links)
    if(is.na(linkmatch)) stop("Invalid 'link' function.")
    link <- links[linkmatch]
    
    # derivative of link
    if(link == "logit"){
      lid <- function(eta) exp(eta) / (1 + exp(eta))^2
    } else if(link == "probit"){
      lid <- function(eta) dnorm(eta)
    } else if(link == "cauchit"){
      lid <- function(eta) 1 / (pi * (eta^2 + 1))
    } else if(link == "cloglog"){
      lid <- function(eta) exp(eta - exp(eta))
    } else if (link == "identity"){
      lid <- function(eta) rep(1, times = length(eta))
    } else if(link == "log"){
      lid <- function(eta) exp(eta)
    } else if(link == "sqrt"){
      lid <- function(eta) 2 * eta
    } else if(link == "1/mu^2"){
      lid <- function(eta) 0 - 1 / (2 * eta^(3/2))
    } else if(link == "inverse"){
      lid <- function(eta) 0 - 1 / eta^2
    }
    
    return(lid)
    
  }