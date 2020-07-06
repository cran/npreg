NegBin <-
  function(theta = NULL, link = "log"){
    
    # check theta
    if(is.null(theta)){
      fixed.theta <- FALSE
    } else {
      fixed.theta <- TRUE
      theta <- as.numeric(theta[1])
      if(theta <= 0 | is.infinite(theta)) stop("'theta' must satisfy:  0 < theta < Inf")
      if(is.na(theta) | is.nan(theta)) stop("NA and NaN values are no allowed for 'theta'")
    }
    
    # make link
    linktemp <- substitute(link)
    if (!is.character(linktemp))
      linktemp <- deparse(linktemp)
    if (linktemp %in% c("log", "identity", "sqrt")) 
      stats <- make.link(linktemp)
    else if (is.character(link)) {
      stats <- make.link(link)
      linktemp <- link
    }
    else {
      if (inherits(link, "link-glm")) {
        stats <- link
        if (!is.null(stats$name)) 
          linktemp <- stats$name
      }
      else stop(gettextf("\"%s\" link not available for negative binomial family; available links are \"identity\", \"log\" and \"sqrt\"", 
                         linktemp))
    }
    
    # assign .Theta
    .Theta <- theta
    env <- new.env(parent = .GlobalEnv)
    assign(".Theta", theta, envir = env)
    
    # variance function
    variance <- function(mu) mu + mu^2/.Theta
    
    # validmu function
    validmu <- function(mu) all(mu > 0)
    
    # dev.resids function
    dev.resids <- function(y, mu, wt) {
      2 * wt * (y * log(pmax(1, y)/mu) - (y + .Theta) * log((y + .Theta)/(mu + .Theta)))
    }
    
    # aic function 
    aic <- function(y, n, mu, wt, dev) {
      term <- (y + .Theta) * log(mu + .Theta) - y * log(mu) + 
        lgamma(y + 1) - .Theta * log(.Theta) + lgamma(.Theta) - 
        lgamma(.Theta + y)
      2 * sum(term * wt)
    }
    
    # initialize expression
    initialize <- expression({
      if (any(y < 0)) stop("negative values not allowed for the negative binomial family")
      n <- rep(1, nobs)
      mustart <- y + (y == 0)/6
    })
    
    # simulate function
    simfun <- function(object, nsim) {
      ftd <- fitted(object)
      n <- nsim * length(ftd)
      k <- if (length(n) > 1L) 
        length(n)
      else n
      rpois(k, (ftd * rgamma(k, .Theta))/.Theta)
    }
    
    # logLik function
    logLik <- function(y, n, mu, wt, dev){
      term <- (y + .Theta) * log(mu + .Theta) - y * log(mu) + 
        lgamma(y + 1) - .Theta * log(.Theta) + lgamma(.Theta) - 
        lgamma(.Theta + y)
      -(term * wt)
    }
    
    # canpar function
    canpar <- function(mu){
      log(mu) - log(.Theta + mu)
    }
    
    # cumulant function
    cumulant <- function(mu){
      - .Theta * (log(.Theta) - log(.Theta + mu))
    }
    
    # assign environment
    environment(variance) <- environment(validmu) <- 
      environment(dev.resids) <- environment(aic) <- 
      environment(simfun) <- environment(logLik) <- 
      environment(canpar) <- environment(cumulant) <- env
    
    # return results
    structure(list(family = "NegBin", link = linktemp, linkfun = stats$linkfun, 
                   linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids, 
                   aic = aic, mu.eta = stats$mu.eta, initialize = initialize, 
                   validmu = validmu, valideta = stats$valideta, simulate = simfun,
                   logLik = logLik, canpar = canpar, cumulant = cumulant, 
                   canonical = FALSE, theta = theta, fixed.theta = fixed.theta), 
              class = "family")
    
  } # end NegBin