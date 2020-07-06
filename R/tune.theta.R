tune.theta <-
  function(logtheta, y, n, mu, wt){
    
    theta <- exp(logtheta)
    probs <- theta / (theta + mu)
    aterm <- sum(wt * lgamma(theta)) / n
    bterm <- sum(wt * lgamma(theta + y)) / n
    cterm <- sum(wt * theta * log(probs)) / n
    dterm <- sum(wt * y * log(1 - probs)) / n
    val <- aterm - bterm - cterm - dterm
    return(val)
    
  }