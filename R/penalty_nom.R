penalty_nom <- 
  function(x, K = NULL){
    # Nominal Smoothing Spline Penalty
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Update: 2019-04-04
    
    if(is.null(K)) K <- length(unique(x))
    outer(x, x, FUN = "==") - 1/K
    
  }