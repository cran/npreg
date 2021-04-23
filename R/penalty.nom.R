penalty.nom <- 
  function(x, K = NULL){
    # Nominal Smoothing Spline Penalty
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Update: 2021-04-09
    
    if(is.null(K)) K <- length(unique(x))
    outer(x, x, FUN = "==") - 1/K
    
  }