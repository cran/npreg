penalty_nom <-
  function(x, K = NULL){
    penalty.nom(x, K)
  }

penalty_ord <-
  function(x, K = NULL, xlev = NULL){
    penalty.ord(x, K, xlev)
  }

penalty_poly <-
  function(x, m = 2, xmin = min(x), xmax = max(x), 
           periodic = FALSE, rescale = FALSE, bernoulli = TRUE){
    penalty.poly(x, m, xmin, xmax, periodic, rescale, bernoulli)
  }

penalty_sph <-
  function(x, m = 2){
    penalty.sph(x, m)
  }

penalty_tps <-
  function(x, m = 2, rk = TRUE){
    penalty.tps(x, m, rk)
  }