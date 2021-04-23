basis_nom <- 
  function(x, knots, K = NULL, intercept = FALSE, ridge = FALSE){
    basis.nom(x, knots, K, intercept, ridge)
  }

basis_ord <- 
  function(x, knots, K = NULL, intercept = FALSE, ridge = FALSE){
    basis.ord(x, knots, K, intercept, ridge)
  }

basis_poly <- 
  function(x, knots, m = 2, d = 0, xmin = min(x), xmax = max(x), 
           periodic = FALSE, rescale = FALSE, intercept = FALSE, 
           bernoulli = TRUE, ridge = FALSE){
    basis.poly(x, knots, m, d, xmin, xmax, 
               periodic, rescale, intercept, 
               bernoulli, ridge)
  }

basis_sph <- 
  function(x, knots, m = 2, rescale = TRUE, intercept = FALSE, ridge = FALSE){
    basis.sph(x, knots, m, rescale, intercept, ridge)
  }

basis_tps <- 
  function(x, knots, m = 2, rk = TRUE, intercept = FALSE, ridge = FALSE){
    basis.tps(x, knots, m, rk, intercept, ridge)
  }