number2color <- 
  function(x, colors, ncol = 21, equidistant = TRUE, xmin = min(x), xmax = max(x)){
    # Map Numbers to Colors
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Update: 2021-11-04
    
    # convert x to numeric
    x <- as.numeric(x)
    
    # check x ranges
    if(any(x < xmin)) warning("Input 'x' contains values less than 'xmin'.")
    if(any(x > xmax)) warning("Input 'x' contains values greater than 'xmax'.")
    
    # check colors
    if(missing(colors)){
      colors <- c("darkblue", rainbow(12)[c(9, 8, 7, 5, 3, 2, 1)], "darkred")
    }
    
    # number of colors
    ncol <- as.integer(ncol[1])
    
    # get 'ncol' colors using colorRampPalette()
    colors <- colorRampPalette(colors)(ncol)
    
    # define breaks
    eps <- sqrt(.Machine$double.eps)
    if(equidistant){
      breaks <- seq(xmin - eps, xmax + eps, length.out = ncol + 1L)
    } else {
      xtemp <- x
      xtemp[which.min(x)] <- xmin - eps
      xtemp[which.max(x)] <- xmax + eps
      breaks <- quantile(xtemp, probs = seq(0, 1, length.out = ncol + 1L))
    }
    
    # bin x values
    xbin <- .bincode(x = x, breaks = breaks, include.lowest = TRUE)
    
    # map to colors
    colors[xbin]
    
  } # end number2color