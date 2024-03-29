color.legend <-
  function(zlim, side = 4, col = NULL, ncol = NULL, zlab = "z", 
           zline = 2.5, box = TRUE, zcex = 1, ...){
    # Color Legend for Image Plot
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Update: 2024-03-29
    
    
    #########***#########   INITIAL CHECKS   #########***#########
    
    # zlim
    zlim <- as.numeric(zlim)
    if(length(zlim) != 2L) stop("Input 'zlim' must be a vector of length two:  zlim = range(z)")
    
    # check side
    side <- as.integer(side[1])
    if(side < 1L | side > 4L) stop("Input 'side' must be an integer: 1, 2, 3, 4")
    
    # check col
    if(is.null(col)){
      col <- c("#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7",
               "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f")
    } else {
      col <- as.character(col)
    }
    
    # check ncol
    if(is.null(ncol)){
      ncol <- length(col)
    } else {
      ncol <- as.integer(ncol[1])
      if(ncol < 1L) stop("Input 'ncol' must be a positive integer.")
    }
    
    # create color palette
    col <- colorRampPalette(col)(ncol)
    
    
    #########***#########   PLOTTING   #########***#########
    
    # zseq for plotting
    zseq <- seq(zlim[1], zlim[2], length.out = ncol)
    
    # draw legend
    if(any(side == c(2, 4))){
      image(x = 1, y = zseq,
            z = matrix(zseq, nrow = 1, ncol = ncol), 
            xlab = "", ylab = "", col = col, axes = FALSE)
    } else {
      image(x = zseq, y = 1,
            z = matrix(zseq, nrow = ncol, ncol = 1), 
            xlab = "", ylab = "", col = col, axes = FALSE)
    }
    
    # box?
    if(box) box()
    
    # add axis
    axis(side, ...)
    mtext(zlab, side = side, line = zline, cex = zcex, ...)
    
  }