plot.sm <-
  function(x, terms = x$terms, se = TRUE, n = 201, intercept = FALSE,
           ask = prod(par("mfcol")) < length(terms) && dev.interactive(),
           zero.line = TRUE, zero.lty = 3, zero.col = "black", ncolor = 21, 
           colors = NULL, rev = FALSE, zlim = NULL, lty.col = NULL, 
           legend.xy = "top", main = NULL, xlab = NULL, ylab = NULL, ...){
    # Plot Effects of Smooth Model Fits
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Update: 2024-03-29
    
    
    #########***#########   INITIAL CHECKS   #########***#########
    
    ### check x
    if(!inherits(x, "sm")) stop("Input 'x' must be of class 'sm'.")
    
    ### check terms
    if(is.null(terms)){
      terms <- x$terms
      tid <- 1:length(x$terms)
    } else {
      terms <- as.character(terms)
      tid <- match(terms, x$terms)
      natid <- is.na(tid)
      if(any(natid)) stop("Input 'terms' contains terms not included in the model:\n  ", paste(terms[which(natid)], collapse = ", "))
    }
    nterms0 <- length(x$terms)
    
    ### check se
    se <- as.logical(se[1])
    if(!any(se == c(TRUE, FALSE))) stop("Input 'se' must be TRUE or FALSE")
    
    ### check n
    n <- as.integer(n[1])
    if(n <= 3) stop("Inpit 'n' must be a integer greater than or equal to 3.")
    
    ### check intercept
    intercept <- as.logical(intercept[1])
    if(!any(intercept == c(TRUE, FALSE))) stop("Input 'intercept' must be TRUE or FALSE")
    
    ### check zero.line
    zero.line <- as.logical(zero.line[1])
    if(!any(zero.line == c(TRUE, FALSE))) stop("Input 'zero.line' must be TRUE or FALSE")
    
    ### check zero.lty
    zero.lty <- as.integer(zero.lty[1])
    if(zero.lty <= 0 | zero.lty > 6) stop("Input 'zero.lty' must be a positive integer satisfying:  1 <= zero.lty <= 6")
    
    ### check zero.col
    zero.col <- zero.col[1]
    
    ### check ncolor
    ncolor <- as.integer(ncolor[1])
    if(ncolor < 2L) stop("Input 'ncolor' must be a positive integer (>= 2)")
    
    ### check rev
    rev <- as.logical(rev[1])
    if(!any(rev == c(TRUE, FALSE))) stop("Input 'rev' must be TRUE or FALSE")
    
    ### check colors
    if(is.null(colors)){
      colors <- c("#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7",
                  "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f")
      col <- colorRampPalette(colors)(ncolor)
    } else if(length(colors) == 1L){
      col <- hcl.colors(ncolor, palette = colors, rev = rev)
    } else {
      col <- colorRampPalette(colors)(ncolor)
    }
    
    ### check legend.xy
    legend.xy <- as.character(legend.xy[1])
    
    ### check main and xlab
    inmain <- inxlab <- inylab <- NULL
    if(length(terms) == 1L){
      if(!is.null(main)) inmain <- as.character(main[1])
      if(!is.null(xlab)) inxlab <- as.character(xlab[1])
      if(!is.null(ylab)) inylab <- as.character(ylab[1])
    }
    
    ### reset par on exit
    oldplt <- par()$plt
    oldnew <- par()$new
    on.exit(par(plt = oldplt, new = oldnew))
    
    
    #########***#########   SETUP AND PLOT   #########***#########
    
    ### effect type (main or interaction?)
    nterms <- length(terms)
    effect <- rep(NA, nterms)
    effvar <- vector("list", nterms)
    for(i in 1:nterms){
      effvar[[i]] <- strsplit(terms[i], split = ":")[[1]]
      effect[i] <- length(effvar[[i]])
    } # end for(i in 1:nterms)
    
    ### define types vectors
    factypes <- c("nom", "ord", "ran")
    polytypes <- c("lin", "cub", "qui", "per", "per.lin", "per.cub", "per.qui")
    tpstypes <- c("tps", "tps.lin", "tps.cub", "tps.qui")
    sphtypes <- c("sph", "sph.2", "sph.3", "sph.4")
    
    ### variable names and classes
    nvars <- length(x$types)
    varnames <- names(x$types)
    varclass <- rep(NA, nvars)
    for(j in 1:nvars) varclass[j] <- class(x$data[,varnames[j]])[1]
    
    ### interactive plot setup
    if(i > 1L && ask){
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }
    
    
    #########***#########   PLOT TERMS   #########***#########
    
    ### plotting...
    for(i in 1:nterms){
      
      ## main effect
      if(effect[i] == 1L){
        
        # get variable type
        varid <- match(effvar[[i]], varnames)
        itype <- x$types[varid]
        
        # parametric effects (convert to others)
        if(itype == "par"){
          if(varclass[varid] == "factor"){
            itype <- "nom"
          } else if(varclass[varid] == "ordered"){
            itype <- "ord"
          } else {
            itype <- "cub"
          }
        } # end if(itype == "par")
        
        # 1-dim thin-plate (convert to others)
        if(any(itype == tpstypes) && ncol(x$data[,terms[i]]) == 1L) itype <- "cub"
        
        # nominal and ordinal effects
        if(any(itype == c("nom", "ord"))){
          xlev <- x$specs$xlev[[varid]]
          xlev <- factor(xlev, levels = xlev, ordered = ifelse(itype == "ord", TRUE, FALSE))
          newdata <- data.frame(xlev)
          names(newdata) <- terms[i]
          yhat <- predict(x, newdata = newdata, terms = terms[i], se.fit = se,
                          intercept = intercept)
          dev.hold()
          if(nterms0 == 1L){
            ylab <- substitute(expression(hat(eta) * "("  * tt * ")"), list(tt = terms[i]))
          } else {
            ylab <- substitute(expression(hat(eta)[ii] * "("  * tt * ")"), list(ii = tid[i], tt = terms[i]))
          }
          if(se){
            plotci(x = xlev, y = yhat$fit, se = yhat$se.fit, 
                   xlab = "", ylab = "", 
                   main = ifelse(is.null(inmain), paste("Main effect of", terms[i]), inmain), ...)
            mtext(text = ifelse(is.null(inxlab), terms[i], inxlab), side = 1, line = 2.5, ...)
            mtext(text = ifelse(is.null(inylab), eval(ylab), inylab), side = 2, line = 2.5, ...)
          } else {
            plotci(x = xlev, y = yhat$fit, xlab = "", 
                   ylab = "", main = ifelse(is.null(inmain), paste("Main effect of", terms[i]), inmain), ...)
            mtext(text = ifelse(is.null(inxlab), terms[i], inxlab), side = 1, line = 2.5, ...)
            mtext(text = ifelse(is.null(inylab), eval(ylab), inylab), side = 2, line = 2.5, ...)
          }
          if(zero.line) abline(h = 0, lty = zero.lty, col = zero.col)
          dev.flush()
        } # end if(any(itype == c("nom", "ord")))
        
        # random intercept
        if(itype == "ran"){
          xlev <- x$specs$xlev[[varid]]
          xlev <- factor(xlev, levels = xlev)
          newdata <- data.frame(xlev)
          names(newdata) <- terms[i]
          yhat <- predict(x, newdata = newdata, terms = terms[i], intercept = intercept)
          dev.hold()
          if(nterms0 == 1L){
            ylab <- substitute(expression(hat(eta) * "("  * tt * ")"), list(tt = terms[i]))
          } else {
            ylab <- substitute(expression(hat(eta)[ii] * "("  * tt * ")"), list(ii = tid[i], tt = terms[i]))
          }
          qqnorm(yhat, xlab = "", ylab = "",
                 main = ifelse(is.null(inmain), paste0(terms[i], " Normal Q-Q Plot"), inmain), ...)
          qqline(yhat)
          mtext(text = ifelse(is.null(inxlab), "Theoretical Quantiles", inxlab), side = 1, line = 2.5, ...)
          mtext(text = ifelse(is.null(inylab), eval(ylab), inylab), side = 2, line = 2.5, ...)
          dev.flush()
        } # end if(itype == "ran")
        
        # polynomial
        if(any(itype == polytypes)){
          xrng <- x$specs$xrng[[varid]]
          newdata <- data.frame(seq(xrng[1], xrng[2], length.out = n))
          names(newdata) <- terms[i]
          yhat <- predict(x, newdata = newdata, terms = terms[i], se.fit = se,
                          intercept = intercept)
          dev.hold()
          if(nterms0 == 1L){
            ylab <- substitute(expression(hat(eta) * "("  * tt * ")"), list(tt = terms[i]))
          } else {
            ylab <- substitute(expression(hat(eta)[ii] * "("  * tt * ")"), list(ii = tid[i], tt = terms[i]))
          }
          if(se){
            plotci(x = newdata[,1], y = yhat$fit, se = yhat$se.fit,
                   xlab = "", ylab = "", 
                   main = ifelse(is.null(inmain), paste("Main effect of", terms[i]), inmain), ...)
            mtext(text = ifelse(is.null(inxlab), terms[i], inxlab), side = 1, line = 2.5, ...)
            mtext(text = ifelse(is.null(inylab), eval(ylab), inylab), side = 2, line = 2.5, ...)
          } else {
            plotci(x = newdata[,1], y = yhat$fit, xlab = "", 
                   ylab = "", main = ifelse(is.null(inmain), paste("Main effect of", terms[i]), inmain), ...)
            mtext(text = ifelse(is.null(inxlab), terms[i], inxlab), side = 1, line = 2.5, ...)
            mtext(text = ifelse(is.null(inylab), eval(ylab), inylab), side = 2, line = 2.5, ...)
          }
          if(zero.line) abline(h = 0, lty = zero.lty, col = zero.col)
          dev.flush()
        } # end if(any(itype == polytypes))
        
        # 2-dim thin-plate
        if(any(itype == tpstypes) && ncol(x$data[,terms[i]]) == 2L){
          xname <- colnames(x$data[,terms[i]])
          x1rng <- x$specs$xrng[[varid]][,1]
          x2rng <- x$specs$xrng[[varid]][,2]
          x1seq <- seq(x1rng[1], x1rng[2], length.out = sqrt(n))
          x2seq <- seq(x2rng[1], x2rng[2], length.out = sqrt(n))
          newdata <- as.matrix(expand.grid(x1seq, x2seq))
          colnames(newdata) <- colnames(x$data[,terms[i]])
          newdata <- list(newdata)
          names(newdata) <- terms[i]
          yhat <- predict(x, newdata = newdata, terms = terms[i], se.fit = se,
                          intercept = intercept)
          zmat <- matrix(yhat$fit, nrow = length(x1seq), ncol = length(x2seq))
          dev.hold()
          par(plt = c(0.125, 0.8, oldplt[3:4]))
          image(x = x1seq, y = x2seq, z = zmat, col = col,
                xlab = "",  ylab = "", 
                main = ifelse(is.null(inmain), paste("Main effect of", terms[i]), inmain), ...)
          mtext(text = ifelse(is.null(inxlab), paste(terms[i], xname[1], sep = "."), inxlab), side = 1, line = 2.5, ...)
          mtext(text = ifelse(is.null(inylab), paste(terms[i], xname[2], sep = "."), inylab), side = 2, line = 2.5, ...)
          if(nterms0 == 1L){
            zlab <- substitute(hat(eta) * "("  * tt * ")", list(tt = terms[i]))
          } else {
            zlab <- substitute(hat(eta)[ii] * "("  * tt * ")", list(ii = tid[i], tt = terms[i]))
          }
          par(plt = c(0.85, 0.9, oldplt[3:4]), new = TRUE)
          color.legend(zlim = range(zmat), col = col, ncol = ncolor, zlab = zlab, ...)
          dev.flush()
          par(plt = oldplt, new = oldnew)
        } # end if(any(itype == tpstypes) && ncol(x$data[,terms[i]]) == 2L)
        
        # 3-dim thin-plate
        if(any(itype == tpstypes) && ncol(x$data[,terms[i]]) == 3L){
          warning("Plots of 3-dimensional thin-plate splines are not supported.")
        } # end if(any(itype == tpstypes) && ncol(x$data[,terms[i]]) == 3L)
        
        # spherical spline
        if(any(itype == sphtypes)){
          xname <- colnames(x$data[,terms[i]])
          x1rng <- x$specs$xrng[[varid]][,1]
          x2rng <- x$specs$xrng[[varid]][,2]
          x1seq <- seq(x1rng[1], x1rng[2], length.out = sqrt(n))
          x2seq <- seq(x2rng[1], x2rng[2], length.out = sqrt(n))
          newdata <- as.matrix(expand.grid(x1seq, x2seq))
          colnames(newdata) <- colnames(x$data[,terms[i]])
          newdata <- list(newdata)
          names(newdata) <- terms[i]
          yhat <- predict(x, newdata = newdata, terms = terms[i], se.fit = se,
                          intercept = intercept)
          zmat <- matrix(yhat$fit, nrow = length(x1seq), ncol = length(x2seq))
          dev.hold()
          par(plt = c(0.125, 0.8, oldplt[3:4]))
          image(x = x1seq, y = x2seq, z = zmat, col = col,
                xlab = "",  ylab = "", 
                main = ifelse(is.null(inmain), paste("Main effect of", terms[i]), inmain), ...)
          mtext(text = ifelse(is.null(inxlab), paste(terms[i], xname[1], sep = "."), inxlab), side = 1, line = 2.5, ...)
          mtext(text = ifelse(is.null(inylab), paste(terms[i], xname[2], sep = "."), inylab), side = 2, line = 2.5, ...)
          if(nterms0 == 1L){
            zlab <- substitute(hat(eta) * "("  * tt * ")", list(tt = terms[i]))
          } else {
            zlab <- substitute(hat(eta)[ii] * "("  * tt * ")", list(ii = tid[i], tt = terms[i]))
          }
          par(plt = c(0.85, 0.9, oldplt[3:4]), new = TRUE)
          color.legend(zlim = range(zmat), col = col, ncol = ncolor, zlab = zlab, ...)
          dev.flush()
          par(plt = oldplt, new = oldnew)
        } # end if(any(itype == tpstypes) && ncol(x$data[,terms[i]]) == 2L)
        
      } # end if(effect[i] == 1L)
      
      ## 2-way interaction effect
      if(effect[i] == 2L){
        
        # get variable type
        varid <- match(effvar[[i]], varnames)
        itype <- x$types[varid]
        
        # parametric effects (convert to others)
        for(j in 1:2){
          if(itype[j] == "par"){
            if(varclass[varid[j]] == "factor"){
              itype[j] <- "nom"
            } else if(varclass[varid[j]] == "ordered"){
              itype[j] <- "ord"
            } else {
              itype[j] <- "cub"
            }
          } # end if(itype[j] == "par")
        } # end for(j in 1:2)
        
        # 1-dim thin-plate (convert to others)
        for(j in 1:2){
          if(any(itype[j] == tpstypes) && ncol(x$data[,effvar[[i]]][j]) == 1L) itype[j] <- "cub"
        }
        
        # factor by factor
        if(any(itype[1] == factypes) && any(itype[2] == factypes)){
          xlev1 <- x$specs$xlev[[varid[1]]]
          xlev1 <- factor(xlev1, levels = xlev1, ordered = ifelse(itype[1] == "ord", TRUE, FALSE))
          xlev2 <- x$specs$xlev[[varid[2]]]
          xlev2 <- factor(xlev2, levels = xlev2, ordered = ifelse(itype[2] == "ord", TRUE, FALSE))
          newdata <- expand.grid(xlev1, xlev2)
          names(newdata) <- effvar[[i]]
          yhat <- predict(x, newdata = newdata, terms = terms[i], intercept = intercept)
          zmat <- matrix(yhat, nrow = length(xlev1), ncol = length(xlev2))
          dev.hold()
          par(plt = c(0.125, 0.8, oldplt[3:4]))
          image(x = as.integer(xlev1), y = as.integer(xlev2),
                z = zmat, xlab = "", ylab = "",
                main = ifelse(is.null(inmain), paste("Interaction effect of", terms[i]), inmain), 
                col = col, axes = FALSE, ...)
          mtext(text = ifelse(is.null(inxlab), effvar[[i]][1], inxlab), side = 1, line = 2.5, ...)
          mtext(text = ifelse(is.null(inylab), effvar[[i]][2], inylab), side = 2, line = 2.5, ...)
          axis(1, at = as.integer(xlev1), labels = xlev1, ...)
          axis(2, at = as.integer(xlev2), labels = xlev2, ...)
          box()
          zlab <- substitute(hat(eta)[ii] * "("  * xx * ", " * yy * ")", 
                             list(ii = tid[i], xx = effvar[[i]][1], yy = effvar[[i]][2]))
          par(plt = c(0.85, 0.9, oldplt[3:4]), new = TRUE)
          color.legend(zlim = range(zmat), col = col, ncol = ncolor, zlab = zlab, ...)
          dev.flush()
          par(plt = oldplt, new = oldnew)
        } # end if(any(itype[1] == factypes) && any(itype[2] == factypes))
        
        # polynomial by factor
        if(any(itype[1] == polytypes) && any(itype[2] == factypes)){
          x1rng <- x$specs$xrng[[varid[1]]]
          x1seq <- seq(x1rng[1], x1rng[2], length.out = n)
          xlev2 <- x$specs$xlev[[varid[2]]]
          xlev2 <- factor(xlev2, levels = xlev2, ordered = ifelse(itype[2] == "ord", TRUE, FALSE))
          ylab <- substitute(expression(hat(eta)[ii] * "("  * xx * ", " * yy * ")"), 
                             list(ii = tid[i], xx = effvar[[i]][1], yy = effvar[[i]][2]))
          yhat <- vector("list", length(xlev2))
          for(j in 1:length(xlev2)){
            newdata <- data.frame(x1seq, xlev2[j])
            colnames(newdata) <- effvar[[i]]
            yhat[[j]] <- predict(x, newdata = newdata, terms = terms[i], intercept = intercept)
          }
          ylim <- range(unlist(yhat)) * c(1, 1.2)
          dev.hold()
          for(j in 1:length(xlev2)){
            plotci(x = newdata[,1], y = yhat[[j]], 
                   add = ifelse(j > 1, TRUE, FALSE), col.ci = NA, ylim = ylim,
                   xlab = "", ylab = "", lty = j, col = j,
                   main = ifelse(is.null(inmain), paste("Interaction effect of", terms[i]), inmain), ...)
            mtext(text = ifelse(is.null(inxlab), effvar[[i]][1], inxlab), side = 1, line = 2.5, ...)
            mtext(text = ifelse(is.null(inylab), eval(ylab), inylab), side = 2, line = 2.5, ...)
          }
          legend(legend.xy, legend = paste(effvar[[i]][2], xlev2, sep = " = "), 
                 lty = 1:length(xlev2), col = 1:length(xlev2), lwd = 2, bty = "n")
          dev.flush()
        } # end if(any(itype[1] == polytypes) && any(itype[2] == factypes))
        
        # factor by polynomial (same as previous)
        if(any(itype[1] == factypes) && any(itype[2] == polytypes)){
          xlev1 <- x$specs$xlev[[varid[1]]]
          xlev1 <- factor(xlev1, levels = xlev1, ordered = ifelse(itype[1] == "ord", TRUE, FALSE))
          x2rng <- x$specs$xrng[[varid[2]]]
          x2seq <- seq(x2rng[1], x2rng[2], length.out = n)
          ylab <- substitute(expression(hat(eta)[ii] * "("  * xx * ", " * yy * ")"), 
                             list(ii = tid[i], xx = effvar[[i]][1], yy = effvar[[i]][2]))
          yhat <- vector("list", length(xlev1))
          for(j in 1:length(xlev1)){
            newdata <- data.frame(xlev1[j], x2seq)
            colnames(newdata) <- effvar[[i]]
            yhat[[j]] <- predict(x, newdata = newdata, terms = terms[i], intercept = intercept)
          }
          ylim <- range(unlist(yhat)) * c(1, 1.2)
          dev.hold()
          for(j in 1:length(xlev1)){
            plotci(x = newdata[,2], y = yhat[[j]], 
                   add = ifelse(j > 1, TRUE, FALSE), col.ci = NA, ylim = ylim,
                   xlab = "", ylab = "", lty = j, col = j,
                   main = ifelse(is.null(inmain), paste("Interaction effect of", terms[i]), inmain), ...)
            mtext(text = ifelse(is.null(inxlab), effvar[[i]][2], inxlab), side = 1, line = 2.5, ...)
            mtext(text = ifelse(is.null(inylab), eval(ylab), inylab), side = 2, line = 2.5, ...)
          }
          legend(legend.xy, legend = paste(effvar[[i]][1], xlev1, sep = " = "), 
                 lty = 1:length(xlev1), col = 1:length(xlev1), lwd = 2, bty = "n")
          dev.flush()
        } # end if(any(itype[1] == factypes) && any(itype[2] == polytypes))
        
        # polynomial by polynomial
        if(any(itype[1] == polytypes) && any(itype[2] == polytypes)){
          x1rng <- x$specs$xrng[[varid[1]]]
          x2rng <- x$specs$xrng[[varid[2]]]
          x1seq <- seq(x1rng[1], x1rng[2], length.out = sqrt(n))
          x2seq <- seq(x2rng[1], x2rng[2], length.out = sqrt(n))
          newdata <- expand.grid(x1seq, x2seq)
          colnames(newdata) <- effvar[[i]]
          yhat <- predict(x, newdata = newdata, terms = terms[i],
                          intercept = intercept)
          zmat <- matrix(yhat, nrow = length(x1seq), ncol = length(x2seq))
          dev.hold()
          par(plt = c(0.125, 0.8, oldplt[3:4]))
          image(x = x1seq, y = x2seq, z = zmat, col = col,
                xlab = "", ylab = "", 
                main = ifelse(is.null(inmain), paste("Interaction effect of", terms[i]), inmain), ...)
          mtext(text = ifelse(is.null(inxlab), effvar[[i]][1], inxlab), side = 1, line = 2.5, ...)
          mtext(text = ifelse(is.null(inylab), effvar[[i]][2], inylab), side = 2, line = 2.5, ...)
          zlab <- substitute(hat(eta)[ii] * "("  * xx * ", " * yy * ")", 
                             list(ii = tid[i], xx = effvar[[i]][1], yy = effvar[[i]][2]))
          par(plt = c(0.85, 0.9, oldplt[3:4]), new = TRUE)
          color.legend(zlim = range(zmat), col = col, 
                       ncol = ncolor, zlab = zlab, ...)
          dev.flush()
          par(plt = oldplt, new = oldnew)
        } # end if(any(itype[1] == polytypes) && any(itype[2] == polytypes))
        
      } # end if(effect[i] == 2L)
      
      ## 3-way interaction effect
      if(effect[i] == 3L){
        warning("Plots of 3-way interaction effects are not supported.")
      } # end if(effect[i] == 3L)
      
    } # end for(i in 1:nterms)
    
  } # end plot.sm