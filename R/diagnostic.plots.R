diagnostic.plots <- 
  function(x, which = c(1, 2, 3, 5), 
           caption = list("Residuals vs Fitted", 
                          "Normal Q-Q", "Scale-Location", 
                          "Cook's distance", "Residuals vs Leverage", 
                          "Cook's dist vs Variance ratio"), 
           panel = if (add.smooth) function(x, y, ...) panel.smooth(x, y, iter = iter.smooth, ...) else points, sub.caption = NULL, 
           main = "", ask = prod(par("mfcol")) < length(which) && dev.interactive(), 
           ..., id.n = 3, labels.id = names(residuals(x)), cex.id = 0.75, cex.pt = 1,
           qqline = TRUE, cook.levels = c(0.5, 1), add.smooth = getOption("add.smooth"), 
           iter.smooth = if (isGlm) 0 else 3, label.pos = c(4, 2), cex.caption = 1, 
           cex.oma.main = 1.25, cex.lab = 1, line.lab = 3, xlim = NULL, ylim = NULL){
    
    cex.lab <- cex.lab * par("cex")
    
    # define dropInf function
    dropInf <- function(x, h) {
      if (any(isInf <- h >= 1)) {
        warning(gettextf("not plotting observations with leverage one:\n  %s", 
                         paste(which(isInf), collapse = ", ")), call. = FALSE, 
                domain = NA)
        x[isInf] <- NaN
      }
      x
    }
    
    # check x
    xc <- class(x)
    if(!any(xc == c("ss", "sm", "gsm")))
      stop("Input 'x' must be an object of class 'ss' or 'sm' or 'gsm'")
    
    # check which
    which <- unique(as.integer(which))
    if(any(which < 1L) | any(which > 6L))
      stop("Input 'which' must be a subset of 1:6")
    nwhich <- length(which)
    
    # check xlim
    if(is.null(xlim)){
      xlim <- vector("list", nwhich)
    } else {
      if(nwhich == 1L) xlim <- list(unlist(xlim))
      xlim <- as.list(xlim)
      if(length(xlim) != nwhich) stop("length(which) must equal length(xlim)")
      for(i in 1:nwhich){
        if(!is.null(xlim[[i]])){
          xlim[[i]] <- as.numeric(xlim[[i]])
          if(length(xlim[[i]]) != 2L) stop("xlim[[", i, "]] must contain a vector of length 2 of the form c(xmin, xmax)")
          if(xlim[[i]][1] >= xlim[[i]][2]) stop("xlim[[", i, "]] must contain a vector of the form c(xmin, xmax) with xmin < xmax")
        }
      }
    }
    xlimits <- vector("list", 6L)
    xlimits[which] <- xlim
    
    # check ylim
    if(is.null(ylim)){
      ylim <- vector("list", nwhich)
    } else {
      if(nwhich == 1L) ylim <- list(unlist(ylim))
      ylim <- as.list(ylim)
      if(length(ylim) != nwhich) stop("length(which) must equal length(ylim)")
      for(i in 1:nwhich){
        if(!is.null(ylim[[i]])){
          ylim[[i]] <- as.numeric(ylim[[i]])
          if(length(ylim[[i]]) != 2L) stop("ylim[[", i, "]] must contain a vector of length 2 of the form c(ymin, ymax)")
          if(ylim[[i]][1] >= ylim[[i]][2]) stop("ylim[[", i, "]] must contain a vector of the form c(ymin, ymax) with ymin < ymax")
        }
      }
    }
    ylimits <- vector("list", 6L)
    ylimits[which] <- ylim
    
    # is gsm?
    isGlm <- inherits(x, "gsm")
    
    # which plots to show
    show <- rep(FALSE, 6)
    show[which] <- TRUE
    
    # get residuals
    r <- if (isGlm) residuals(x, type = "pearson") else residuals(x)
    
    # get n
    n <- length(r)
    
    # get linear predictors
    yh <- if(isGlm) x$linear.predictors else fitted(x)
    
    # get weights
    w <- weights(x)
    wind <- w != 0
    r <- r[wind]
    yh <- yh[wind]
    w <- w[wind]
    labels.id <- labels.id[wind]
    
    # need 2:6?
    if (any(show[2L:6L])) {
      s <- if(isGlm) sqrt(x$dispersion) else x$sigma
      hii <- hatvalues(x)
      if (any(show[4L:6L])) cook <- cooks.distance(x, hat = hii)
    }
    
    # need c(2,3,5)?
    if (any(show[c(2L, 3L, 5L)])) {
      ylab5 <- ylab23 <- ifelse(isGlm, 
                                "Std. Pearson resid.",
                                "Standardized residuals")
      r.w <- sqrt(w) * r
      rsp <- rs <- dropInf(if (isGlm) 
        rstandard(x, type = "pearson")
        else r.w / (s * sqrt(1 - hii)), hii)
    }
    
    # need 5:6?
    if (any(show[5L:6L])) {
      r.hat <- range(hii, na.rm = TRUE)
      isConst.hat <- all(r.hat == 0) || diff(r.hat) < 1e-10 * mean(hii, na.rm = TRUE)
    }
    
    # need c(1, 3)?
    if (any(show[c(1L, 3L)])) 
      l.fit <- ifelse(isGlm, "Predicted values", "Fitted values")
    
    # check id.n
    if (is.null(id.n)) {
      id.n <- 0
    } else {
      id.n <- as.integer(id.n)
      if (id.n < 0L || id.n > n) 
        stop(gettextf("'id.n' must be in {1,..,%d}", n), domain = NA)
    }
    
    # non-zero id.n?
    if (id.n > 0L) {
      if (is.null(labels.id)) labels.id <- paste(1L:n)
      iid <- 1L:id.n
      show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
      if (any(show[2L:3L])) 
        show.rs <- sort.list(abs(rs), decreasing = TRUE)[iid]
      text.id <- function(x, y, ind, adj.x = TRUE) {
        labpos <- if (adj.x) 
          label.pos[1 + as.numeric(x > mean(range(x)))]
        else 3
        text(x, y, labels.id[ind], cex = cex.id, xpd = TRUE, 
             pos = labpos, offset = 0.25)
      }
    }
    
    # caption function
    getCaption <- function(k) if (length(caption) < k) 
      NA_character_
    else as.graphicsAnnot(caption[[k]])
    if (is.null(sub.caption)) {
      cal <- x$call
      if (!is.na(m.f <- match("formula", names(cal)))) {
        cal <- cal[c(1, m.f)]
        names(cal)[2L] <- ""
      }
      cc <- deparse(cal, 80)
      nc <- nchar(cc[1L], "c")
      abbr <- length(cc) > 1 || nc > 75
      sub.caption <- if (abbr) 
        paste(substr(cc[1L], 1L, min(75L, nc)), "...")
      else cc[1L]
    }
    
    # setup fig
    one.fig <- prod(par("mfcol")) == 1
    if (ask) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }
    
    # residuals vs fitted values
    if (show[1L]) {
      xlim <- if(is.null(xlimits[[1]])) extendrange(yh, f = 0.001) else xlimits[[1]]
      if(is.null(ylimits[[1]])){
        ylim <- range(r, na.rm = TRUE)
        if (id.n > 0) ylim <- extendrange(r = ylim, f = 0.08)
      } else {
        ylim <- ylimits[[1]]
      }
      dev.hold()
      plot(yh, r, xlab = "", ylab = "", main = main, 
           xlim = xlim, ylim = ylim, type = "n", ...)
      mtext(text = l.fit, side = 1, line = line.lab, cex = cex.lab, ...)
      mtext(text = "Residuals", side = 2, line = line.lab, cex = cex.lab, ...)
      panel(yh, r, cex = cex.pt, ...)
      if (one.fig) 
        title(sub = sub.caption, ...)
      mtext(getCaption(1), 3, 0.25, cex = cex.caption)
      if (id.n > 0) {
        y.id <- r[show.r]
        y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
        text.id(yh[show.r], y.id, show.r)
      }
      abline(h = 0, lty = 3, col = "gray")
      dev.flush()
    }
    
    # normal QQ plot
    if (show[2L]) {
      xlim <- if(is.null(xlimits[[2]])) extendrange(qqnorm(rs, plot.it = FALSE)$x, f = 0.001) else xlimits[[2]]
      if(is.null(ylimits[[2]])){
        ylim <- range(rs, na.rm = TRUE)
        ylim[2L] <- ylim[2L] + diff(ylim) * 0.075
      } else {
        ylim <- ylimits[[2]]
      }
      dev.hold()
      qq <- qqnorm(rs, main = main, xlim = xlim, ylim = ylim, 
                   xlab = "", ylab = "", cex = cex.pt, ...)
      mtext(text = "Theoretical Quantiles", side = 1, line = line.lab, cex = cex.lab, ...)
      mtext(text = ylab23, side = 2, line = line.lab, cex = cex.lab, ...)
      if (qqline) 
        qqline(rs, lty = 3, col = "gray50")
      if (one.fig) 
        title(sub = sub.caption, ...)
      mtext(getCaption(2), 3, 0.25, cex = cex.caption)
      if (id.n > 0) 
        text.id(qq$x[show.rs], qq$y[show.rs], show.rs)
      dev.flush()
    }
    
    # scale-location
    if (show[3L]) {
      sqrtabsr <- sqrt(abs(rs))
      ylim <- if(is.null(ylimits[[3]])) c(0, max(sqrtabsr, na.rm = TRUE)) else ylimits[[3]]
      yl <- as.expression(substitute(sqrt(abs(YL)), list(YL = as.name(ylab23))))
      yhn0 <- if (is.null(w)) 
        yh
      else yh[w != 0]
      xlim <- if(is.null(xlimits[[3]])) extendrange(yhn0, f = 0.001) else xlimits[[3]]
      dev.hold()
      plot(yhn0, sqrtabsr, xlab = "", ylab = "", main = main, 
           xlim = xlim, ylim = ylim, type = "n", ...)
      mtext(text = l.fit, side = 1, line = line.lab, cex = cex.lab, ...)
      mtext(text = yl, side = 2, line = line.lab, cex = cex.lab, ...)
      panel(yhn0, sqrtabsr, cex = cex.pt, ...)
      if (one.fig) 
        title(sub = sub.caption, ...)
      mtext(getCaption(3), 3, 0.25, cex = cex.caption)
      if (id.n > 0) 
        text.id(yhn0[show.rs], sqrtabsr[show.rs], show.rs)
      dev.flush()
    }
    
    # Cook's distance
    if (show[4L]) {
      if (id.n > 0) {
        show.r <- order(-cook)[iid]
        ymx <- cook[show.r[1L]] * 1.075
      }
      else ymx <- max(cook, na.rm = TRUE)
      dev.hold()
      xlim <- if(is.null(xlimits[[4]])) extendrange(c(1, length(cook)), f = 0.001) else xlimits[[4]]
      ylim <- if(is.null(ylimits[[4]])) c(0, ymx) else ylimits[[4]]
      plot(cook, type = "h", xlim = xlim, ylim = ylim, main = main, 
           xlab = "", ylab = "", cex = cex.pt, ...)
      mtext(text = "Obs. number", side = 1, line = line.lab, cex = cex.lab, ...)
      mtext(text = "Cook's distance", side = 2, line = line.lab, cex = cex.lab, ...)
      if (one.fig) 
        title(sub = sub.caption, ...)
      mtext(getCaption(4), 3, 0.25, cex = cex.caption)
      if (id.n > 0) 
        text.id(show.r, cook[show.r], show.r, adj.x = FALSE)
      dev.flush()
    }
    
    # residuals vs leverages
    if (show[5L]) {
      if(is.null(ylimits[[5]])){
        ylim <- range(rsp, na.rm = TRUE)
        if (id.n > 0) {
          ylim <- extendrange(r = ylim, f = 0.08)
          show.rsp <- order(-cook)[iid]
        }
      } else {
        ylim <- ylimits[[5]]
      }
      if (id.n > 0) {
        #ylim <- extendrange(r = ylim, f = 0.08)
        show.rsp <- order(-cook)[iid]
      }
      do.plot <- TRUE
      xx <- hii
      xx[xx >= 1] <- NA
      xlim <- if(is.null(xlimits[[5]])) c(0, max(xx, na.rm = TRUE)) else xlimits[[5]]
      dev.hold()
      plot(xx, rsp, xlim = xlim, ylim = ylim, main = main, xlab = "", 
           ylab = "", type = "n", ...)
      mtext(text = "Leverage", side = 1, line = line.lab, cex = cex.lab, ...)
      mtext(text = ylab5, side = 2, line = line.lab, cex = cex.lab, ...)
      panel(xx, rsp, cex = cex.pt, ...)
      abline(h = 0, v = 0, lty = 3, col = "gray")
      if (one.fig) 
        title(sub = sub.caption, ...)
      if (length(cook.levels)) {
        p <- x$df
        usr <- par("usr")
        hh <- seq.int(min(r.hat[1L], r.hat[2L]/100), 
                      usr[2L], length.out = 101)
        for (crit in cook.levels) {
          cl.h <- sqrt(crit * p * (1 - hh)/hh)
          lines(hh, cl.h, lty = 2, col = 2)
          lines(hh, -cl.h, lty = 2, col = 2)
        }
        legend("bottomleft", legend = "Cook's distance", 
               lty = 2, col = 2, bty = "n")
        xmax <- min(0.99, usr[2L])
        ymult <- sqrt(p * (1 - xmax)/xmax)
        aty <- sqrt(cook.levels) * ymult
        axis(4, at = c(-rev(aty), aty), labels = paste(c(rev(cook.levels), 
                                                         cook.levels)), mgp = c(0.25, 0.25, 0), las = 2, 
             tck = 0, cex.axis = cex.id, col.axis = 2)
      }
      dev.flush()
      if (do.plot) {
        mtext(getCaption(5), 3, 0.25, cex = cex.caption)
        if (id.n > 0) {
          y.id <- rsp[show.rsp]
          y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
          text.id(xx[show.rsp], y.id, show.rsp)
        }
      }
    }
    
    # expression("Cook's dist vs Leverage  " * s[ii]/(1 - s[ii])))
    if (show[6L]) {
      g <- dropInf(hii/(1 - hii), hii)
      ymx <- max(cook, na.rm = TRUE) * 1.025
      xlim <- if(is.null(xlimits[[6]])) extendrange(c(0, max(g, na.rm = TRUE)), f = 0.001) else xlimits[[6]]
      ylim <- if(is.null(ylimits[[6]])) c(0, ymx) else ylimits[[6]]
      dev.hold()
      plot(g, cook, xlim = xlim, ylim = ylim, 
           main = main, ylab = "", 
           xlab = "", type = "n", ...)
      mtext(text = expression("Variance ratio:  " * s[ii] / (1 - s[ii]) ), side = 1, line = line.lab + 0.2, cex = cex.lab, ...)
      mtext(text = "Cook's distance", side = 2, line = line.lab, cex = cex.lab, ...)
      panel(g, cook, cex = cex.pt, ...)
      #athat <- pretty(hii)
      #axis(1, at = athat/(1 - athat), labels = paste(athat))
      if (one.fig) 
        title(sub = sub.caption, ...)
      p <- x$df
      #bval <- pretty(sqrt(p * cook/g), 5)
      bval <- pretty(abs(rstandard(x)), 5)
      usr <- par("usr")
      xmax <- usr[2L]
      ymax <- usr[4L]
      for (i in seq_along(bval)) {
        bi2 <- bval[i]^2
        if (p * ymax > bi2 * xmax) {
          xi <- xmax + strwidth(" ")/3
          yi <- bi2 * xi/p
          #abline(0, bi2, lty = 2)
          lines(c(0, xi), c(0, yi), lty = 2)
          text(xi, yi, paste(bval[i]), adj = 0, xpd = TRUE, cex = cex.id * 0.8)
        }
        else {
          yi <- ymax - 1.5 * strheight(" ")
          xi <- p * yi/bi2
          lines(c(0, xi), c(0, yi), lty = 2)
          text(xi, ymax - 0.8 * strheight(" "), paste(bval[i]), 
               adj = 0.5, xpd = TRUE, cex = cex.id * 0.8)
        }
      }
      mtext(getCaption(6), 3, 0.25, cex = cex.caption)
      if (id.n > 0) {
        show.r <- order(-cook)[iid]
        text.id(g[show.r], cook[show.r], show.r)
      }
      dev.flush()
    }
    
    # reset and exit
    if (!one.fig && par("oma")[3L] >= 1) 
      mtext(sub.caption, outer = TRUE, cex = cex.oma.main)
    invisible()
    
  } # end diagnostic.plots