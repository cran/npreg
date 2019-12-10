summary.sm <- 
  function(object, ...){
    # summary method for class "sm"
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2019-09-24
    
    # get residuals
    resid <- object$data[,1] - object$fitted.values
    
    # get diagnostics
    fit.terms <- scale(predict(object, type = "terms"), scale = FALSE)
    nterms <- ncol(fit.terms)
    Cmat <- matrix(0, nrow = nterms, ncol = nterms)
    for(i in 1:nterms){
      for(j in 1:i){
        Cmat[i,j] <- Cmat[j,i] <- sum(fit.terms[,i] * fit.terms[,j]) / sqrt(sum(fit.terms[,i]^2) * sum(fit.terms[,j]^2))
      }
    }
    kappa <- sqrt(diag(psolve(Cmat)))
    fitc <- rowSums(fit.terms)
    pi <- as.numeric(crossprod(fit.terms, fitc) / sum(fitc^2))
    names(kappa) <- names(pi) <- object$terms
    
    # adjusted R-squared
    adjRsq <- 1 - object$sigma^2 / var(object$data[,1])
    
    # some basic info
    nobs <- length(object$fitted.values)
    nullindx <- 1:object$nsdf
    erdf <- nobs - object$df
    
    # parametric coefficients table
    p.coef <- object$coefficients[nullindx]
    p.ster <- sqrt(rowSums(object$cov.sqrt[nullindx,,drop=FALSE]^2))
    p.tval <- p.coef / p.ster
    p.pval <- 2 * (1 - pt(abs(p.tval), df = erdf))
    p.table <- data.frame(p.coef, p.ster, p.tval, p.pval)
    rownames(p.table) <- names(p.coef)
    colnames(p.table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    
    # smooth term names
    s.coef <- object$coefficients[-nullindx]
    snames <- names(s.coef)
    mxknot <- max(sapply(object$specs$knots, function(x) nrow(as.matrix(x))))
    for(m in mxknot:1){
      snames <- gsub(paste0(".knot.", m), "", snames, fixed = TRUE)
    }
    
    # smooth term degrees of freedom and SS
    ss <- s.df <- rep(0, nterms)
    if(object$specs$tprk){
      
      Xmats <- predict(object, combine = FALSE, design = TRUE)
      Xpsvd <- svd(Xmats$X$p, nv = 0)
      Xmats$X$s <- Xmats$X$s - Xpsvd$u %*% crossprod(Xpsvd$u, Xmats$X$s)
      ss <- sum((psolve(Xmats$X$s %*% object$cov.sqrt[-nullindx,]) %*% (Xmats$X$s %*% s.coef))^2) * object$sigma^2
      pis <- rep(0, nterms)
      for(i in 1:nterms){
        fitX <- predict(object, terms = object$terms[i],
                        combine = FALSE, design = TRUE)
        fitX$X$s <- fitX$X$s - Xpsvd$u %*% crossprod(Xpsvd$u, fitX$X$s)
        pis[i] <- sum((fitX$X$s %*% s.coef)^2)
        s.df[i] <- sum(rowSums((fitX$X$s %*% object$cov.sqrt[-nullindx,] / object$sigma)^2))
      }
      pis <- pis / sum(pis)
      ss <- ss * pis
      
    } else {
      
      Xmats <- predict(object, combine = FALSE, design = TRUE)
      Xpsvd <- svd(Xmats$X$p, nv = 0)
      Xmats$X$s <- Xmats$X$s - Xpsvd$u %*% crossprod(Xpsvd$u, Xmats$X$s)
      colnames(Xmats$X$s) <- snames
      nsdf <- object$nsdf
      for(i in 1:nterms){
        sindx <- which(snames == object$terms[i])
        X <- Xmats$X$s[,sindx]
        ss[i] <- sum((psolve(X %*% object$cov.sqrt[nsdf + sindx,]) %*% (X %*% s.coef[sindx]))^2) * object$sigma^2
        s.df[i] <- sum(rowSums((X %*% object$cov.sqrt[nsdf + sindx,] / object$sigma)^2))
      }
      
    } # end if(object$specs$tprk)
    
    # smooth coefficients table
    s.Fstat <- (pmax(ss, 0) / s.df) / object$sigma^2
    s.pvals <- 1 - pf(s.Fstat, df1 = s.df, df2 = erdf)
    s.table <- rbind(cbind(s.df, ss, ss / s.df, s.Fstat, s.pvals),
                     c(erdf, erdf * object$sigma^2, object$sigma^2, NA, NA))
    colnames(s.table) <- c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")
    rownames(s.table) <- c(paste0("s(",object$terms, ")"), "Residuals")
    
    # overall F stat
    numdf <- object$df - 1L
    demdf <- nobs - object$df
    coef <- object$coefficients[-1]
    chisq <- sum((psolve(object$cov.sqrt[-1,]) %*% coef)^2)
    Fstat <- chisq / numdf
    fstatistic <- c(Fstat, numdf, demdf)
    names(fstatistic) <- c("value", "numdf", "demdf")
    
    # return results
    res <- list(residuals = resid, 
                fstatistic = fstatistic,
                p.table = as.data.frame(p.table),
                s.table = as.data.frame(s.table),
                sigma = object$sigma,
                r.squared = object$r.squared,
                adj.r.squared = adjRsq, 
                kappa = kappa, pi = pi, 
                call = object$call)
    class(res) <- "summary.sm"
    return(res)
    
  } # end summary.sm

print.summary.sm <-
  function(x, digits = max(3, getOption("digits") - 3), 
           signif.stars = getOption("show.signif.stars"), ...){
    
    # set digits
    oldoptions <- options(digits = digits)
    on.exit(options(oldoptions))
    
    # print call
    cat("\nCall:\n")
    print(x$call)
    
    # print residuals
    cat("\nResiduals:\n")
    rsum <- summary(x$residuals)[-4]
    names(rsum) <- c("Min", "1Q", "Median", "3Q", "Max")
    print(rsum)
    
    # check if significance codes are needed
    if(signif.stars){
      pcodes <- rev(c("***", "**", "*", ".", ""))
    }
    
    # parametric effect table
    cat("\nApprox. Signif. of Parametric Effects:\n")
    coef <- as.data.frame(x$p.table)
    if(signif.stars){
      pstars <- symnum(x$p.table[,4],
                       cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                       symbols = c("***", "**", "*", ".", " "))
      coefnames <- names(coef)
      coef <- cbind(coef, pstars)
      names(coef) <- c(coefnames, "")
    }
    #coef[,4] <- ifelse(coef[,4] < 2e-16, "<2e-16", round(coef[,4], digits))
    print(coef)
    if(signif.stars){
      cat("---\n")
      cat("Signif. codes: ", attr(pstars, "legend"), "\n")
    }
    
    # smooth effect table
    cat("\nApprox. Signif. of Nonparametric Effects:\n")
    smoo <- as.data.frame(x$s.table)
    if(signif.stars){
      pstars <- symnum(x$s.table[,5],
                       cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                       symbols = c("***", "**", "*", ".", " "), na = NA)
      smoonames <- names(smoo)
      smoo <- cbind(smoo, pstars)
      colnames(smoo) <- c(smoonames, "")
    }
    #smoo[,5] <- ifelse(smoo[,5] < 2e-16, "<2e-16", round(smoo[,5], digits))
    print(as.matrix(smoo), na.print = "", quote = FALSE, right = TRUE)
    if(signif.stars){
      cat("---\n")
      cat("Signif. codes: ", attr(pstars, "legend"), "\n")
    }
    
    cat("\nResidual standard error:", round(x$sigma, 3), "on", round(x$fstatistic[3],1 ), "degrees of freedom")
    cat("\nMultiple R-squared:  ", round(x$r.squared,4), 
        ",    Adjusted R-squared:  ", round(x$adj.r.squared, 4), sep = "")
    pval <- 1 - pf(x$fstatistic[1], x$fstatistic[2], x$fstatistic[3])
    if(pval < 2e-16) {
      pval <- "<2e-16"
    } else if(pval < 0.0001) {
      pval <- "<.0001"
    } 
    cat("\nF-statistic:", x$fstatistic[1], "on", x$fstatistic[2], "and", x$fstatistic[3], "DF,  p-value:", pval, "\n\n")
    
  } # end print.summary.sm
