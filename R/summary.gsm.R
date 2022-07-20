summary.gsm <- 
  function(object, ...){
    # summary method for class "gsm"
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2022-05-25
    
    # convert binomial (if response is a factor)
    if(object$family$family == "binomial" && any(class(object$data[,1])[1] == c("factor", "ordered"))){
      object$data[,1] <- ifelse(object$data[,1] == levels(object$data[,1])[1], 0, 1)
    }
    
    # get deviance residuals
    nobs <- length(object$linear.predictors)
    fitted.values <- object$family$linkinv(object$linear.predictors)
    weighted <- ifelse(any(colnames(object$data) == "(weights)"), TRUE, FALSE)
    wt <- if(weighted) object$data$'(weights)' else 1
    wsqrt <- sqrt(wt)
    rsign <- sign(object$data[,1] - fitted.values)
    dev.resid <- rsign * sqrt(object$family$dev.resids(object$data[,1], fitted.values, wt))
    
    # get diagnostics
    fit.terms <- scale(predict(object, type = "terms"), scale = FALSE)
    nterms <- ncol(fit.terms)
    Cmat <- crossprod(fit.terms) / tcrossprod(sqrt(colSums(fit.terms^2)))
    kappa <- sqrt(diag(psolve(Cmat)))
    fitc <- rowSums(fit.terms)
    pi <- as.numeric(crossprod(fit.terms, fitc) / sum(fitc^2))
    names(kappa) <- names(pi) <- object$terms
    
    # adjusted R-squared
    adjRsq <- 1 - var(object$data[,1] - fitted.values) * (nobs - 1L) / (var(object$data[,1]) * (nobs - object$df))
    
    # deviance explained
    dev.expl <- (object$null.deviance - object$deviance) / object$null.deviance
    
    # some basic info
    nobs <- length(object$linear.predictors)
    nullindx <- 1:object$nsdf
    erdf <- nobs - object$df
    
    # parametric coefficients table
    p.coef <- object$coefficients[nullindx]
    p.ster <- sqrt(rowSums(object$cov.sqrt[nullindx,,drop=FALSE]^2))
    p.tval <- p.coef / p.ster
    if(object$family$family == "gaussian"){
      p.pval <- 2 * (1 - pt(abs(p.tval), df = erdf))
      p.names <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    } else {
      p.pval <- 2 * (1 - pnorm(abs(p.tval)))
      p.names <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    }
    p.table <- data.frame(p.coef, p.ster, p.tval, p.pval)
    rownames(p.table) <- names(p.coef)
    colnames(p.table) <- p.names
    
    # smooth term names
    s.coef <- object$coefficients[-nullindx]
    snames <- names(s.coef)
    mxknot <- max(sapply(object$specs$knots, function(x) nrow(as.matrix(x))))
    for(m in mxknot:1){
      snames <- gsub(paste0(".knot.", m), "", snames, fixed = TRUE)
    }
    
    # smooth term degrees of freedom and SS
    ss <- s.df <- rep(0, nterms)
    mueta <- object$family$mu.eta(object$linear.predictors)
    vsqrt <- as.numeric(sqrt(mueta^2/object$family$variance(fitted.values)))
    if(object$specs$tprk){
      
      Xmats <- predict(object, combine = FALSE, design = TRUE)
      Xmats$X$p <- wsqrt * vsqrt * Xmats$X$p
      Xmats$X$s <- wsqrt * vsqrt * Xmats$X$s
      Xpsvd <- svd(Xmats$X$p, nv = 0)
      Xmats$X$s <- Xmats$X$s - Xpsvd$u %*% crossprod(Xpsvd$u, Xmats$X$s)
      
      ss <- sum((psolve(Xmats$X$s %*% object$cov.sqrt[-nullindx,,drop=FALSE]) %*% (Xmats$X$s %*% s.coef))^2)
      pis <- rep(0, nterms)
      for(i in 1:nterms){
        fitX <- predict(object, terms = object$terms[i],
                        combine = FALSE, design = TRUE)
        fitX$X$s <- wsqrt * vsqrt * fitX$X$s - Xpsvd$u %*% crossprod(Xpsvd$u, wsqrt * vsqrt * fitX$X$s)
        s.df[i] <- sum(rowSums((fitX$X$s %*% object$cov.sqrt[-nullindx,,drop=FALSE])^2)) / object$dispersion
        pis[i] <- sum((fitX$X$s %*% s.coef)^2)
      }
      pis <- pis / sum(pis)
      ss <- ss * pis
      
    } else {
      
      Xmats <- predict(object, combine = FALSE, design = TRUE)
      Xmats$X$p <- wsqrt * vsqrt * Xmats$X$p
      Xmats$X$s <- wsqrt * vsqrt * Xmats$X$s
      Xpsvd <- svd(Xmats$X$p, nv = 0)
      Xmats$X$s <- Xmats$X$s - Xpsvd$u %*% crossprod(Xpsvd$u, Xmats$X$s)
      
      colnames(Xmats$X$s) <- snames
      nsdf <- object$nsdf
      for(i in 1:nterms){
        sindx <- which(snames == object$terms[i])
        X <- Xmats$X$s[,sindx,drop=FALSE]
        s.df[i] <- sum(rowSums((X %*% object$cov.sqrt[nsdf + sindx,,drop=FALSE])^2)) / object$dispersion
        ss[i] <- sum((psolve(X %*% object$cov.sqrt[nsdf + sindx,,drop=FALSE]) %*% (X %*% s.coef[sindx]))^2)
      }
      
    } # end if(object$specs$tprk)
    
    # smooth coefficients table
    fam.Ftest <- c("gaussian", "Gamma", "inverse.gaussian")
    if(object$family$family %in% fam.Ftest){
      ss <- ss * object$dispersion
      s.Fstat <- (pmax(ss, 0) / s.df) / object$dispersion
      s.pvals <- 1 - pf(s.Fstat, df1 = pmax(s.df, 1), df2 = erdf)
      s.table <- rbind(cbind(s.df, ss, ss / s.df, s.Fstat, s.pvals),
                       c(erdf, object$deviance, object$dispersion, NA, NA))
      colnames(s.table) <- c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")
    } else {
      s.pvals <- 1 - pchisq(ss, df = pmax(s.df, 1))
      s.table <- rbind(cbind(s.df, ss, s.pvals),
                       c(erdf, object$deviance, NA))
      colnames(s.table) <- c("Df", "Chi Sq", "p.value")
    }
    rownames(s.table) <- c(paste0("s(",object$terms, ")"), "Residuals")
    
    # overall F stat
    numdf <- object$df - 1L
    demdf <- nobs - object$df
    coef <- object$coefficients[-1]
    chisq <- sum((psolve(object$cov.sqrt[-1,,drop=FALSE]) %*% coef)^2)
    Fstat <- chisq / numdf
    fstatistic <- c(Fstat, numdf, demdf)
    names(fstatistic) <- c("value", "numdf", "demdf")
    
    # return results
    res <- list(residuals = dev.resid, 
                dev.expl = dev.expl,
                p.table = as.data.frame(p.table),
                s.table = as.data.frame(s.table),
                dispersion = object$dispersion,
                r.squared = object$r.squared,
                adj.r.squared = as.numeric(adjRsq), 
                kappa = kappa, pi = pi, 
                call = object$call,
                family = object$family)
    class(res) <- "summary.gsm"
    return(res)
    
  } # end summary.gsm

print.summary.gsm <-
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
    fam.Ftest <- c("gaussian", "Gamma", "inverse.gaussian")
    if(signif.stars){
      pstars <- symnum(x$s.table[,ifelse(x$family$family %in% fam.Ftest, 5, 3)],
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
    
    cat("\nDispersion Parameter: ", x$dispersion, 
        ",   Deviance Explained: ", x$dev.expl, sep = "")
    cat("\n  Multiple R-squared: ", x$r.squared, 
        ",   Adjusted R-squared: ", x$adj.r.squared, "\n\n",sep = "")
    
  } # end print.summary.gsm
