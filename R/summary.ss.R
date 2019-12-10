summary.ss <- 
  function(object, ...){
    # summary method for class "ss"
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2019-12-08
    
    # get residuals
    if(is.null(object$data)){
      fitted <- resid <- NA
    } else {
      fitted <- predict(object, x = object$data$x)$y
      resid <- object$data$y - fitted
    }
    
    # adjusted R-squared
    object$r.squared <- adjRsq <- NA
    if(!is.null(object$data)) {
      object$r.squared <- cor(object$data$y, fitted)^2
      adjRsq <- 1 - object$sigma^2 / var(object$data$y)
    }
    
    # some basic info
    nobs <- object$fit$n
    if(object$fit$periodic){
      nullindx <- 1L
    } else {
      nullindx <- 1:object$fit$m
    }
    erdf <- nobs - object$df
    
    # parametric coefficients table
    p.coef <- object$fit$coef[nullindx]
    p.ster <- sqrt(rowSums(object$fit$cov.sqrt[nullindx,,drop=FALSE]^2))
    p.tval <- p.coef / p.ster
    p.pval <- 2 * (1 - pt(abs(p.tval), df = erdf))
    p.table <- data.frame(p.coef, p.ster, p.tval, p.pval)
    rownames(p.table) <- names(p.coef)
    colnames(p.table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    
    # smooth term degrees of freedom and SS
    s.coef <- object$fit$coef[-nullindx]
    ss <- sum((psolve(object$fit$cov.sqrt[-nullindx,]) %*% s.coef)^2) * object$sigma^2
    s.df <- object$df - object$fit$m
    
    # smooth coefficients table
    s.Fstat <- (ss / s.df) / object$sigma^2
    s.pvals <- 1 - pf(s.Fstat, df1 = s.df, df2 = erdf)
    s.table <- rbind(c(s.df, ss, ss / s.df, s.Fstat, s.pvals),
                     c(erdf, erdf * object$sigma^2, object$sigma^2, NA, NA))
    colnames(s.table) <- c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")
    rownames(s.table) <- c("s(x)", "Residuals")
    
    # overall F stat
    numdf <- object$df - 1L
    demdf <- erdf
    coef <- object$fit$coef[-1]
    chisq <- sum((psolve(object$fit$cov.sqrt[-1,]) %*% coef)^2)
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
                call = object$call)
    class(res) <- "summary.ss"
    return(res)
    
  } # end summary.ss


print.summary.ss <-
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
    
    # fit information
    cat("\nResidual standard error:", x$sigma, "on", x$fstatistic[3], "degrees of freedom\n")
    cat("Multiple R-squared:  ", x$r.squared, ",    Adjusted R-squared:  ", x$adj.r.squared,"\n", sep = "")
    pvalue <- pf(x$fstatistic[1], df1 = x$fstatistic[2], df2 = x$fstatistic[3], lower.tail = FALSE)
    cat("F-statistic:", x$fstatistic[1], "on", x$fstatistic[2], "and", x$fstatistic[3], "DF,  p-value:",
        ifelse(pvalue < 2e-16, "<2e-16", pvalue), "\n\n")
    
  } # end print.summary.ss
