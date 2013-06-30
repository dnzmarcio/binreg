# binreg package for R (http://www.R-project.org)
# Copyright (C) 2013 Bernardo dos Santos, Adriano Polpo, 
#                    Carlos A. de B. Pereira.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

binreg <- function(formula, link = c("Aranda-Ordaz", "Weibull", "CWeibull", "Prentice", "Stukel"), 
                   data, subset, start = NULL, tol = 1e-4, iterlim = 5000, na.action){
  
  call <- match.call()
  if(class(formula) != "formula")
    stop("Invalid formula.")
  if(missing(data)) 
    data <- environment(formula)
  if(!missing(subset))
    data <- subset(data, subset)
  X    <- model.matrix(formula, data)
  resp <- model.frame(formula, data)[,1]
  resp <- factor(resp)

  if (length(resp[resp == 0])+length(resp[resp==1]) != length(resp)) {
    stop("Response must be a binary variable.")
  }
  
  if(nlevels(resp) > 2)
    warning("Number of levels greater than 2. The first one will be regressed against the others.")
  if(nlevels(resp) < 2)
    stop("Response has fewer than two levels.")
  resp <- ifelse(resp == levels(resp)[1], 0, 1)
  
  link <- match.arg(link)
  i    <- 0
  aic  <- 0
  
  if(link == "Aranda-Ordaz"){
    a    <- 1
    fit  <- glm(formula, binomial(aranda(a)), data, start = start)
    bet  <- coef(fit)
    cand <- nlm(.lik.aranda, a, beta = bet, y = resp, X = X)
    
    while (((abs(cand$estimate - a) > tol) || (abs((cand$estimate - a)/a) > tol)) && (abs(fit$aic - aic) > tol)){
      aic  <- fit$aic
      a    <- cand$estimate
      fit  <- glm(formula, binomial(aranda(a)), data, start = bet)
      bet  <- coef(fit)
      cand <- nlm(.lik.aranda, a, beta = bet, y = resp, X = X, hessian = TRUE)
      i    <- i + 1
      
      if(i == iterlim){
        warning("Maximum iteration reached. Solution may not be optimal.")
        break
      }
    }
    a <- cand$estimate
    fit <- glm(formula, binomial(aranda(a)), data, start = bet, na.action = na.action)
    fit$link.par <- a
    fit$link.err <- as.numeric(sqrt(1/cand$hessian))
    names(fit$link.par) <- names(fit$link.err) <- "alpha"
    fit$aic <- fit$aic + 2
  }
  
  if(link == "Weibull") {
    gamma <- 3.6
    fit  <- glm(formula, binomial(weibull(gamma)), data, start = start)
    bet  <- coef(fit)
    cand <- nlm(.lik.weibull, gamma, beta = bet, y = resp, X = X)
    flag <- 0.1
    
    while(((abs(cand$estimate - gamma) > tol) || (abs((cand$estimate - gamma)/gamma) > tol)) && (abs(fit$aic - aic) > tol)){
      if ((i >= 10) && (abs(fit$aic - aic) < tol) &&
        (all(flag[(length(flag)-10):length(flag)] > 0))) {
        warning("Link parameter is too large. Complementary log-log estimate returned.")
        fit <- glm(formula, binomial("cloglog"), data, na.action = na.action)
        fit$formula <- formula
        fit$call <- call
        return(fit)
      }
      
      aic     <- fit$aic
      gamma   <- cand$estimate
      fit     <- glm(formula, binomial(weibull(gamma)), data, start = bet)
      bet     <- coef(fit)
      cand    <- nlm(.lik.weibull, gamma, beta = bet, y = resp, X = X, hessian = TRUE)
      i       <- i + 1
      flag[i] <- (cand$estimate - gamma)
      
      if(i == iterlim){
        warning("Maximum iteration reached. Solution may not be optimal.")
        break
      }
    }
  
    gamma <- cand$estimate
    fit <- glm(formula, binomial(weibull(gamma)), data, start = bet, na.action = na.action)
    fit$link.par <- gamma
    fit$link.err <- as.numeric(sqrt(1/cand$hessian))
    names(fit$link.par) <- names(fit$link.err) <- "gamma"
    fit$aic <- fit$aic + 2
  }

  if(link == "CWeibull"){
    gamma <- 3.6
    fit  <- glm(formula, binomial(cweibull(gamma)), data, start = start)
    bet  <- coef(fit)
    cand <- nlm(.lik.cweibull, gamma, beta = bet, y = resp, X = X)
    flag <- 0.1
    
    while(((abs(cand$estimate - gamma) > tol) || (abs((cand$estimate - gamma)/gamma) > tol)) && (abs(fit$aic - aic) > tol)){
      if ((i >= 10) && (abs(fit$aic - aic) < tol) &&
        (all(flag[(length(flag)-10):length(flag)] > 0))){
        warning("Link parameter is too large. Log-log estimate returned.")
        fit <- glm(formula, binomial("loglog"), data, na.action = na.action)
        fit$formula <- formula
        fit$call <- call
        return(fit)
      }
      
      aic     <- fit$aic
      gamma   <- cand$estimate
      fit     <- glm(formula, binomial(cweibull(gamma)), data, start = bet)
      bet     <- coef(fit)
      cand    <- nlm(.lik.cweibull, gamma, beta = bet, y = resp, X = X, hessian = TRUE)
      i       <- i + 1
      flag[i] <- (cand$estimate - gamma)
      
      if(i == iterlim){
        warning("Maximum iteration reached. Solution may not be optimal.")
        break
      }
    }
  
    gamma <- cand$estimate
    fit <- glm(formula, binomial(cweibull(gamma)), data, start = bet, na.action = na.action)
    fit$link.par <- gamma
    fit$link.err <- as.numeric(sqrt(1/cand$hessian))
    names(fit$link.par) <- names(fit$link.err) <- "gamma"
    fit$aic <- fit$aic + 2
  }
  
  if(link == "Prentice"){
    m    <- c(1, 1)
    fit  <- glm(formula, binomial(prentice(m)), data, start = start)
    bet  <- coef(fit)
    cand <- nlm(.lik.prentice, m, beta = bet, y = resp, X = X)
    flag <- c(0.1, 0.1)
    
    while(((any(abs(cand$estimate - m) > tol)) || (any((abs(cand$estimate - m)/m) > tol))) && (abs(fit$aic - aic) > tol)){
      if ((i >= 10) && (abs(fit$aic - aic) < tol) &&
        ((all(flag[(length(flag)-10):length(flag),1] > 0)) &&
        (all(flag[(length(flag)-10):length(flag),2] > 0)))){
        warning("Link parameters are too large. Probit estimate returned.")
        fit <- glm(formula, binomial("probit"), data, na.action = na.action)
        fit$formula <- formula
        fit$call <- call
        return(fit)
      }
      
      aic  <- fit$aic
      m    <- cand$estimate
      fit  <- glm(formula, binomial(prentice(m)), data, start = bet)
      bet  <- coef(fit)
      cand <- nlm(.lik.prentice, m, beta = bet, y = resp, X = X, hessian = TRUE)
      i    <- i + 1
      flag <- rbind(flag, (cand$estimate - m))
      
      if(i == iterlim){
        warning("Maximum iteration reached. Solution may not be optimal.")
        break
      }
    }
    
    m <- cand$estimate
    fit <- glm(formula, binomial(prentice(m)), data, start = bet, na.action = na.action)
    fit$link.par <- m
    fit$link.err <- as.numeric(sqrt(diag(solve(cand$hessian))))
    names(fit$link.par) <- names(fit$link.err) <- c("m1", "m2")
    fit$aic <- fit$aic + 4
  }
  
  if(link == "Stukel"){
    a    <- c(0, 0)
    fit  <- glm(formula, binomial(stukel(a)), data, start = start)
    bet  <- coef(fit)
    cand <- nlm(.lik.stukel, a, beta = bet, y = resp, X = X)
    
    while(((any(abs(cand$estimate - a) > tol)) || (any(abs((cand$estimate - a)/a) > tol))) && (abs(fit$aic - aic) > tol)){
      aic  <- fit$aic
      a    <- cand$estimate
      fit  <- glm(formula, binomial(stukel(a)), data, start = bet)
      bet  <- coef(fit)
      cand <- nlm(.lik.stukel, a, beta = bet, y = resp, X = X, hessian = TRUE)
      i    <- i + 1
      
      if(i == iterlim){
        warning("Maximum iteration reached. Solution may not be optimal.")
        break
      }
    }
    
    a <- cand$estimate
    fit <- glm(formula, binomial(stukel(a)), data, start = bet, na.action = na.action)
    fit$link.par <- a
    fit$link.err <- as.numeric(sqrt(diag(solve(cand$hessian))))
    names(fit$link.par) <- names(fit$link.err) <- c("a1", "a2")
    fit$aic <- fit$aic + 4
  }
  
  fit$iterations <- i
  fit$formula <- formula
  fit$call <- call
  class(fit) <- c("binreg", class(fit))
  fit
}


#############################################
################## Utility ##################
#############################################

logLik.binreg <- function(object, ...){
  if (length(list(...))) 
    warning("extra arguments discarded")
  
  link <- strtrim(object$family$link, 1)
  p <- object$rank + switch(link, A = 1, W = 1, C = 1, P = 2, S = 2)
  
  val <- p - object$aic/2
  attr(val, "nobs") <- sum(!is.na(object$residuals))
  attr(val, "df") <- p
  class(val) <- "logLik"
  val
}


summary.binreg <- function(object, ...){
  ans <- NextMethod("summary")
  link <- strtrim(object$family$link, 1)
  
  if(link == "A"){
    zval <- (object$link.par - 1)/object$link.err
    pval <- pchisq(zval^2, 1, lower.tail = FALSE)
    ans$coefficients <- rbind(ans$coefficients, alpha = c(object$link.par, object$link.err, zval, pval))
  }
  
  if(link == "W"){
    zval <- (object$link.par - 3.50215)/object$link.err
    pval <- pchisq(zval^2, 1, lower.tail = FALSE)
    ans$coefficients <- rbind(ans$coefficients, gamma = c(object$link.par, object$link.err, zval, pval))
  }
  
  if(link == "C"){
    zval <- (object$link.par - 3.50215)/object$link.err
    pval <- pchisq(zval^2, 1, lower.tail = FALSE)
    ans$coefficients <- rbind(ans$coefficients, gamma = c(object$link.par, object$link.err, zval, pval))
  }

  if(link == "P"){
    zval <- (object$link.par - 1)/object$link.err
    pval <- pchisq(zval^2, 1, lower.tail = FALSE)
    ans$coefficients <- rbind(ans$coefficients, m1 = c(object$link.par[1], object$link.err[1], zval[1], pval[1]),
                                                m2 = c(object$link.par[2], object$link.err[2], zval[2], pval[2]))
  }
  
  if(link == "S"){
    zval <- object$link.par/object$link.err
    pval <- pchisq(zval^2, 1, lower.tail = FALSE)
    ans$coefficients <- rbind(ans$coefficients, a1 = c(object$link.par[1], object$link.err[1], zval[1], pval[1]),
                                                a2 = c(object$link.par[2], object$link.err[2], zval[2], pval[2]))
  }
  
  return(ans)
}

print.binreg <- function (x, digits = max(3, getOption("digits") - 3), ...){
  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients")
    if (is.character(co <- x$contrasts)) 
      cat("  [contrasts: ", apply(cbind(names(co), co), 
                                  1L, paste, collapse = "="), "]")
    cat(":\n")
    print.default(format(x$coefficients, digits = digits), 
                  print.gap = 2, quote = FALSE)
  }
  else cat("No coefficients\n\n")
  
  cat("\nLink Parameters:\n")
  print.default(format(x$link.par, digits = digits), print.gap = 2, quote = FALSE)
  
  cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ", x$df.residual, "Residual\n")
  if (nzchar(mess <- naprint(x$na.action))) 
    cat("  (", mess, ")\n", sep = "")
  cat("Null Deviance:\t   ", format(signif(x$null.deviance, digits)), 
      "\nResidual Deviance:", format(signif(x$deviance, digits)), 
      "\tAIC:", format(signif(x$aic, digits)), "\n")
  invisible(x)
}

#compareIC <- function(...){
#  names <- as.character(match.call())[-1]
#  
#  objects <- list(...)
#  
#  df <- sapply(objects, function(x) attr(logLik(x), "df"))
#  llik <- sapply(objects, logLik)
#  aic <- sapply(objects, AIC)
#  bic <- sapply(objects, BIC)
#  
#  data.frame(df = df, logLik = llik, AIC = aic, BIC = bic, row.names = names)
#}

#LRTlogit <- function(object){
#  if(!inherits(object, "binreg"))
#    stop("The object is not from class binreg")
#  
#  dname <- deparse(substitute(object))
#  logit <- glm(object$formula, binomial, object$data)
#  lr <- 2*as.numeric(logLik(object) - logLik(logit))
#  df <- attr(logLik(object), "df") - attr(logLik(logit), "df")
#  pval <- pchisq(lr, df, lower.tail = FALSE)
#  method <- "Likelihood ratio test for Binary Regression with parametric link"
#  
#  RVAL <- list(statistic = c(LR = lr), parameter = c(df = df), p.value = pval,
#               method = method, data.name = dname)
#  class(RVAL) <- "htest"
#  return(RVAL)
#}




