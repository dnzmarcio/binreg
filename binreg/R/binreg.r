require(LaplacesDemon)

#############################################
############ Auxiliary functions ############
#############################################

pprentice <- function(q, m = c(1, 1)) {
  if(length(m) < 2)
    stop("m = ", deparse(substitute(m)), " should have exactly two elements.")
  if(length(m) > 2)
    warning("Length of m = ", deparse(substitute(m)), "is greater than 2. Only the first two elements will be used.")
  
  return( pbeta(plogis(q), m[1], m[2]) )
}

qprentice <- function(p, m = c(1, 1)) {
  if(length(m) < 2)
    stop("m = ", deparse(substitute(m)), " should have exactly two elements.")
  if(length(m) > 2)
    warning("Length of m = ", deparse(substitute(m)), "is greater than 2. Only the first two elements will be used.")
  
  return( qlogis(qbeta(p, m[1], m[2])) )
}

dprentice <- function(x, m = c(1, 1)) {
  if(length(m) < 2)
    stop("m = ", deparse(substitute(m)), " should have exactly two elements.")
  if(length(m) > 2)
    warning("Length of m = ", deparse(substitute(m)), "is greater than 2. Only the first two elements will be used.")
  
  exp(x*m[1])*((1 + exp(x))^(-(m[1] + m[2])))/beta(m[1], m[2])
}


hstukel <- function(eta, alpha){
  if(length(alpha) < 2)
    stop("alpha = ", deparse(substitute(alpha)), " should have exactly two elements.")
  if(length(alpha) > 2)
    warning("Length of alpha = ", deparse(substitute(alpha)), "is greater than 2. Only the first two elements will be used.")
  
  h <- function(eta){
    if(eta > 0){
      if(alpha[1] > 0) out <- (exp(alpha[1]*eta) - 1)/alpha[1]
      else if(alpha[1] == 0) out <- eta
      else out <- -log(1 - alpha[1]*eta)/alpha[1]
    } else {
      if(alpha[2] > 0) out <- -(exp(-alpha[2]*eta) - 1)/alpha[2]
      else if(alpha[2] == 0) out <- eta
      else out <- log(1 + alpha[2]*eta)/alpha[2]
    }
    return(out)
  }
  
  return(sapply(eta, h))
}


dhstukel <- function(eta, alpha){
  if(length(alpha) < 2)
    stop("alpha = ", deparse(substitute(alpha)), " should have exactly two elements.")
  if(length(alpha) > 2) 
    warning("Length of alpha = ", deparse(substitute(alpha)), "is greater than 2. Only the first two elements will be used")
  
  h <- function(eta){
    if(eta > 0){
      if(alpha[1] > 0) out <- exp(alpha[1]*eta)
      else if(alpha[1] == 0) out <- 1
      else out <- 1/(1 - alpha[1]*eta)
    } else {
      if(alpha[2] > 0) out <- exp(-alpha[2]*eta)
      else if(alpha[2] == 0) out <- 1
      else out <- 1/(1 + alpha[2]*eta)
    }
    return(out)
  }
  
  return(sapply(eta, h))
}


invhstukel <- function(y, alpha){
  if(length(alpha) < 2)
    stop("alpha = ", deparse(substitute(alpha)), " should have exactly two elements.")
  if(length(alpha) > 2) 
    warning("Length of alpha = ", deparse(substitute(alpha)), "is greater than 2. Only the first two elements will be used")
  
  h <- function(y){
    if(y > 0){
      if(alpha[1] > 0) out <- log(y*alpha[1] + 1)/alpha[1]
      else if(alpha[1] == 0) out <- y
      else out <- -(exp(-y*alpha[1]) - 1)/alpha[1]
    } else {
      if(alpha[2] > 0) out <- -log(1 - y*alpha[2])/alpha[2]
      else if(alpha[2] == 0) out <- y
      else out <- (exp(y*alpha[2]) - 1)/alpha[2]
    }
  }
  
  return(sapply(y, h)) 
}

#############################################
########### Likelihood functions ############
#############################################


likaranda <- function(alpha, beta, y, X) {
  if(alpha > 0){
    eta <- X %*% beta
    aranda <- ifelse(alpha*exp(eta) > -1, 1 - (alpha*exp(eta) + 1)^(-1/alpha), 1)
    aranda <- pmin(pmax(aranda, .Machine$double.eps), 1 - .Machine$double.eps)
    return(-sum(dbern(y, aranda, log = TRUE)))
  } else{
    return(1e100)
  }
}


likweibull = function(gama, beta, y, X){
  if (gama > 0) {
    eta <- X %*% beta
    eta <- pmax(eta, 0)
    weibull <- 1 - exp(-eta^gama)
    weibull <- pmin(pmax(weibull, .Machine$double.eps), 1 - .Machine$double.eps)
    return(-sum(dbern(y, weibull, log = TRUE)))
  } else {
    return(1e100)
  }
}


likprentice <- function(m = c(1, 1), beta, y, X) {
  if(all(m > 0)){
    eta <- X %*% beta
    prentice <- pprentice(eta, m)
    prentice <- pmin(pmax(prentice, .Machine$double.eps), 1 - .Machine$double.eps)
    return(-sum(dbern(y, prentice, log = TRUE)))
  } else {
    return(1e100)
  }
}


likstukel <- function(alpha, beta, y, X) {
  eta <- X %*% beta
  stukel <- plogis(hstukel(eta, alpha))
  stukel <- pmin(pmax(stukel, .Machine$double.eps), 1 - .Machine$double.eps)
  return(-sum(dbern(y, stukel, log = TRUE)))
}

#############################################
############### Link functions ##############
#############################################

aranda <- function(alpha = 1) {
  linkfun <- function(mu)
    log( ((1 - mu)^(-alpha) - 1)/alpha )
  linkinv <- function(eta)
    ifelse(alpha*exp(eta) > -1, 1 - (alpha*exp(eta) + 1)^(-1/alpha), 1)
  mu.eta <- function(eta)
    ifelse(alpha*exp(eta) > -1, exp(eta)*(alpha*exp(eta) + 1)^(-1/alpha - 1), 0)
  
  link <- paste("Aranda-Ordaz, alpha =", alpha)
  valideta <- function(eta) TRUE
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta, name = link), 
            class = "link-glm")
}


weibull <- function(gama){
  linkfun <- function(mu)
    return((-log(1 - mu))^(1/gama))
  linkinv <- function(eta) {
    eta <- pmax(eta, 0)
    return(1 - exp(-(eta^gama)))
  }
  mu.eta <- function(eta) {
    eta <- pmax(eta, 0) 
    return(gama*eta^(gama - 1)*exp(-eta^gama))
  }
  link <- paste("Weibull, gamma =", gama)
  valideta <- function(eta)
    return(all(eta > 0))
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta, name = link), 
            class = "link-glm")
}


prentice <- function(m = c(1, 1)) {
  linkfun <- function(mu)  qprentice(mu,m)
  linkinv <- function(eta) pprentice(eta,m)
  mu.eta  <- function(eta) dprentice(eta,m)
  link     <- paste("Prentice, m1 = ", m[1],"  -- m2 = ", m[2], sep = "")
  valideta <- function(eta) TRUE
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta, name = link), 
            class = "link-glm")
}


stukel <- function(alpha = c(0, 0)) {
  linkfun  <- function(mu)  invhstukel(qlogis(mu), alpha)
  linkinv  <- function(eta) plogis(hstukel(eta, alpha))
  mu.eta   <- function(eta) dlogis(hstukel(eta, alpha))*dhstukel(eta, alpha)
  link     <- paste("Stukel, a1 = ", alpha[1],"  -- a2 = ",alpha[2], sep = "")
  valideta <- function(eta) TRUE
  structure(list(linkfun = linkfun, linkinv = linkinv, 
                 mu.eta = mu.eta, valideta = valideta, name = link), 
            class = "link-glm")
}



############################################
################# Wrapper ##################
############################################

binreg <- function(formula, link = c("Aranda-Ordaz", "Weibull", "Prentice", "Stukel"), 
                   data, subset, start = NULL, tol = 1e-4, iterlim = 5000, na.action){
  
  call <- match.call()
  if(class(formula) != "formula")
    stop("Invalid formula")
  if(missing(data)) 
    data <- environment(formula)
  if(!missing(subset))
    data <- subset(data, subset)
  X    <- model.matrix(formula, data)
  resp <- model.frame(formula, data)[,1]
  resp <- factor(resp)
  
  if(nlevels(resp) > 2)
    warning("Number of levels greater than 2. The first one will be regressed against the others")
  if(nlevels(resp) < 2)
    stop("Response has fewer than two levels")
  resp <- ifelse(resp == levels(resp)[1], 0, 1)
  
  link <- match.arg(link)
  i    <- 0
  aic  <- 0
  
  if(link == "Aranda-Ordaz"){
    a    <- 1
    fit  <- glm(formula, binomial(aranda(a)), data, start = start)
    bet  <- coef(fit)
    cand <- nlm(likaranda, a, beta = bet, y = resp, X = X)
    
    while (((abs(cand$estimate - a) > tol) || (abs((cand$estimate - a)/a) > tol)) && (abs(fit$aic - aic) > tol)){
      aic  <- fit$aic
      a    <- cand$estimate
      fit  <- glm(formula, binomial(aranda(a)), data, start = bet)
      bet  <- coef(fit)
      cand <- nlm(likaranda, a, beta = bet, y = resp, X = X, hessian = TRUE)
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
  
  if(link == "Weibull"){
    gama <- 3.6
    fit  <- glm(formula, binomial(weibull(gama)), data, start = start)
    bet  <- coef(fit)
    cand <- nlm(likweibull, gama, beta = bet, y = resp, X = X)
    flag <- 0.1
    
    while(((abs(cand$estimate - gama) > tol) || (abs((cand$estimate - gama)/gama) > tol)) && (abs(fit$aic - aic) > tol)){
      if ((i >= 10) && (abs(fit$aic - aic) < tol) &&
        (all(flag[(length(flag)-10):length(flag)] > 0))){
        warning("Link parameter is too large. Complementary log-log estimate returned.")
        fit <- glm(formula, binomial("cloglog"), data, na.action = na.action)
        fit$formula <- formula
        fit$call <- call
        return(fit)
      }
      
      aic     <- fit$aic
      gama    <- cand$estimate
      fit     <- glm(formula, binomial(weibull(gama)), data, start = bet)
      bet     <- coef(fit)
      cand    <- nlm(likweibull, gama, beta = bet, y = resp, X = X, hessian = TRUE)
      i       <- i + 1
      flag[i] <- (cand$estimate - gama)
      
      if(i == iterlim){
        warning("Maximum iteration reached. Solution may not be optimal.")
        break
      }
    }
  
  gama <- cand$estimate
  fit <- glm(formula, binomial(weibull(gama)), data, start = bet, na.action = na.action)
  fit$link.par <- gama
  fit$link.err <- as.numeric(sqrt(1/cand$hessian))
  names(fit$link.par) <- names(fit$link.err) <- "gamma"
  fit$aic <- fit$aic + 2
  }
  
  if(link == "Prentice"){
    m    <- c(1, 1)
    fit  <- glm(formula, binomial(prentice(m)), data, start = start)
    bet  <- coef(fit)
    cand <- nlm(likprentice, m, beta = bet, y = resp, X = X)
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
      cand <- nlm(likprentice, m, beta = bet, y = resp, X = X, hessian = TRUE)
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
    cand <- nlm(likstukel, a, beta = bet, y = resp, X = X)
    
    while(((any(abs(cand$estimate - a) > tol)) || (any(abs((cand$estimate - a)/a) > tol))) && (abs(fit$aic - aic) > tol)){
      aic  <- fit$aic
      a    <- cand$estimate
      fit  <- glm(formula, binomial(stukel(a)), data, start = bet)
      bet  <- coef(fit)
      cand <- nlm(likstukel, a, beta = bet, y = resp, X = X, hessian = TRUE)
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
  p <- object$rank + switch(link, A = 1, W = 1, P = 2, S = 2)
  
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

compareIC <- function(...){
  names <- as.character(match.call())[-1]
  
  objects <- list(...)
  
  df <- sapply(objects, function(x) attr(logLik(x), "df"))
  llik <- sapply(objects, logLik)
  aic <- sapply(objects, AIC)
  bic <- sapply(objects, BIC)
  
  data.frame(df = df, logLik = llik, AIC = aic, BIC = bic, row.names = names)
}

LRTlogit <- function(object){
  if(!inherits(object, "binreg"))
    stop("The object is not from class binreg")
  
  dname <- deparse(substitute(object))
  logit <- glm(object$formula, binomial, object$data)
  lr <- 2*as.numeric(logLik(object) - logLik(logit))
  df <- attr(logLik(object), "df") - attr(logLik(logit), "df")
  pval <- pchisq(lr, df, lower.tail = FALSE)
  method <- "Likelihood ratio test for Binary Regression with parametric link"
  
  RVAL <- list(statistic = c(LR = lr), parameter = c(df = df), p.value = pval,
               method = method, data.name = dname)
  class(RVAL) <- "htest"
  return(RVAL)
}



dataLD <- function(formula, link = c("logit", "probit", "cloglog", "Aranda-Ordaz", "Weibull", "Prentice", "Stukel"), data, subset){
  formula <- formula(formula)
  if(missing(data)) 
    data <- environment(formula)
  if(!missing(subset))
    data <- subset(data, subset)
  X    <- model.matrix(formula, data)
  resp <- model.frame(formula, data)[,1]
  resp <- factor(resp)
  
  if(nlevels(resp) > 2)
    warning("Number of levels greater than 2. The first one will be regressed against the others.")
  if(nlevels(resp) < 2)
    stop("Response has fewer than two levels.")
  resp <- ifelse(resp == levels(resp)[1], 0, 1)
  
  link <- match.arg(link)
  
  if(link == "logit"){
    parm.names <- as.parm.names(list(beta = rep(0, ncol(X))))
    model <- "ModelLogit()"
  }
  
  if(link == "probit"){
    parm.names <- as.parm.names(list(beta = rep(0, ncol(X))))
    model <- "ModelProbit()"
  }
  
  if(link == "cloglog"){
    parm.names <- as.parm.names(list(beta = rep(0, ncol(X))))
    model <- "ModelCloglog()"
  }
  
  if(link == "Aranda-Ordaz"){
    parm.names <- as.parm.names(list(beta = rep(0, ncol(X)), alpha = 0))
    model <- "ModelAranda()"
  }
  
  if(link == "Weibull"){
    parm.names <- as.parm.names(list(beta = rep(0, ncol(X)), gamma = 0))
    model <- "ModelWeibull()"
  }
  
  if(link == "Prentice"){
    parm.names <- as.parm.names(list(beta = rep(0, ncol(X)), m = c(0, 0)))
    model <- "ModelPrentice()"
  }
  
  if(link == "Stukel"){
    parm.names <- as.parm.names(list(beta = rep(0, ncol(X)), a = c(0, 0)))
    model <- "ModelStukel()"
  }
  
  dataset <- list(J = ncol(X), X = X, mon.names = "LP", parm.names = parm.names, y = resp)
  
  cat("Demonic dataset is ready!\n\n")
  print(data.frame(variable = colnames(X), parameter = parm.names[1:ncol(X)]), row.names = F)
  cat("\nLaplace's Demon suggests running the following code:\n\n")
  cat("fit = LaplacesDemon(", model, ", data, GIV(", model, ", data))\n\n", sep = "")
  cat("Remember to substitute 'data' with the object returned by this function!")
  
  invisible(dataset)
}


##############################################
########## Bayesian Model functions ##########
##############################################

ModelLogit <- function(parm, data, beta.mean = 0, beta.var = 1000){
  Model <- function(parm, data){
    ###parameter
    beta <- parm[1:data$J]
    ###priors
    beta.prior <- dnormv(beta, beta.mean, beta.var, log = TRUE)
    ###logLik
    eta <- data$X %*% beta
    mu <- make.link("logit")$linkinv(eta)
    LL <- sum(dbern(data$y, mu, log = TRUE))
    ###log-posterior
    LP <- LL + sum(beta.prior)
    Modelout <- list(LP = LP, Dev = -2*LL, Monitor = LP, mu = eta, parm = parm)
    return(Modelout)
  }
  return(Model)
}


ModelProbit <- function(parm, data, beta.mean = 0, beta.var = 1000){
  Model <- function(parm, data){
    ###parameter
    beta <- parm[1:data$J]
    ###priors
    beta.prior <- dnormv(beta, beta.mean, beta.var, log = TRUE)
    ###logLik
    eta <- data$X %*% beta
    mu <- make.link("probit")$linkinv(eta)
    LL <- sum(dbern(data$y, mu, log = TRUE))
    ###log-posterior
    LP <- LL + sum(beta.prior)
    Modelout <- list(LP = LP, Dev = -2*LL, Monitor = LP, mu = eta, parm = parm)
    return(Modelout)
  }
  return(Model)
}


ModelCloglog <- function(parm, data, beta.mean = 0, beta.var = 1000){
  Model <- function(parm, data){
    ###parameter
    beta <- parm[1:data$J]
    ###priors
    beta.prior <- dnormv(beta, beta.mean, beta.var, log = TRUE)
    ###logLik
    eta <- data$X %*% beta
    mu <- make.link("cloglog")$linkinv(eta)
    LL <- sum(dbern(data$y, mu, log = TRUE))
    ###log-posterior
    LP <- LL + sum(beta.prior)
    Modelout <- list(LP = LP, Dev = -2*LL, Monitor = LP, mu = eta, parm = parm)
    return(Modelout)
  }
  return(Model)
}


ModelAranda <- function(parm, data, beta.mean = 0, beta.var = 1000,
                        alpha.scale = 1e-10){
   Model <- function(parm, data){
    ###parameter
    beta <- parm[1:data$J]
    alpha <- parm[data$J + 1]
    ###priors
    beta.prior <- dnormv(beta, beta.mean, beta.var, log = TRUE)
    alpha.prior <- dhalfnorm(alpha, alpha.scale, log = TRUE)
    ###logLik
    LL <- -likaranda(alpha, beta, data$y, data$X)
    ###log-posterior
    LP <- LL + sum(beta.prior) + alpha.prior 
    Modelout <- list(LP = LP, Dev = -2*LL, Monitor = LP, mu = data$X %*% beta, parm = parm)
    return(Modelout)
   }
   return(Model)
}



ModelWeibull <- function(parm, data, beta.mean = 0, beta.var = 1000, 
                        gamma.scale = 1e-10){
  Model <- function(parm, data){
    ###parameter
    beta <- parm[1:data$J]
    gamma <- parm[data$J + 1]
    ###priors
    beta.prior <- dnormv(beta, beta.mean, beta.var, log = TRUE)
    gamma.prior <- dhalfnorm(gamma, gamma.scale, log = TRUE)
    ###logLik
    LL <- -likweibull(gamma, beta, data$y, data$X)
    ###log-posterior
    LP <- LL + sum(beta.prior) + gamma.prior 
    Modelout <- list(LP = LP, Dev = -2*LL, Monitor = LP, mu = data$X %*% beta, parm = parm)
    return(Modelout)
  }
  return(Model)
}


ModelPrentice <- function(parm, data, beta.mean = 0, beta.var = 1000, 
                          m.scale = 1e-10){
  Model <- function(parm, data){
    ###parameter
    beta <- parm[1:data$J]
    m <- parm[c(data$J + 1, data$J + 2)]
    ###priors
    beta.prior <- dnormv(beta, beta.mean, beta.var, log = TRUE)
    m.prior <- dhalfnorm(m, m.scale, log = TRUE)
    ###logLik
    LL <- -likprentice(m, beta, data$y, data$X)
    ###log-posterior
    LP <- LL + sum(beta.prior) + sum(m.prior)
    Modelout <- list(LP = LP, Dev = -2*LL, Monitor = LP, mu = data$X %*% beta, parm = parm)
    return(Modelout)
  }
  return(Model)
}


ModelStukel <- function(parm, data, beta.mean = 0, beta.var = 1000, 
                        alpha.mean = 0, alpha.var = 1000){
  Model <- function(parm, data){
    ###parameter
    beta <- parm[1:data$J]
    alpha <- parm[c(data$J + 1, data$J + 2)]
    ###priors
    beta.prior <- dnormv(beta, beta.mean, beta.var, log = TRUE)
    alpha.prior <- dnormv(alpha, alpha.mean, alpha.var, log = TRUE)
    ###logLik
    LL <- -likstukel(alpha, beta, data$y, data$X)
    ###log-posterior
    LP <- LL + sum(beta.prior) + sum(alpha.prior)
    Modelout <- list(LP = LP, Dev = -2*LL, Monitor = LP, mu = data$X %*% beta, parm = parm)
    return(Modelout)
  }
  return(Model)
}