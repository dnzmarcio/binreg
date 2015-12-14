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

dataLD <- function(formula, 
                   link = c("logit", "probit", "cloglog", "Aranda-Ordaz", 
                            "Weibull", "CWeibull", "Prentice", "Stukel"), data, subset){
  formula <- formula(formula)
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
  
  if(link == "CWeibull"){
    parm.names <- as.parm.names(list(beta = rep(0, ncol(X)), gamma = 0))
    model <- "ModelCWeibull()"
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
  cat("Remember to substitute 'data' with the object returned by this function!\n")
  
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
    Modelout <- list(LP = LP, Dev = -2*LL, Monitor = LP, yhat = mu, parm = parm)
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
    Modelout <- list(LP = LP, Dev = -2*LL, Monitor = LP, yhat = mu, parm = parm)
    return(Modelout)
  }
  return(Model)
}

ModelCloglog <- function(parm, data, beta.mean = 0, beta.var = 1000) {
  Model <- function(parm, data) {
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
    Modelout <- list(LP = LP, Dev = -2*LL, Monitor = LP, yhat = mu, parm = parm)
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
    if (alpha > 0.001) {
      alpha.prior <- dhalfnorm(alpha, alpha.scale, log = TRUE)
    } else {
      alpha.prior <- log(0)
    }
    ###logLik
    LL <- -.lik.aranda(alpha, beta, data$y, data$X)
    ###log-posterior
    LP <- LL + sum(beta.prior) + alpha.prior
    mu <- paranda(data$X %*% beta,alpha)
    Modelout <- list(LP = LP, Dev = -2*LL, Monitor = LP, yhat = mu, parm = parm)
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
    if (gamma > 0.001) {
      gamma.prior <- dhalfnorm(gamma, gamma.scale, log = TRUE)
    } else {
      gamma.prior <- log(0)
    }
    ###logLik
    LL <- -.lik.weibull(gamma, beta, data$y, data$X)
    ###log-posterior
    LP <- LL + sum(beta.prior) + gamma.prior 
    mu <- pweib(data$X %*% beta,gamma)
    Modelout <- list(LP = LP, Dev = -2*LL, Monitor = LP, yhat = mu, parm = parm)
    return(Modelout)
  }
  return(Model)
}

ModelCWeibull <- function(parm, data, beta.mean = 0, beta.var = 1000, 
                          gamma.scale = 1e-10){
  Model <- function(parm, data){
    ###parameter
    beta <- parm[1:data$J]
    gamma <- parm[data$J + 1]
    ###priors
    beta.prior <- dnormv(beta, beta.mean, beta.var, log = TRUE)
    if (gamma > 0.001) {
      gamma.prior <- dhalfnorm(gamma, gamma.scale, log = TRUE)
    } else {
      gamma.prior <- log(0)
    }
    ###logLik
    LL <- -.lik.cweibull(gamma, beta, data$y, data$X)
    ###log-posterior
    LP <- LL + sum(beta.prior) + gamma.prior 
    mu <- pcweib(data$X %*% beta,gamma)
    Modelout <- list(LP = LP, Dev = -2*LL, Monitor = LP, yhat = mu, parm = parm)
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
    if (all(m > 0)) {
      m.prior <- dhalfnorm(m, m.scale, log = TRUE)
    } else {
      m.prior <- log(0)
    }
    ###logLik
    LL <- -.lik.prentice(m, beta, data$y, data$X)
    ###log-posterior
    LP <- LL + sum(beta.prior) + sum(m.prior)
    if (all(m > 0)) {
      mu <- pprentice(data$X %*% beta,m)
    } else {
      mu <- rep(0,length(data$X %*% beta))
    }
    Modelout <- list(LP = LP, Dev = -2*LL, Monitor = LP, yhat = mu, parm = parm)
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
    LL <- -.lik.stukel(alpha, beta, data$y, data$X)
    ###log-posterior
    LP <- LL + sum(beta.prior) + sum(alpha.prior)
    mu <- pstukel(data$X %*% beta,alpha)
    Modelout <- list(LP = LP, Dev = -2*LL, Monitor = LP, yhat = mu, parm = parm)
    return(Modelout)
  }
  return(Model)
}

