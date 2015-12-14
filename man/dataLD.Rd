% File binreg/man/dataLD.Rd
% Part of the binreg package for R (http://www.R-project.org)
% Copyright (C) 2013 Bernardo dos Santos, Adriano Polpo, 
%                    Carlos A. de B. Pereira.
% Distributed under GPL 3 or later

\name{dataLD}

%%%
\alias{dataLD}

\title{Transform formula to LaplacesDemon data.}
\description{
  Transform a formula of binary regression model in data to be used with LaplacesDemon function.
}
\usage{
 dataLD(formula, link, data, subset)
}
\arguments{
  \item{formula}{an object of class 'formula' (or one that can be coerced to
                 that class): a symbolic description of the model to be
                 fitted.}
  \item{link}{the link function to be used. The options are: 'logit', 'probit'
              'cloglog', 'Aranda-Ordaz', 'Weibull', 'CWeibull', 'Stukel', 'Prentice'.}
  \item{data}{an optional data frame, list or environment (or object
              coercible by 'as.data.frame' to a data frame) containing the
              variables in the model. If not found in 'data', the
              variables are taken from 'environment(formula)'.}
  \item{subset}{an optional vector specifying a subset of observations to be
                used in the fitting process.}
}
\details{
  It is a function to facilitate de use of LaplacesDemon function. 
  Transform a formula of binary regression model in data to be used with 
  LaplacesDemon function, and suggest an command to run the LaplacesDemon.
}
\value{
\item{J}{number of covariates.}
\item{mon.names}{name of the variable to monitored.}
\item{y}{binary response variable.}
\item{X}{covariates of the model.}
\item{parm.names}{name of the parameters.}
}
\references{
  Statisticat, LLC. (2013). \emph{LaplacesDemon: Complete Environment for
  Bayesian Inference}. CRAN. R package version 13.03.04.
  http://cran.r-project.org/web/packages/LaplacesDemon/index.html
}
\examples{
  set.seed(666)
  x <- rnorm(0,1)
  y <- rbern(length(x),pnorm(0.5+2*x))
  dataLD(y ~ x, "cloglog")
}
