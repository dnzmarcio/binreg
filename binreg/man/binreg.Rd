% File binreg/man/binreg.Rd
% Part of the binreg package for R (http://www.R-project.org)
% Copyright (C) 2013 Bernardo dos Santos, Adriano Polpo, 
%                    Carlos A. de B. Pereira.
% Distributed under GPL 3 or later

\name{binreg}

%%%
\alias{binreg}

\title{Binary regression models.}
\description{
  Performs the Maximum Likelihood Estimation (MLE) of some models for binary data.
}
\usage{
 binreg(formula, link, data, subset, start, tol, iterlim, na.action)
}
\arguments{
  \item{formula}{an object of class 'formula' (or one that can be coerced to
                 that class): a symbolic description of the model to be
                 fitted.}
  \item{link}{the link function to be used. The options are: 'Aranda-Ordaz',
              'Weibull', 'CWeibull', 'Stukel', 'Prentice'.}
  \item{data}{an optional data frame, list or environment (or object
              coercible by 'as.data.frame' to a data frame) containing the
              variables in the model. If not found in 'data', the
              variables are taken from 'environment(formula)'.}
  \item{subset}{an optional vector specifying a subset of observations to be
                used in the fitting process.}
  \item{start}{initial guess to model's parameters. Generally it is not necessary.}
  \item{tol}{tolerance error used in the estimation.}
  \item{iterlim}{maximum number of iteration to be performed.}
  \item{na.action}{a function which indicates what should happen when the data
                   contain 'NA's. More details see 'glm' function.}
}
\details{
  It is a function to facilitate de use of LaplacesDemon function. 
  Transform a formula of binary regression model in data to be used with 
  LaplacesDemon function, and suggest an command to run the LaplacesDemon.
}
\value{
  'binreg' returns an object of class inheriting from 'glm'. For more details
  see 'glm' function.
}
\seealso{
\code{\link{glm}}
}
\examples{
  set.seed(666)
  x <- rnorm(0,1)
  y <- rbern(length(x),pnorm(0.5+2*x))
  fit <- binreg(y ~ x, "Aranda")
  summary(fit)
}
