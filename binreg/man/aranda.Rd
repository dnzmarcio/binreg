% File binreg/man/aranda.Rd
% Part of the binreg package for R (http://www.R-project.org)
% Copyright (C) 2013 Bernardo dos Santos, Adriano Polpo, 
%                    Carlos A. de B. Pereira.
% Distributed under GPL 3 or later

\name{aranda}

%%%
\alias{daranda}
\alias{paranda}
\alias{qaranda}
\alias{raranda}

\title{The Aranda distribution}
\description{
  Density, distribution function, quantile function and random generation 
  function for the Aranda distribution.
}
\usage{
 daranda(x, alpha)
 paranda(q, alpha)
 qaranda(p, alpha)
 raranda(n, alpha)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.}
  \item{alpha}{shape parameter.}
}
\details{
  The Aranda distribution with parameter alpha has density given by

  f(x) = exp(x)*((alpha*exp(x) + 1)^(-1/alpha - 1))

  for alpha > 0.
}
\value{
  `daranda' gives the density, `paranda' gives the distribution 
  function, `qaranda' gives the quantile function, and `raranda'
  generates random deviates.
}
\references{
  Aranda-Ordaz, F. J. (1981) \emph{On two families of transformations to 
  additivity for binary response data}. Biometrika, 68(2):357–363.
}
\examples{
  set.seed(1234)
  alpha <- 2
  x <- seq(-5,5,0.01)
  y <- paranda(x,alpha)
  plot(x,y,type="l",xlab="x",ylab="paranda")
}
