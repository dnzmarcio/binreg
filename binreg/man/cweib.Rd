% File binreg/man/cweib.Rd
% Part of the binreg package for R (http://www.R-project.org)
% Copyright (C) 2013 Bernardo dos Santos, Adriano Polpo, 
%                    Carlos A. de B. Pereira.
% Distributed under GPL 3 or later

\name{cweib}

%%%
\alias{dcweib}
\alias{pcweib}
\alias{qcweib}
\alias{rcweib}

\title{The Reflected Weibull distribution, for binary regression}
\description{
  Density, distribution function, quantile function and random generation 
  function for the Reflected Weibull distribution with one parameter.
}
\usage{
 dcweib(x, gamma)
 pcweib(q, gamma)
 qcweib(p, gamma)
 rcweib(n, gamma)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.}
  \item{gamma}{shape parameter.}
}
\details{
  The Reflected Weibull distribution with parameter gamma has density given by

  f(x) = gamma*((-x)^(gamma - 1))*exp(-x^gamma)

  for x < 0 and gamma > 0.
}
\value{
  `dcweib' gives the density, `pcweib' gives the distribution 
  function, `qcweib' gives the quantile function, and `rcweib'
  generates random deviates.
}
\references{
  Rinne, H. (2009). \emph{The Weibull distribution : a handbook}. Boca Raton, CRC
  Press.
}
\examples{
  set.seed(1234)
  gamma <- 2
  x <- seq(-5,-0.01,0.01)
  y <- pcweib(x,gamma)
  plot(x,y,type="l",xlab="x",ylab="pcweib")
}
