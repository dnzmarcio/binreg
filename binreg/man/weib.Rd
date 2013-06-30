% File binreg/man/weib.Rd
% Part of the binreg package for R (http://www.R-project.org)
% Copyright (C) 2013 Bernardo dos Santos, Adriano Polpo, 
%                    Carlos A. de B. Pereira.
% Distributed under GPL 3 or later

\name{weib}

%%%
\alias{dweib}
\alias{pweib}
\alias{qweib}
\alias{rweib}

\title{The Weibull distribution, for binary regression}
\description{
  Density, distribution function, quantile function and random generation 
  function for the Weibull distribution with one parameter.
}
\usage{
 dweib(x, gamma)
 pweib(q, gamma)
 qweib(p, gamma)
 rweib(n, gamma)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.}
  \item{gamma}{shape parameter.}
}
\details{
  The Weibull distribution with parameter gamma has density given by

  f(x) = gamma*(x^(gamma - 1))*exp(-x^gamma)

  for x > 0 and gamma > 0.
}
\value{
  `dweib' gives the density, `pweib' gives the distribution 
  function, `qweib' gives the quantile function, and `rweib'
  generates random deviates.
}
\references{
  Caron, R. and Polpo, A. (2009) \emph{Binary data regression: Weibull 
  distribution}. In: Aip Conference Proceedings, volume 1193, ed. 
  Goggans, P. M. and Chan, C. Y.
}
\examples{
  set.seed(1234)
  gamma <- 2
  x <- seq(0.01,5,0.01)
  y <- pweib(x,gamma)
  plot(x,y,type="l",xlab="x",ylab="pweib")
}
