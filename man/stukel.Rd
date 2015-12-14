% File binreg/man/stukel.Rd
% Part of the binreg package for R (http://www.R-project.org)
% Copyright (C) 2013 Bernardo dos Santos, Adriano Polpo, 
%                    Carlos A. de B. Pereira.
% Distributed under GPL 3 or later

\name{pstukel}

%%%
\alias{dstukel}
\alias{pstukel}
\alias{qstukel}
\alias{rstukel}

\title{The Stukel distribution}
\description{
  Density, distribution function, quantile function and random generation 
  function for the Prenctice distribution.
}
\usage{
 dstukel(x, alpha)
 pstukel(q, alpha)
 qstukel(p, alpha)
 rstukel(n, alpha)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.}
  \item{alpha}{vector of size two with shape parameters.}
}
\details{
  The Stukel distribution with parameter m has distribution function given by

  F(x) = 1/(1 + exp(-h(x,alpha)))

  for alpha[1] > 0 and alpha[2] > 0, h(x,alpha) is

  h(x,alpha) = 
        (exp(alpha[1]*x) - 1)/alpha[1]   if (x >= 0) and (alpha[1] > 0)
        x                                if (x >= 0) and (alpha[1] = 0)
        -log(1 - alpha[1]*x)/alpha[1]    if (x >= 0) and (alpha[1] < 0)
        -(exp(-alpha[2]*x) - 1)/alpha[2] if (x < 0) and (alpha[2] > 0)
        x                                if (x < 0) and (alpha[2] = 0)
        log(1 + alpha[2]*eta)/alpha[2]   if (x < 0) and (alpha[2] < 0)
}
\value{
  'dstukel' gives the density, 'pstukel' gives the distribution 
  function, 'qstukel' gives the quantile function, and 'rstukel'
  generates random deviates.
}
\references{
  Stukel, T. A. (1988) \emph{Generalized logistic models}. Journal of 
  the American Statistical Association, pp. 426-431.
}
\examples{
  set.seed(1234)
  alpha <- c(2,0.5)
  x <- seq(-5,5,0.01)
  y <- pstukel(x,alpha)
  plot(x,y,type="l",xlab="x",ylab="pstukel")
}
