% File binreg/man/prentice.Rd
% Part of the binreg package for R (http://www.R-project.org)
% Copyright (C) 2013 Bernardo dos Santos, Adriano Polpo, 
%                    Carlos A. de B. Pereira.
% Distributed under GPL 3 or later

\name{prentice}

%%%
\alias{dprentice}
\alias{pprentice}
\alias{qprentice}
\alias{rprentice}

\title{The Prentice distribution}
\description{
  Density, distribution function, quantile function and random generation 
  function for the Prenctice distribution.
}
\usage{
 dprentice(x, m)
 pprentice(q, m)
 qprentice(p, m)
 rprentice(n, m)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.}
  \item{m}{vector of size two with shape parameters.}
}
\details{
  The Prentice distribution with parameter m has density given by

  f(x) = exp(x*m[1])*((1+exp(x))^(-(m[1]+m[2])))/B(m[1],m[2])

  for m[1] > 0 and m[2] > 0, where B(a,b) is the Beta function.
}
\value{
     'dprentice' gives the density, 'pprentice' gives the distribution 
     function, 'qprentice' gives the quantile function, and 'rprentice'
     generates random deviates.
}
\references{
  Prentice, R. L. (1975) \emph{Discrimination among some parametric models}.
  Biometrika, 62(3):607-614.
}
\examples{
  set.seed(1234)
  m <- c(2,0.5)
  x <- seq(-5,5,0.01)
  y <- pprentice(x,m)
  plot(x,y,type="l",xlab="x",ylab="pprentice")
}
