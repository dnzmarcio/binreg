% File binreg/man/b.laplace.Rd
% Part of the binreg package for R (http://www.R-project.org)
% Copyright (C) 2013 Bernardo dos Santos, Adriano Polpo, 
%                    Carlos A. de B. Pereira.
% Distributed under GPL 3 or later

\name{ModelLogit}

%%%
\alias{ModelLogit}
\alias{ModelProbit}
\alias{ModelCloglog}
\alias{ModelAranda}
\alias{ModelWeibull}
\alias{ModelCWeibull}
\alias{ModelPrentice}
\alias{ModelStukel}

\title{Binary models for Bayesian analysis.}
\description{
  Binary models to be used with LaplacesDemon function to generate a sample from
  the posterior distribution of the model's parameters.
}
\usage{
  It is not to be directly used, but only with LaplacesDemon function.
}
\details{
  It is a function with the binary regression models to facilitate de use of 
  LaplacesDemon function, and to perform Bayesian analysis. 
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
  data <- dataLD(y ~ x, "cloglog")
  fit <- LaplacesDemon(ModelCloglog(), data, GIV(ModelCloglog(), data),
                       Iterations=1000,Thinning=2)
}
