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

.lik.aranda <- function(alpha, beta, y, X) {
  if(alpha > 0){
    eta <- X %*% beta
    aranda <- ifelse(alpha*exp(eta) > -1, 1 - (alpha*exp(eta) + 1)^(-1/alpha), 1)
    aranda <- pmin(pmax(aranda, .Machine$double.eps), 1 - .Machine$double.eps)
    return(-sum(dbern(y, aranda, log = TRUE)))
  } else{
    return(1e100)
  }
}


.lik.weibull = function(gama, beta, y, X){
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

.lik.prentice <- function(m = c(1, 1), beta, y, X) {
  if(all(m > 0)){
    eta <- X %*% beta
    prentice <- pprentice(eta, m)
    prentice <- pmin(pmax(prentice, .Machine$double.eps), 1 - .Machine$double.eps)
    return(-sum(dbern(y, prentice, log = TRUE)))
  } else {
    return(1e100)
  }
}

.lik.stukel <- function(alpha, beta, y, X) {
  eta <- X %*% beta
  stukel <- pstukel(eta, alpha)
  stukel <- pmin(pmax(stukel, .Machine$double.eps), 1 - .Machine$double.eps)
  return(-sum(dbern(y, stukel, log = TRUE)))
}

