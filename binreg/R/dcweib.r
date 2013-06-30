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

# Reflected Weibull density function
dcweib <- function(x, gamma = 1) {
  if(length(gamma) != 1)
    stop("gamma = ", deparse(substitute(gamma)), " should have exactly one element.")

  x <- pmin(x, 0) 
  return(gamma*(-x)^(gamma - 1)*exp(-(-x)^gamma))
}

# Reflected Weibull distribution function
pcweib <- function(q, gamma = 1) {
  if(length(gamma) != 1)
    stop("gamma = ", deparse(substitute(gamma)), " should have exactly one element.")

  q <- pmin(q, 0)
  return(exp(-((-q)^gamma)))
}

# Reflected Weibull quantile function
qcweib <- function(p, gamma = 1) {
  if(length(gamma) != 1)
    stop("gamma = ", deparse(substitute(gamma)), " should have exactly one element.")

  return(-(-log(p))^(1/gamma))
}

# random generate from Reflected Weibull distribution
rcweib <- function(n, gamma = 1) {
  if(length(gamma) != 1)
    stop("gamma = ", deparse(substitute(gamma)), " should have exactly one element.")

  return(-qweib(runif(n),gamma))
}

