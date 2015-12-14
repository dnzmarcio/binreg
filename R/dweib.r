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

# Weibull density function
dweib <- function(x, gamma = 1) {
  if(length(gamma) != 1)
    stop("gamma = ", deparse(substitute(gamma)), " should have exactly one element.")

  x <- pmax(x, 0) 
  return(gamma*x^(gamma - 1)*exp(-x^gamma))
}

# Weibull distribution function
pweib <- function(q, gamma = 1) {
  if(length(gamma) != 1)
    stop("gamma = ", deparse(substitute(gamma)), " should have exactly one element.")

  q <- pmax(q, 0)
  return(1 - exp(-(q^gamma)))
}

# Weibull quantile function
qweib <- function(p, gamma = 1) {
  if(length(gamma) != 1)
    stop("gamma = ", deparse(substitute(gamma)), " should have exactly one element.")

  return((-log(1 - p))^(1/gamma))
}

# random generate from Weibull distribution
rweib <- function(n, gamma = 1) {
  if(length(gamma) != 1)
    stop("gamma = ", deparse(substitute(gamma)), " should have exactly one element.")

  return(qweib(runif(n),gamma))
}

