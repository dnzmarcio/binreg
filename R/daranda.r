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

# Aranda density function
daranda <- function(x, alpha = 1) {
  if(length(alpha) != 1)
    stop("alpha = ", deparse(substitute(alpha)), " should have exactly one element.")

  return( ifelse(alpha*exp(x) > -1, exp(x)*(alpha*exp(x) + 1)^(-1/alpha - 1), 0) )
}

# Aranda distribution function
paranda <- function(q, alpha = 1) {
  if(length(alpha) != 1)
    stop("alpha = ", deparse(substitute(alpha)), " should have exactly one element.")

  return( ifelse(alpha*exp(q) > -1, 1 - (alpha*exp(q) + 1)^(-1/alpha), 1) )
}

# Aranda quantile function
qaranda <- function(p, alpha = 1) {
  if(length(alpha) != 1)
    stop("alpha = ", deparse(substitute(alpha)), " should have exactly one element.")

  return( log(((1-p)^(-alpha)-1)/alpha) )
}

# random generate from Aranda distribution
raranda <- function(n, alpha = 1) {
  if(length(alpha) != 1)
    stop("alpha = ", deparse(substitute(alpha)), " should have exactly one element.")

  return(qaranda(runif(n),alpha))
}

