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

# Prentice density function
dprentice <- function(x, m = c(1, 1)) {
  if(length(m) < 2)
    stop("m = ", deparse(substitute(m)), " should have exactly two elements.")
  if(length(m) > 2)
    warning("Length of m = ", deparse(substitute(m)), "is greater than 2. Only the first two elements will be used.")
  
  return(exp(x*m[1])*((1 + exp(x))^(-(m[1] + m[2])))/beta(m[1], m[2]))
}

# Prentice distribution function
pprentice <- function(q, m = c(1, 1)) {
  if(length(m) < 2)
    stop("m = ", deparse(substitute(m)), " should have exactly two elements.")
  if(length(m) > 2)
    warning("Length of m = ", deparse(substitute(m)), "is greater than 2. Only the first two elements will be used.")
  
  return( pbeta(plogis(q), m[1], m[2]) )
}

# Prentice quantile function
qprentice <- function(p, m = c(1, 1)) {
  if(length(m) < 2)
    stop("m = ", deparse(substitute(m)), " should have exactly two elements.")
  if(length(m) > 2)
    warning("Length of m = ", deparse(substitute(m)), "is greater than 2. Only the first two elements will be used.")
  
  return( qlogis(qbeta(p, m[1], m[2])) )
}

# random generate from Prentice distribution
rprentice <- function(n, m = c(1, 1)) {
  if(length(m) < 2)
    stop("m = ", deparse(substitute(m)), " should have exactly two elements.")
  if(length(m) > 2)
    warning("Length of m = ", deparse(substitute(m)), "is greater than 2. Only the first two elements will be used.")

  return(qprentice(runif(n),m))
}



