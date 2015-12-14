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

# Stukel density function
dstukel <- function(x,alpha = c(0,0)) {
  if(length(alpha) < 2)
    stop("alpha = ", deparse(substitute(alpha)), " should have exactly two elements.")
  if(length(alpha) > 2)
    warning("Length of alpha = ", deparse(substitute(alpha)), "is greater than 2. Only the first two elements will be used.")

  return( dlogis(.hstukel(x, alpha))*.dhstukel(x, alpha) )
}

# Stukel distribution function
pstukel <- function(q,alpha = c(0,0)) {
  if(length(alpha) < 2)
    stop("alpha = ", deparse(substitute(alpha)), " should have exactly two elements.")
  if(length(alpha) > 2)
    warning("Length of alpha = ", deparse(substitute(alpha)), "is greater than 2. Only the first two elements will be used.")

  return( plogis(.hstukel(q, alpha)) )
}

# Stukel quantile function
qstukel <- function(p,alpha = c(0,0)) {
  if(length(alpha) < 2)
    stop("alpha = ", deparse(substitute(alpha)), " should have exactly two elements.")
  if(length(alpha) > 2)
    warning("Length of alpha = ", deparse(substitute(alpha)), "is greater than 2. Only the first two elements will be used.")
  
  return( .invhstukel(qlogis(p), alpha) )
}

# random generate from Stukel distribution
rstukel <- function(n,alpha = c(0,0)) {
  if(length(alpha) < 2)
    stop("alpha = ", deparse(substitute(alpha)), " should have exactly two elements.")
  if(length(alpha) > 2)
    warning("Length of alpha = ", deparse(substitute(alpha)), "is greater than 2. Only the first two elements will be used.")
  
  return(qstukel(runif(n),alpha))
}


############################################
.hstukel <- function(eta, alpha) {
  if(length(alpha) < 2)
    stop("alpha = ", deparse(substitute(alpha)), " should have exactly two elements.")
  if(length(alpha) > 2)
    warning("Length of alpha = ", deparse(substitute(alpha)), "is greater than 2. Only the first two elements will be used.")
  
  h <- function(eta){
    if(eta > 0){
      if(alpha[1] > 0) out <- (exp(alpha[1]*eta) - 1)/alpha[1]
      else if(alpha[1] == 0) out <- eta
      else out <- -log(1 - alpha[1]*eta)/alpha[1]
    } else {
      if(alpha[2] > 0) out <- -(exp(-alpha[2]*eta) - 1)/alpha[2]
      else if(alpha[2] == 0) out <- eta
      else out <- log(1 + alpha[2]*eta)/alpha[2]
    }
    return(out)
  }
  
  return(sapply(eta, h))
}

.dhstukel <- function(eta, alpha) {
  if(length(alpha) < 2)
    stop("alpha = ", deparse(substitute(alpha)), " should have exactly two elements.")
  if(length(alpha) > 2) 
    warning("Length of alpha = ", deparse(substitute(alpha)), "is greater than 2. Only the first two elements will be used")
  
  h <- function(eta){
    if(eta > 0){
      if(alpha[1] > 0) out <- exp(alpha[1]*eta)
      else if(alpha[1] == 0) out <- 1
      else out <- 1/(1 - alpha[1]*eta)
    } else {
      if(alpha[2] > 0) out <- exp(-alpha[2]*eta)
      else if(alpha[2] == 0) out <- 1
      else out <- 1/(1 + alpha[2]*eta)
    }
    return(out)
  }
  
  return(sapply(eta, h))
}

.invhstukel <- function(y, alpha) {
  if(length(alpha) < 2)
    stop("alpha = ", deparse(substitute(alpha)), " should have exactly two elements.")
  if(length(alpha) > 2) 
    warning("Length of alpha = ", deparse(substitute(alpha)), "is greater than 2. Only the first two elements will be used")
  
  h <- function(y){
    if(y > 0){
      if(alpha[1] > 0) out <- log(y*alpha[1] + 1)/alpha[1]
      else if(alpha[1] == 0) out <- y
      else out <- -(exp(-y*alpha[1]) - 1)/alpha[1]
    } else {
      if(alpha[2] > 0) out <- -log(1 - y*alpha[2])/alpha[2]
      else if(alpha[2] == 0) out <- y
      else out <- (exp(y*alpha[2]) - 1)/alpha[2]
    }
  }
  
  return(sapply(y, h)) 
}

