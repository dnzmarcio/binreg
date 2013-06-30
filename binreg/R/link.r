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

# Link: aranda
aranda <- function(alpha = 1) {
  linkfun <- function(mu)  qaranda(mu,alpha)
  linkinv <- function(eta) paranda(eta,alpha)
  mu.eta  <- function(eta) daranda(eta,alpha)
  link <- paste("Aranda-Ordaz, alpha =", alpha)
  valideta <- function(eta) TRUE
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta, name = link), 
            class = "link-glm")
}

# Link: weibull
weibull <- function(gamma){
  linkfun <- function(mu)  qweib(mu,gamma)
  linkinv <- function(eta) pweib(eta,gamma)
  mu.eta  <- function(eta) dweib(eta,gamma)
  link <- paste("Weibull, gamma =", gamma)
  valideta <- function(eta)
    return(all(eta > 0))
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta, name = link), 
            class = "link-glm")
}

# Link: prentice
prentice <- function(m = c(1, 1)) {
  linkfun <- function(mu)  qprentice(mu,m)
  linkinv <- function(eta) pprentice(eta,m)
  mu.eta  <- function(eta) dprentice(eta,m)
  link     <- paste("Prentice, m1 = ", m[1],"  -- m2 = ", m[2], sep = "")
  valideta <- function(eta) TRUE
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta, name = link), 
            class = "link-glm")
}

# Link: stukel
stukel <- function(alpha = c(0, 0)) {
  linkfun  <- function(mu)  qstukel(mu, alpha)
  linkinv  <- function(eta) pstukel(eta, alpha)
  mu.eta   <- function(eta) dstukel(eta, alpha)
  link     <- paste("Stukel, a1 = ", alpha[1],"  -- a2 = ",alpha[2], sep = "")
  valideta <- function(eta) TRUE
  structure(list(linkfun = linkfun, linkinv = linkinv, 
                 mu.eta = mu.eta, valideta = valideta, name = link), 
            class = "link-glm")
}

