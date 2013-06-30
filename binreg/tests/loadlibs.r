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

## This code is used for testing purposes. The binreg library does not
## depend on it for any of its functionalities

installpacks <- function(loc=NULL,repos="http://stat.ethz.ch/CRAN/") {
  ## set the repository to use
  options(repos=repos)
  ## install the packages
  install.packages("LaplacesDemon",lib=loc)
#  install.packages("coda",lib=loc)
#  install.packages("mcmc",lib=loc)
#  install.packages("R.methodsS3",lib=loc)
#  install.packages("R.oo",lib=loc)
#  install.packages("R.utils",lib=loc)
  ## this following line install the binreg package itself, so nothing else is needed.
  ## For testing, sometimes it is better to work without installing it for a while...
  ##      install.packages('./binreg_version.tar.gz',repos=NULL,type="source")
}

loadlibs <- function(libdir=NULL) {
  w <- options("warn")
  options("warn" = -1)
  if (require("binreg",quietly=TRUE)==FALSE) {
    library("LaplacesDemon",lib.loc=libdir)
#    library("coda",lib.loc=libdir)
#    library("mcmc",lib.loc=libdir)
#    library("R.methodsS3",lib.loc=libdir)
#    library("R.oo",lib.loc=libdir)
#    library("R.utils",lib.loc=libdir)
    source("../R/binreg.r")
    source("../R/dprentice.r")
    source("../R/dstukel.r")
    source("../R/daranda.r")
    source("../R/dweib.r")
    source("../R/lik.r")
    source("../R/link.r")
    source("../R/bayes.r")
  } else {
    library("binreg")
  }
  options("warn" = w[[1]])
}
