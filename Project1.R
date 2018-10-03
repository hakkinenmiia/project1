##################################################
#### TG-GATE Open database #######################

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("affy")
browseVignettes("affy")
### R code from vignette source 'vim.Rnw'

###################################################
### code chunk number 1: vim.Rnw:42-43
###################################################

library(affy)

###################################################
### code chunk number 2: <
###################################################
getNrowForCEL <- function(x) max(getPosXForCEL(x))
getNcolForCEL <- function() max(getPosYForCEL())

import.celfile <- function(celfile, ...) {
  cel.nrow <- getNrowForCEL(celfile)
  cel.ncol <- getNcolForCEL(celfile)
  x <- matrix(NA, nr=cel.nrow, nc=cel.ncol)
  cel.intensities <- getIntensitiesForCEL(celfile)
  cel.posx <- getPosXForCEL(celfile) # +1 if indexing starts at 0 (like in .CEL)
  cel.posy <- getPosYForCEL(celfile) # idem
  x[cbind(cel.posx, cel.posy)] <- cel.intensities
  mycdfName <- whatcdf("aCELfile.CEL")
  myCel <- new("Cel", exprs=x, cdfName=mycdfName)
  return(myCel)
}
import.celfile("003017549001.CEL")



#AffyBatch object
#size of arrays=834x834 features (18 kb)
#cdf=Rat230_2 (31099 affyids)
#number of samples=1
#number of genes=31099
#annotation=rat2302
