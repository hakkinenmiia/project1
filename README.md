##### project1
### Genomics and Bioinformatics
### TG-GATE Open database ######

#### try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("affy")
browseVignettes("affy")
### R code from vignette source 'vim.Rnw'
### code chunk number 1: vim.Rnw:42-43
###################################################

library(affy)

####TO IMPORT .CEL FILE DATA ####
name.the.data <- read.affybatch("")
name.of.matrix = exprs(name.the.data)  #can also use  intensity(name.the.data)
#Both exprs and intensity do the same thing (note you got identical results!). 
#They extract the matrix of expression values, which is in some sense 'the contents of a .cel file'.
