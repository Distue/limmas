detach("package:limmas")
remove.packages("limmas")
library(devtools)
setwd("~/bioc/limmas")
build("limmas", binary = T)
# install.packages("limmas", repos=NULL)
# R CMD install limmas_0.3.tgz
#promptPackage("limmas") 

library(limmas)
load("limmas/data/gefh1inhib.rda")
data <- gefh1inhib
pheno <- gefh1inhib.pheno
pData(pheno)
originalNamesCol(pheno)
pheno.without.control <- pheno[,-c(1:3)]
