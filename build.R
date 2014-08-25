detach("package:limmas")
remove.packages("limmas")
library(devtools)
setwd("~/bioc/limmas")
build("limmas", binary = T)
# install.packages("limmas", repos=NULL)
# R CMD install limmas_0.3.tgz
#promptPackage("limmas") 
