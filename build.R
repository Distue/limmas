detach("package:limmas")
remove.packages("limmas")
library(devtools)
setwd("~/bioc/limmas")
build("limmas", binary=T)
install.packages("limmas", repos=NULL)

#promptPackage("limmas") 
