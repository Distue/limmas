##' Reads in a pheno file text file.
##'
##' @title read.pheno
##' @param file fileame
##' @param sampleNamesCol column describing the sample names
##' @param originalNamesCol column describing the original names
##' @return SmartAnnotatedDataFrame pheno information
##' @export
read.pheno <- function(file, sampleNamesCol="", originalNamesCol="intensityCols", header=T, row.names=1, sep="\t", ...) {
   pheno <- read.delim(file, header=header, row.names=row.names, sep=sep, ...)
   rownames(pheno) <- gsub(" ", ".", rownames(pheno))
   return(SmartAnnotatedDataFrame(originalNamesCol=originalNamesCol, sampleNamesCol=sampleNamesCol, pheno))
}

