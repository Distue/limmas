# ------------------------------------------------------------------------------------
# Functions for reading in data and creating creating ExpressionSets
# ------------------------------------------------------------------------------------

# ------------------------------------------
# public functions
# ------------------------------------------
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

##' Read in a maxQuant proteinGroups file
##'  
##' @title read.maxQuant
##' @param file fileame
##' @param pheno pheno information 
##' @param idcolumn name of column containing the protein id
##' @param splitIds locial if the ids should be split 
##' @return ExpressionSet protein intensities and annotation as ExpressionSet
##' @export
read.maxQuant <- function(file, pheno, idcolumn="Protein.IDs", splitIds = F, ...) {
   if(!is.character(file)) {
      stop("file is not a character string")
   }
   if(!file.exists(file)) {
      stop(paste("file '", file, "' does not exist or is not accessable.", sep=""))
   }
   if(!is(pheno, "SmartAnnotatedDataFrame")) {
      stop("pheno an SmartAnnotatedDataFrame.")
   }
   if(!is.character(idcolumn)) {
      stop("idcolumn is not a character string.")
   }
   if(!is.logical(splitIds) && !is.character(splitIds)) {
      stop("splitIds has to be a character symbol")
   }
   
   data <- read.delim(file, header=T, stringsAsFactors=F, ...)
   exprset <- createExpressionSetFromMaxQuant(data, pheno, idcolumn=idcolumn)
   
   if(!is.logical(splitIds)) {
      annotation(exprset) <- unlist(lapply(strsplit(as.character(annotation(exprset)), splitIds), function(x) { if(length(x)>0) { x[[1]] } else { ""} }))
   }
   
   return(exprset)
}


##' Creates an ExpressionSet from the proteinGroups data frame
##'  
##' @title createExpressionSetFromMaxQuant
##' @param file fileame
##' @param pheno pheno information 
##' @param idcolumn name of column containing the protein id
##' @param annotation feature annotations
##' @return ExpressionSet protein intensities and annotation as ExpressionSet
##' @export
createExpressionSetFromMaxQuant <- function(data, pheno, idcolumn="Protein.IDs", annotation=as.character(data[,idcolumn]), ...) {
   if(!is.data.frame(data)) {
      stop("data is not a data frame.")
   }
   if(!is(pheno, "SmartAnnotatedDataFrame")) {
      stop("pheno is not a SmartAnnotatedDataFrame.")
   }
   
   rowsNotInData <- getOriginalNames(pheno)[!getOriginalNames(pheno) %in% colnames(data)]
   
   if(length(rowsNotInData) > 0) {
      stop(paste("pheno rownames not in data:", paste(rowsNotInData, sep=","), sep=" "))
   }
   
   # select intensities
   intensities <- data[,getOriginalNames(pheno)]
   colnames(intensities) <- getSampleNames(pheno)
   
   # select feature data
   features <- data[,!colnames(data) %in% getOriginalNames(pheno)]
   
   features <- new("AnnotatedDataFrame", data=features)
   
   return(createExpressionSet(as.matrix(intensities), pheno, features, annotation))
}

##' Create an expression set expression set and annotated data frame
##'  
##' @title createExpressionSet
##' @param intensities intensity matrix
##' @param pheno pheno information 
##' @param features feature meta data
##' @param annotation feature annotations
##' @return ExpressionSet protein intensities and annotation as ExpressionSet
##' @export
createExpressionSet <- function(intensities, pheno, features, annotation, ...) {
   if(!(is(pheno, "SmartAnnotatedDataFrame") || is(pheno, "AnnotatedDataFrame"))) {
      stop("pheno is not an object of class AnnotatedDataFrame or SmartAnnotatedDataFrame.")
   }
   if(!is(features, "AnnotatedDataFrame")) {
      stop("feature is not an object of class AnnotatedDataFrame")
   }
   if(!is.character(annotation)) {
      stop("annotation is not a vector of characters.")
   }
   if(is.data.frame(intensities)) {
      intensities <- as.matrix(intensities)
   }
   if(!is.matrix(intensities)) {
      stop("intensities is not a matrix or a data frame.")
   }
   if(is(pheno, "SmartAnnotatedDataFrame")) {
      pheno <- getAnnotatedDataFrame(pheno)
   }
   return(new("ExpressionSet", exprs=intensities, phenoData=pheno, featureData=features, annotation=annotation, ...))
}


# ------------------------------------------
# private functions
# ------------------------------------------

