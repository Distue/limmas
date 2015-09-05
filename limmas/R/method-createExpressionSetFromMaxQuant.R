
##' Creates an ExpressionSet from the proteinGroups data frame
##'
##' @title createExpressionSetFromMaxQuant
##' @param file fileame
##' @param pheno pheno information
##' @param idcolumn name of column containing the protein id
##' @param annotation feature annotations
##' @return ExpressionSet protein intensities and annotation as ExpressionSet
##' @importClassesFrom Biobase ExpressionSet AnnotatedDataFrame
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

