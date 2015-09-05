##' Create an expression set expression set and annotated data frame
##'
##' @title createExpressionSet
##' @param intensities intensity matrix
##' @param pheno pheno information
##' @param features feature meta data
##' @param annotation feature annotations
##' @importClassesFrom Biobase ExpressionSet AnnotatedDataFrame
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

