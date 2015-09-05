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
