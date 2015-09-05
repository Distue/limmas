##' get the column name containing the original sample names
##'
##' @importClassesFrom Biobase AnnotatedDataFrame
##' @return column name
##' @export
setMethod("originalNamesCol", "SmartAnnotatedDataFrame", function(object) {
   return(slot(object, "originalNamesCol"))
})

##' set the column name containing the original sample names
##'
##' @importClassesFrom Biobase AnnotatedDataFrame
##' @param value name
##' @export
setReplaceMethod("originalNamesCol", "SmartAnnotatedDataFrame", function(object, value){
   slot(object, "originalNamesCol") <- value
   validObject(object)
   return(object)
})

##' get the column name containing the sample names
##'
##' @importClassesFrom Biobase AnnotatedDataFrame
##' @return column name
##' @export
setMethod("sampleNamesCol", "SmartAnnotatedDataFrame", function(object) {
   return(slot(object, "sampleNamesCol"))
})

##' set the column name containing the sample names
##'
##' @importClassesFrom Biobase AnnotatedDataFrame
##' @param value name
##' @export
setReplaceMethod("sampleNamesCol", "SmartAnnotatedDataFrame", function(object, value){
   slot(object, "sampleNamesCol") <- value
   validObject(object)
   return(object)
})

##' get original sample names
##'
##' @importClassesFrom Biobase AnnotatedDataFrame
##' @return vector of original sample names
##' @export
setMethod("getOriginalNames", "SmartAnnotatedDataFrame", function(object) {
   if(object@originalNamesCol == "") {
      return(rownames(pData(object)))
   } else  {
      return(as.character(pData(object)[,object@originalNamesCol]))
   }
})

##' get new sample names
##'
##' @importClassesFrom Biobase AnnotatedDataFrame
##' @return vector of newsample names
##' @export
setMethod("getSampleNames", "SmartAnnotatedDataFrame", function(object) {
   if(object@sampleNamesCol == "") {
      return(rownames(pData(object)))
   } else  {
      return(as.character(pData(object)[,object@sampleNamesCol]))
   }
})


##' get AnnotatedDataFrame
##'
##' @importClassesFrom Biobase AnnotatedDataFrame
##' @return AnnotatedDataFrame
##' @export
setMethod("getAnnotatedDataFrame", "SmartAnnotatedDataFrame", function(object) {
   data <- pData(object)
   # do nothing if the original and the sample names are the same as rownames of pheno
   if(!(object@originalNamesCol == "" && object@sampleNamesCol == "")) {

      # if the sample names are rownames and the original names are in a column, remove the column,
      # unless 1 column would be left
      # (data.frames of column size of 1 will be converted to bloody vectors)
      if(object@sampleNamesCol == "" && !object@originalNamesCol == "") {
         if(ncol(data) > 2) {
            data <- data[,!colnames(data) %in% object@originalNamesCol]
         }

      # if the sample names are in a column and the original names are as row names, overwrite row names and remove
      } else if(!object@sampleNamesCol == "" && object@originalNamesCol == "") {
         rownames(data) <- data[,object@sampleNamesCol]
         data <- data[,!colnames(data) %in% object@sampleNamesCol ]

      # if both are in columns
      } else {
         rownames(data) <- data[,object@sampleNamesCol]
         data <- data[,!colnames(data) %in% c(object@originalNamesCol, object@sampleNamesCol)]
      }
   }

   return(new("AnnotatedDataFrame", data))
})



##' Overwriting the [ operator so that factors are releveled
##'
##' @importClassesFrom Biobase AnnotatedDataFrame
##' @export
setMethod("[", "SmartAnnotatedDataFrame", function(x, i, j, ..., drop = TRUE){
   if (missing(drop)) {
      drop = FALSE
   } else if (drop) {
      stop("'AnnotatedDataFrame' does not support drop = TRUE")
   }

   # missing j
   if (missing(j)) {
      mD <- x@varMetadata
      pD <- x@data[i, , drop = drop]
   } else {
      mD <- x@varMetadata[j, , drop = drop]

      # missing i, j available
      if (missing(i)) {
         pD <- x@data[, j, drop = drop]
      # missing i, missing j
      } else {
         pD <- x@data[i, j, drop = drop]
      }
   }

   # columns selecter
   if(!missing(i)) {
      # relevel the factors
      pD <- as.data.frame(apply(pD, 2, function(y) {
         if(is.factor(y)) {
            y <- relevel(y)
         }

         return(y)
      }))

      # assure that the names are still character strings
      if(sampleNamesCol(x) != "" && !is.null(sampleNamesCol(x))) {
         x@data[,sampleNamesCol(x)] <-  as.character(x@data[,sampleNamesCol(x)])
      }

      if(originalNamesCol(x) != "" && !is.null(originalNamesCol(x))) {
         x@data[,originalNamesCol(x)] <- as.character(x@data[,originalNamesCol(x)])
      }
   }

   initialize(x, data = pD, varMetadata = mD)
})


##' Constructor
##'
##' @importClassesFrom Biobase AnnotatedDataFrame
##' @param sampleNamesCol column name containing sample names, "" for rownames
##' @param originalNamesCol column name containing original names, "" for rownames
##' @export
SmartAnnotatedDataFrame <- function(sampleNamesCol = "",originalNamesCol = "", ...){
   return(new(Class="SmartAnnotatedDataFrame",
              sampleNamesCol=sampleNamesCol,
              originalNamesCol=originalNamesCol,
              ...))
}

