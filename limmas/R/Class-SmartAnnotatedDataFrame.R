setClass("SmartAnnotatedDataFrame",
         representation(sampleNamesCol="character",
                        originalNamesCol="character"
         ),
         contains="AnnotatedDataFrame",
)


setGeneric("originalNamesCol", function(object) standardGeneric("originalNamesCol"))
#setGeneric("originalNamesCol<-", function(object, value) standardGeneric("originalNamesCol<-"))
setMethod("originalNamesCol", "SmartAnnotatedDataFrame", function(object) {
   slot(object, "originalNamesCol")
})
#setReplaceMethod("originalNamesCol<-", "SmartAnnotatedDataFrame", function(object, value){
#   slot(object, "originalNamesCol") <- value
#   validObject(object)
#   return(object)
#})


setGeneric("sampleNamesCol", function(object) standardGeneric("sampleNamesCol"))
setGeneric("sampleNamesCol<-", function(object, value) standardGeneric("sampleNamesCol<-"))
setMethod(sampleNamesCol, "SmartAnnotatedDataFrame", function(object) slot(object, "sampleNamesCol"))
setReplaceMethod("sampleNamesCol", "SmartAnnotatedDataFrame", function(object, value){
   slot(object, "sampleNamesCol") <- value
   validObject(object)
   return(object)
})


setGeneric("getOriginalNames", function(object) standardGeneric("getOriginalNames"))
setMethod("getOriginalNames", "SmartAnnotatedDataFrame", function(object) {
   if(object@originalNamesCol == "") {
      return(rownames(pData(object)))
   } else  {
      return(as.character(pData(object)[,object@originalNamesCol]))
   }
})

setGeneric("getSampleNames", function(object) standardGeneric("getSampleNames"))
setMethod("getSampleNames", "SmartAnnotatedDataFrame", function(object) {
   if(object@sampleNamesCol == "") {
      return(rownames(pData(object)))
   } else  {
      return(as.character(pData(object)[,object@sampleNamesCol]))
   }
})

# assertions
setValidity("SmartAnnotatedDataFrame", function(object) {
   msg <- NULL

   if(!is.character(object@sampleNamesCol)) {
      msg <- c(msg, "sampleNamesCol is not a character string") 
   }
   if(!is.character(object@originalNamesCol)) {
      msg <- c(msg, "originalNamesCol is not a character string") 
   }
   if(object@sampleNamesCol != "" && !(object@sampleNamesCol %in% colnames(pData(object)))) {
      msg <- c(msg, "sampleNamesCol is not a column name.") 
   } else {
      if(any(make.names(getSampleNames(object)) != getSampleNames(object))) {
         msg <- c(msg, "sampleNames have to be syntactically correct. look up make.names for specification of syntactically valid names.")    
      }
   }
   if(object@originalNamesCol != "" && !(object@originalNamesCol %in% colnames(pData(object)))) {
      msg <- c(msg, "originalNamesCol is not a column name.") 
   }
   if(!nrow(pData(object)) > 1) {
      msg <- c(msg, "pheno has to have more than 1 row.") 
   }
   
   if (is.null(msg)) TRUE else msg
})            



setGeneric("getAnnotatedDataFrame", function(object) standardGeneric("getAnnotatedDataFrame"))
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




# Overwriting the [ operator so that factors are releveled
setMethod("[", "SmartAnnotatedDataFrame", function(x,i,j, ..., drop = TRUE){
   if (missing(drop)) 
      drop = FALSE
   else if (drop) 
      stop("'AnnotatedDataFrame' does not support drop = TRUE")
   if (missing(j)) {
      mD <- x@varMetadata
      pD <- x@data[i, , drop = drop]
   }
   else {
      mD <- x@varMetadata[j, , drop = drop]
      if (missing(i)) 
         pD <- x@data[, j, drop = drop]
      else pD <- x@data[i, j, drop = drop]
   }
   
   if(!missing(i)) {
      # relevel the factors 
      pD <- as.data.frame(apply(pD, 2, function(y) { 
         if(is.factor(y)) {
            y <- relevel(y)
         }
         
         return(y)
      }))
      
      # assure that the names are still character strings
      if(sampleNamesCol(x) != "") {
         x@data[,sampleNamesCol(x)] <-  as.character(x@data[,sampleNamesCol(x)])
      }
      
      if(originalNamesCol(x) != "") {
         x@data[,originalNamesCol(x)] <-  as.character(x@data[,originalNamesCol(x)])
      }
   }
   
   initialize(x, data = pD, varMetadata = mD)  
})


# Friendly constructor
SmartAnnotatedDataFrame <- function(sampleNamesCol = "",originalNamesCol = "", ...){
   return(new(Class="SmartAnnotatedDataFrame",sampleNamesCol=sampleNamesCol,originalNamesCol=originalNamesCol, ...))
}

