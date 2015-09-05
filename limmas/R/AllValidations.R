# assertions
setValidity("CombinedMArrayLM", function(object) {
   msg <- NULL   
   if (is.null(msg)) {
      TRUE
   } else {
      msg
   }
})     


# assertions
setValidity("MImputedExpressionSets", function(object) {
   msg <- NULL
   if(!is.list(object@data) || !all(unlist(lapply(object@data, function(x) is(x,"ExpressionSet"))))) {
      msg <- c(msg, "data is not a list of ExpressionSets")
   }
   if(!is.character(object@groupingCol)) {
      msg <- c(msg, "groupingCol is not a character string") 
   }
   if(!object@groupingCol %in% colnames(pData(object@data[[1]]))) {
      msg <- c(msg, "groupingCol is not in the columnnames of pheno") 
   }
   if(!is.numeric(object@minPresent) && object@minPresent > 0 && object@minPresent < 1) {
      msg <- c(msg, "minPresent has to be numeric and > 0 and < 1")
   }
   if(!is.numeric(object@numberImputations) && object@numberImputations > 0) {
      msg <- c(msg, "numberImputations has to be numeric and > 0")
   }
   if(object@numberImputations != length(object@data)) {
      msg <- c(msg, "numberImputations is not the length of data list")
   }
   if (is.null(msg)) TRUE else msg
})    


# assertions
setValidity("MImputedMArrayLM", function(object) {
   msg <- NULL
   
   if(!all(unlist(lapply(object@data, function(x) is(x, "MArrayLM"))))) {
      msg <- c(msg, "data list can only contain MArrayLM objects") 
   }
   if(!length(object@data) > 1) {
      msg <- c(msg, "data list has to have more than one object") 
   }
   
   if (is.null(msg)) TRUE else msg
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





