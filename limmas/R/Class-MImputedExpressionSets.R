# --------------------------------------------------------
# Class MImputedExpressionSets
# holds ExpressionSets of multiple imputed data
# --------------------------------------------------------
setClass("MImputedExpressionSets",
            representation=representation(
            data="list",
            groupingCol="character",
            minPresent="numeric",
            numberImputations="numeric",
            originalData="ExpressionSet")
         )
            

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

setGeneric("getOriginalData", function(object) standardGeneric("getOriginalData"))
setMethod(f="getOriginalData", signature="MImputedExpressionSets", definition=function(object) {
   return(object@originalData)
})


setGeneric("minPresent", function(object) standardGeneric("minPresent"))
setGeneric("minPresent<-", function(object, value) standardGeneric("minPresent<-"))
setMethod(minPresent, "MImputedExpressionSets", function(object) slot(object, "minPresent"))
setReplaceMethod("minPresent", "MImputedExpressionSets", function(object, value) {
   slot(object, "minPresent") <- value
   validObject(object)
   return(object)
})

setGeneric("data", function(object) standardGeneric("data"))
setGeneric("data<-", function(object, value) standardGeneric("data<-"))
setMethod(data, "MImputedExpressionSets", function(object) slot(object, "data"))
setReplaceMethod("data", "MImputedExpressionSets", function(object, value) {
   slot(object, "data") <- value
   validObject(object)
   return(object)
})

setGeneric("groupingCol", function(object) standardGeneric("groupingCol"))
setGeneric("groupingCol<-", function(object, value) standardGeneric("groupingCol<-"))
setMethod(groupingCol, "MImputedExpressionSets", function(object) slot(object, "groupingCol"))
setReplaceMethod("groupingCol", "MImputedExpressionSets", function(object, value) {
   slot(object, "groupingCol") <- value
   validObject(object)
   return(object)
})
setGeneric("getGroupingCol", function(object) standardGeneric("getGroupingCol"))
setMethod("getGroupingCol", "MImputedExpressionSets", function(object) {
   return(object@groupingCol)
})

setGeneric("numberImputations", function(object) standardGeneric("numberImputations"))
setGeneric("numberImputations<-", function(object, value) standardGeneric("numberImputations<-"))
setMethod(numberImputations, "MImputedExpressionSets", function(object) slot(object, "numberImputations"))
setReplaceMethod("numberImputations", "MImputedExpressionSets", function(object, value) {
   slot(object, "numberImputations") <- value
   validObject(object)
   return(object)
})

#setGeneric("pData", function(object) standardGeneric("pData"))
setMethod("pData", "MImputedExpressionSets", function(object) {
   return(pData(object@data[[1]]))
})

#setGeneric("fData", function(object) standardGeneric("fData"))
setMethod("fData", "MImputedExpressionSets", function(object) {
   return(fData(object@data[[1]]))
})

#setGeneric("annotation", function(object) standardGeneric("annotation"))
setMethod("annotation", "MImputedExpressionSets", function(object) {
   return(annotation(object@data[[1]]))
})

setGeneric("eset", function(object, imputation) standardGeneric("eset"))
setMethod("eset", "MImputedExpressionSets", function(object, imputation) {
   return(object@data[[imputation]])
})

setGeneric("intensities", function(object, imputation) standardGeneric("intensities"))
setMethod("intensities", "MImputedExpressionSets", function(object, imputation) {
   return(exprs(object@data[[imputation]]))
})

setGeneric("limmasFit", function(object, design) standardGeneric("limmasFit"))
setMethod("limmasFit", "MImputedExpressionSets", function(object, design) {
   return(new("MImputedMArrayLM", data=lapply(object@data, function(x) {
      fit  <- lmFit(x, design)
   })))
})


# completeCase filter
setGeneric("completeCases", function(object) standardGeneric("completeCases"))
setMethod("completeCases", "MImputedExpressionSets", function(object) {
   return(filterRows(object, apply(exprs(object@data[[1]]), 1, function(x) !any(is.na(x)))))
})

#setGeneric("filterRows", function(object, filter) standardGeneric("filterRows"))
setMethod("filterRows", "MImputedExpressionSets", function(object, filter) {
   object@data  <- lapply(object@data, function(x) filterRows(x, filter))
   return(object)
})

#setGeneric("filterCols", function(object, filter) standardGeneric("filterCols"))
setMethod("filterCols", "MImputedExpressionSets", function(object, filter) {
   object@data  <- lapply(object@data, function(x) filterCols(x, filter))
   return(object)
})

setGeneric("fillNAs", function(object, value) standardGeneric("fillNAs"))
setMethod("fillNAs", c(object="MImputedExpressionSets", value="numeric"), function(object, value) { 
   object@data  <- lapply(object@data, function(x) fillNAs(x, value))
   return(object)
})



setGeneric("plotExpression", function(object, ID, ...) standardGeneric("plotExpression"))
setMethod("plotExpression", c(object = "MImputedExpressionSets", ID = "character"), function(object, ID, ...) {  
   id <- which(annotation(object) == ID)
   if (length(id) == 0) {
      stop(paste("id '", ID, "' not found", sep=""))
   } else if (length(id) > 1) {
      stop(paste(ID, "has more than 1 entries:", paste(id, sep=","), sep=" "))
   } else {
      plotExpression(object, id, ...) 
   }
})


setMethod("plotExpression", c(object = "MImputedExpressionSets", ID = "numeric"), function(object, ID,
                                                   groupingCol = getGroupingCol(object),
                                                   colors = topo.colors(length(unique(pData(object)[,groupingCol]))),
                                                   samplelabels = rownames(pData(object)),
                                                   ... ) {  
   
   if (!is.character(groupingCol)) {
      stop("groupingCol is not a character string")   
   } 
   if(!groupingCol %in% colnames(pData(object))) {
      stop("groupingCol is not a column name in pData")  
   }
   if(length(ID) > 1) {
      stop(paste("ID has length ", length(ID) ,". Only one ID can be plotted at the same time.", sep=""))
   }
   if(is.null(colors) || is.na(colors)) { 
      colors = topo.colors(length(unique(pData(object)[,groupingCol])))
   }
      
   groups <- as.character(pData(object)[,groupingCol])   
   printTable <- do.call(rbind, lapply(object@data, function(x) exprs(x)[ID,]))

   if(length(samplelabels) != ncol(printTable)) {
      stop("samplelabels is not equal to the length of samples")
   }
   if(length(colors) != length(unique(groups))) {
      stop("length of colors is not equal to the length of groups")  
   }
   
   testPrint <- printTable[1,]
   testPrint[is.na(testPrint)] <- 0
   plot(1:ncol(printTable), testPrint, ylim=c(min(testPrint),
                                                   max(testPrint)), col="white", xaxt="n", 
        main = paste("Expression of ", annotation(object)[ID], sep=""),
        xlab = "", ylab = "log intensities", ...)
   
   axis(1, at=1:ncol(printTable), labels=samplelabels, las=2)
   
   for(groupid in 1:length(unique(groups))) {
      position <- which(groups == unique(groups)[groupid])
      color <- colors[groupid]
      type <- "b"
      pch  <- 1
      if (all(is.na(printTable[i, position]))) {
         printTable[i, position] <- 0
         type <- "p"
         pch  <- "#"
      }
      
      for(i in 1:nrow(printTable)) {
         lines(position, printTable[i, position], type=type, pch=pch, col=color)
      }
   }
})


setGeneric("getAverageExpression", function(object) standardGeneric("getAverageExpression"))
setMethod("getAverageExpression", "MImputedExpressionSets", function(object) {
   mat <- intensities(object, 1)
   for(i in 2:numberImputations(object)) { 
      mat <- mat + intensities(object, i)
   }
   
   mat <- mat / numberImputations(object)
   
   eSet <- eset(object, 1)
   exprs(eSet) <- mat
   return(eSet)
})


setMethod("[", "MImputedExpressionSets", function(x,i,j, ..., drop = TRUE) {
   if (missing(drop)) 
      drop = FALSE
   else if (drop) 
      stop("'MImputedExpressionSets' does not support drop = TRUE")
   
   if (!missing(i)) {
      x <- filterRows(x, i)
   }
   if (!missing(j)) {
      x <- filterCols(x, j)  
   }
   
   return(x)
})

