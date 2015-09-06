# --------------------------------------------------------
# Class MImputedExpressionSets
# holds ExpressionSets of multiple imputed data
# --------------------------------------------------------


# return the original Expression Set wherefrom the imputed values were created
setMethod("getOriginalData", "MImputedExpressionSets", definition=function(object) {
   return(object@originalData)
})

# minimum of present
setMethod(minPresent, "MImputedExpressionSets", function(object) slot(object, "minPresent"))
setReplaceMethod("minPresent", "MImputedExpressionSets", function(object, value) {
   slot(object, "minPresent") <- value
   validObject(object)
   return(object)
})


# setMethod(data, "MImputedExpressionSets", function(object) slot(object, "data"))
#
# setReplaceMethod("data", "MImputedExpressionSets", function(object, value) {
#    slot(object, "data") <- value
#    validObject(object)
#    return(object)
# })


setMethod(groupingCol, "MImputedExpressionSets", function(object) slot(object, "groupingCol"))
setReplaceMethod("groupingCol", "MImputedExpressionSets", function(object, value) {
   slot(object, "groupingCol") <- value
   validObject(object)
   return(object)
})


setMethod("getGroupingCol", "MImputedExpressionSets", function(object) {
   return(object@groupingCol)
})

setMethod(numberImputations, "MImputedExpressionSets", function(object) slot(object, "numberImputations"))
setReplaceMethod("numberImputations", "MImputedExpressionSets", function(object, value) {
   slot(object, "numberImputations") <- value
   validObject(object)
   return(object)
})


##' @importClassesFrom Biobase AnnotatedDataFrame
##' @export
setMethod("pData", "MImputedExpressionSets", function(object) {
   return(pData(object@data[[1]]))
})

##' @importClassesFrom Biobase AnnotatedDataFrame
##' @export
setMethod("fData", "MImputedExpressionSets", function(object) {
   return(fData(object@data[[1]]))
})

##' @importClassesFrom Biobase AnnotatedDataFrame
##' @export
setMethod("annotation", "MImputedExpressionSets", function(object) {
   return(annotation(object@data[[1]]))
})

##' @importClassesFrom Biobase AnnotatedDataFrame
##' @export
setMethod("eset", "MImputedExpressionSets", function(object, imputation) {
   return(object@data[[imputation]])
})

##' @importClassesFrom Biobase AnnotatedDataFrame
##' @export
setMethod("intensities", "MImputedExpressionSets", function(object, imputation) {
   return(exprs(object@data[[imputation]]))
})

##' @importClassesFrom Biobase AnnotatedDataFrame
##' @export
setMethod("limmasFit", "MImputedExpressionSets", function(object, design) {
   return(new("MImputedMArrayLM", data=lapply(object@data, function(x) {
      # This is needed so Limma ignores NAs for building a model
      exprs(x)[is.na(exprs(x))] <- NaN
      fit  <- lmFit(x, design)
   })))
})

##' CompleteCase filter
##'
##' @return returns all complete rows
##'
##' @importClassesFrom Biobase ExpressionSet
##' @export
setMethod("completeCases", "MImputedExpressionSets", function(object) {
   return(filterRows(object, apply(exprs(object@data[[1]]), 1, function(x) !any(is.na(x)))))
})

##' filterRows filter
##'
##' @return filter all rows
##'
##' @importClassesFrom Biobase ExpressionSet
##' @export
setMethod("filterRows", signature = c(object = "MImputedExpressionSets", filter = "logical"), function(object, filter) {
   object@data  <- lapply(object@data, function(x) filterRows(x, filter))
   return(object)
})

##' filterCols filter
##'
##' @return filter all columns
##'
##' @importClassesFrom Biobase ExpressionSet
##' @export
setMethod("filterCols",  c(object = "MImputedExpressionSets", filter = "logical"), function(object, filter) {
   object@data  <- lapply(object@data, function(x) filterCols(x, filter))
   return(object)
})

##' fill missing values with values
##'
##' @param value value to replace NAs
##' @return MIputedExpressionSets with replaced NAs
##' @export
setMethod("fillNAsWithValues", c(object = "MImputedExpressionSets", value = "numeric"), function(object, value) {
   object@data  <- lapply(object@data, function(x) fillNAs(x, value))
   return(object)
})

##' plotExpression
##'
##' function to plot a feature
##' @param ID character identification of feature ID
##' @export
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

##' plotExpression
##'
##' function to plot a feature
##' @param ID character identification of feature ID
##' @param groupingCol grouping Col
##' @param colors colors
##' @param samplelabels labels for plot
##' @export
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
                                                   max(testPrint)+0.3), col="white", xaxt="n",
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

##' get the average expression
##' @importClassesFrom Biobase ExpressionSet
##' @return returns an ExpressionSet with average intensities
##' @export
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

##' selects subset of MImputedExpressionSets
##'
##' @return subset MImputedExpressionSets
##' @export
setMethod("[", "MImputedExpressionSets", function(x, i, j, ..., drop = TRUE) {
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

##' calc group estimations
##' ß
setMethod("calcGroupEstimations", "MImputedExpressionSets", function(object) {
    # for each group, calculate the estimates
    object@groupData <- lapply(levels(pData(object@originalData)[,object@groupingCol]),
                               function(x) {
                                     getGroupData(object@originalData,
                                                  group = x,
                                                  groupCol = object@groupingCol)
                               })

})

##' @export
setMethod("estimateLimits", "MImputedExpressionSets", function(object, design, contrasts) {
   #if(is.null(object@groupData)) {
   calcGroupEstimations(object)
   #}

   contrastMatrix <- makeContrasts(contrasts=contrasts, levels=design)

   contrastMatrix <- makeContrasts(contrasts=contrasts, levels=design)
   apply(contrastMatrix, 2, function(x) {
      if(length(unique(x)) > 2 || length(x[x=="-1"] != 1) || length(x[x=="1"] != 1)) {
         stop("Only two conditions can be compared directly at this point.
               Comparison of multiple conditions will be provided in future
               versions. If you want to contribute to this package, please
               contact the authors. Thank you.")
      }
   })
})




