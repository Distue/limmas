# --------------------------------------------------------
# Class MImputedExpressionSets
# Authors: Thomas Schwarzl <schwarzl@embl.de> with help from Elisa D'Arcangelo
# holds ExpressionSets of multiple imputed data
# --------------------------------------------------------
# Constructor:
##' @rdname MImputedExpressionSets
##' @export
MImputedExpressionSets <- function(data,
                                   groupingCol,
                                   minPresent,
                                   numberImputations,
                                   orginalData,
                                   groupData,
                                   ...){
   return(new(Class            = "MImputedExpressionSets",
              data                 = data,
              groupingCol          = groupingCol,
              minPresent           = minPresent,
              numberImputations    = numberImputations,
              originalData         = originalData,
              groupData            = groupData,
              ...))
}
# --------------------------------------------------------

##' @name nrow
##' @title number of rows
##' @return numeric. number of rows
##' @export
setMethod("nrow", "MImputedExpressionSets", function(object) {
   nrow(eset(object, 1))
})

##' @name getOriginalData
##' @title get orginial data
##' @description return the original Expression Set wherefrom the imputed values were created
##' @export
setMethod("getOriginalData", "MImputedExpressionSets", function(object) {
   return(object@originalData)
})

##' @name minPresent
##' @title set and get minimum present
##' @param value minimal present
##' @description minimum of present
##' @return minPresent
##' @export
setMethod(minPresent, "MImputedExpressionSets", function(object) {
   slot(object, "minPresent")
})


##' @rdname minPresent
##' @export
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

##' @name groupingCol
##' @title get and set grouping columns
##' @description TODO
##'
##' Getter:
##'    groupingCol(x)
##'    getGroupingCol(x)
##' Setter:
##'    setGroupingCol(x, value)
##'    groupingCol(x) <- value
##' @param value value
##' @return grouping column as character
##' @export
setMethod(groupingCol, "MImputedExpressionSets", function(object) {
   slot(object, "groupingCol")
})

##' @rdname groupingCol
##' @export
setReplaceMethod("groupingCol", "MImputedExpressionSets", function(object, value) {
   slot(object, "groupingCol") <- value
   validObject(object)
   return(object)
})

##' @rdname groupingCol
##' @export
setMethod("getGroupingCol", "MImputedExpressionSets", function(object) {
   return(object@groupingCol)
})

##' @name numberImputations
##' @title get or set number of imputations
##' @description
##'
##' Getter
##'   numberImputations(x)
##' Setter
##'   numberImputations(x) <- value
##' @return number of imputations as ßnumeric
##' @export
setMethod(numberImputations, "MImputedExpressionSets", function(object) {
   slot(object, "numberImputations")
})

##' @rdname numberImputations
##' @export
setReplaceMethod("numberImputations", "MImputedExpressionSets", function(object, value) {
   slot(object, "numberImputations") <- value
   validObject(object)
   return(object)
})

##' @name pData
##' @title pData
##' @importClassesFrom Biobase AnnotatedDataFrame
##' @export
setMethod("pData", "MImputedExpressionSets", function(object) {
   return(pData(object@data[[1]]))
})

##' @name fData
##' @title fData
##' @description TODO
##' @importClassesFrom Biobase AnnotatedDataFrame
##' @export
setMethod("fData", "MImputedExpressionSets", function(object) {
   return(fData(object@data[[1]]))
})

##' @name annotation
##' @title annotation
##' @importClassesFrom Biobase AnnotatedDataFrame
##' @export
setMethod("annotation", "MImputedExpressionSets", function(object) {
   return(annotation(object@data[[1]]))
})

##' @name eset
##' @title eset
##' @description desc
##' @importClassesFrom Biobase AnnotatedDataFrame
##' @export
setMethod("eset", "MImputedExpressionSets", function(object, imputation) {
   return(object@data[[imputation]])
})

##' @name intensities
##' @title intensities
##' @description TODO
##' @importClassesFrom Biobase AnnotatedDataFrame
##' @export
setMethod("intensities", "MImputedExpressionSets", function(object, imputation) {
   return(exprs(object@data[[imputation]]))
})

##' @name limmasFit
##' @title limmas fit
##' @description TODO
##' @importClassesFrom Biobase AnnotatedDataFrame
##' @export
setMethod("limmasFit", "MImputedExpressionSets", function(object, design) {
   return(new("MImputedMArrayLM", data=lapply(object@data, function(x) {
      # This is needed so Limma ignores NAs for building a model
      exprs(x)[is.na(exprs(x))] <- NaN
      fit  <- lmFit(x, design)
   })))
})

##' @title CompleteCase filter
##' @name completeCases
##' @return returns all complete rows
##'
##' @importClassesFrom Biobase ExpressionSet
##' @export
setMethod("completeCases", "MImputedExpressionSets", function(object) {
   return(filterRows(object, apply(exprs(object@data[[1]]), 1, function(x) !any(is.na(x)))))
})

##' @name filterRows
##' @title filterRows filter
##' @description TODO
##' @param filter vector for filtering
##' @return filter all rows
##' @importClassesFrom Biobase ExpressionSet
##' @export
setMethod("filterRows",  c(object = "MImputedExpressionSets",
                           filter = "logical"), function(object,
                                                         filter) {
   object@data  <- lapply(object@data, function(x) filterRows(x, filter))
   return(object)
})

##' @name filterCols
##' @title filter columns
##' @description filterCols filter
##' @param filter vector for filtering
##' @return filter all columns
##' @importClassesFrom Biobase ExpressionSet
##' @export
setMethod("filterCols",  c(object = "MImputedExpressionSets",
                           filter = "logical"), function(object,
                                                         filter) {
   object@data  <- lapply(object@data, function(x) filterCols(x, filter))
   return(object)
})

##' @name fillNAsWithValues
##' @title fill NAs with Values
##' @description fill missing values with values
##' @param value value to replace NAs
##' @return MIputedExpressionSets with replaced NAs
##' @export
setMethod("fillNAsWithValues", c(object = "MImputedExpressionSets",
                                 value = "numeric"), function(object,
                                                              value) {
   object@data  <- lapply(object@data, function(x) fillNAs(x, value))
   return(object)
})

##' @rdname plotExpression
##' @export
setMethod("plotExpression", c(object = "MImputedExpressionSets",
                              ID      = "character"        ), function(object,
                                                                       ID,
                                                                       ...) {
   id <- which(annotation(object) == ID)
   if (length(id) == 0) {
      stop(paste("id '", ID, "' not found", sep=""))
   } else if (length(id) > 1) {
      stop(paste(ID, "has more than 1 entries:", paste(id, sep=","), sep=" "))
   } else {
      plotExpression(object, id, ...)
   }
})

##' @name plotExpression
##' @title plot expression
##' @description function to plot a feature
##' @param ID character identification of feature ID
##' @param groupingCol grouping Col
##' @param colors colors
##' @param samplelabels labels for plot
##' @export
setMethod("plotExpression", c(object = "MImputedExpressionSets",
                              ID     = "numeric"), function(object,
                                                            ID,
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


##' @name getAverageExpression
##' @title get average expresion
##' @description get the average expression
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

##' @name [
##' @title [
##' @description selects subset of MImputedExpressionSets
##' @param x x
##' @param i i
##' @param j j
##' @param drop drop
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

##' @name calcGroupEstimations
##' @title calculate group estimations
##' @description calc group estimations
##' @export
setMethod("calcGroupEstimations", "MImputedExpressionSets", function(object) {
    # for each group, calculate the estimates
    object@groupData <- lapply(levels(pData(object@originalData)[,object@groupingCol]),
                               function(x) {
                                     getGroupData(object@originalData,
                                                  group = x,
                                                  groupCol = object@groupingCol)
                               })

})

##' @name estimateLimits
##' @title estimate limits
##' @description TODO
##' @param design design
##' @param contrasts contrasts
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




