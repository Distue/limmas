# Extend the Class ExpressionSet, allowing the functionality to be the same


# ------------------------------------------
# public functions
# ------------------------------------------

# Missingness
# function for calculating percentage of NAs per sample
# data is ExpressionSet
setGeneric("checkMissingness", function(data, ...) standardGeneric("checkMissingness"))
setMethod("checkMissingness", "ExpressionSet", function(data, minIntensity=0) {
   intensities <- getIntensities(data, minIntensity=minIntensity)
   return(apply(intensities, 2, function(x) {
      sum(is.na(x))/length(x)*100
   }))
})

setGeneric("getMaxProbabilityMissingByChance", function(data, groupingCol, numberMissing, ...) standardGeneric("getMaxProbabilityMissingByChance"))
setMethod("getMaxProbabilityMissingByChance", c(data="ExpressionSet", groupingCol="character", numberMissing="numeric"), function(data, groupingCol, numberMissing, ...) {
   if(numberMissing < 1) {
      stop("numberMissing has to be equal or greater than 1")
   }
   if (!any(groupingCol %in% colnames(pData(data)))) {
      stop("groupingCol is not a column name of pData(data)")
   }
   
   missingness <- checkMissingness(data, ...)/100
   
   groups <- as.character(pData(data)[, groupingCol])
   names(groups) <- rownames(pData(data))
   
   groupprobability <- lapply(unique(groups), function(x) {
      prod(sort(missingness[names(groups[groups == x])], decreasing=T)[1:numberMissing])
   })
   
   names(groupprobability) <- unique(groups)
   return(groupprobability)
})

setGeneric("checkMissingness", function(data, ...) standardGeneric("checkMissingness"))
setMethod("checkMissingness", "ExpressionSet", function(data, minIntensity=0){
   intensities <- getIntensities(data, minIntensity=minIntensity)
   return(apply(intensities, 2, function(x) {
      sum(is.na(x))/length(x)*100
   }))
})


# function for checking the number of complete rows
setGeneric("checkCompleteRows", function(data, ...) standardGeneric("checkCompleteRows"))
setMethod("checkCompleteRows", "ExpressionSet", function(data, minIntensity=0) {
   intensities <- getIntensities(data, minIntensity=minIntensity)
   return(sum(apply(intensities, 1, function(x) {
      sum(is.na(x)) == 0
   })))
})

# ------------------------------------------
# private functions
# ------------------------------------------
setGeneric("getIntensities", function(data, minIntensity) standardGeneric("getIntensities"))
setMethod("getIntensities", c(data="ExpressionSet", minIntensity="numeric"), function(data, minIntensity) {
   if(minIntensity < 0) {
      stop("minIntensity cannot be smaller than 0.") 
   }

   intensities <- exprs(data)
   # intensities <= minIntensity are set NA
   intensities[intensities <= minIntensity] <- NA
   
   return(intensities)
})

# -----------------------------------------------------------
# Filter
# -----------------------------------------------------------
setGeneric("filterRows", function(object, filter) standardGeneric("filterRows"))
setMethod("filterRows", "ExpressionSet", function(object, filter) {
   object <- object[filter,]
   fData <- fData(object)[filter,]
   annotation(object) <- annotation(object)[filter]
   return(object)
})

setGeneric("filterCols", function(object, filter) standardGeneric("filterCols"))
setMethod("filterCols", "ExpressionSet", function(object, filter) {
   object <- object[,filter]

   return(object)
})


#Peptide Filter: remove all probes where the peptide count is smaller than or equal to 1
setGeneric("peptideFilter", function(data, ...) standardGeneric("peptideFilter"))
setMethod("peptideFilter", "ExpressionSet", function (data, peptideCutoff=1, peptideColumn="Peptides") {
   if(!class(data)[1] == "ExpressionSet") {
      stop("data is not an object of class ExpressionSet")
   }  
   
   return(filterRows(data, fData(data)[,peptideColumn] > peptideCutoff))
})


# Reverse Filter: removes positive hits from reverse databases
setGeneric("reverseFilter", function(data, ...) standardGeneric("reverseFilter"))
setMethod("reverseFilter", "ExpressionSet", function (data, symbol = "+", reverseColumn="Reverse") {
   rCol <- fData(data)[, reverseColumn]
   rCol[rCol == symbol] <- NA
   return(filterRows(data, !is.na(rCol)))
})

#Contaminant Filter: removes features known to be contaminants
setGeneric("contaminantFilter", function(data, ...) standardGeneric("contaminantFilter"))
setMethod("contaminantFilter", "ExpressionSet", function(data, symbol = "+", contaminantColumn="Contaminant") {
   cCol <- fData(data)[, contaminantColumn]
   cCol[cCol == symbol] <- NA
   
   return(filterRows(data, !is.na(cCol)))
})

setGeneric("completeCases", function(object) standardGeneric("completeCases"))
setMethod("completeCases", "ExpressionSet", function(object) {
   return(filterRows(object, apply(exprs(object), 1, function(x) !any(is.na(x)))))
})

# -----------------------------------------------------------
# Normalization 
# -----------------------------------------------------------

#quantile normalisation
setGeneric("normalizeData", function(data, ...) standardGeneric("normalizeData"))
setMethod("normalizeData", "ExpressionSet", function(data, minIntensity=0, FUN=normalize.ExpressionSet.quantiles, ...) {
   if(!is.function(FUN)) {
      stop("fun has to be a function")
   }
   
   exprs(data)[exprs(data) <= minIntensity] <- NA
   return(FUN(data, ...))
})

# -----------------------------------------------------------
# Transformation 
# -----------------------------------------------------------
setGeneric("transformData", function(data, ...) standardGeneric("transformData"))
setMethod("transformData", "ExpressionSet",
          function(data,
                   FUN = function(x) { exprs(x) <- log(exprs(x)); return(x) }){

   if(!is.function(FUN)) {
      stop("fun has to be a function")
   }
   
   return(FUN(data))
})

# -----------------------------------------------------------
# Scaling 
# -----------------------------------------------------------
setGeneric("scaleData", function(data, scalefactor, ...) standardGeneric("scaleData"))
setMethod("scaleData", "ExpressionSet",
          function(data, scalefactor = 1000, 
                   FUN = function(x) { exprs(x) <- exprs(x) / scalefactor; return(x) }){
             
             if(!is.function(FUN)) {
                stop("fun has to be a function")
             }
             
             return(FUN(data))
          })


# -----------------------------------------------------------
# imputation
# -----------------------------------------------------------

setGeneric("imputeIndependentGroupsWithAmelia", function(data.input, ...) standardGeneric("imputeIndependentGroupsWithAmelia"))
setMethod("imputeIndependentGroupsWithAmelia", "ExpressionSet", function(data.input, minPercentPresent=0, minTotalPresent=1, groupingCol="groups", m=10, empri=.01 * nrow(data.input), ...) {   
   #prepare for Amelia
   extractGroups <- function (data, groups) {
      list_for_amelia <- list()
      for (i in unique(groups)) { 
         list_for_amelia[[i]]  <- exprs(data)[,names(groups)[groups == i]]
      }
      return(list_for_amelia)
   }
   
   
   #for each row: where there are less than 'minPresent' of observed values, 
   #all values are replaced by NA; where there are more than or exactely 'minPresent', 
   #all values and NAs remain unaltered
   correctFalsePositives <- function (groupTables, minPercentPresent, minTotalPresent=0) { 
      return(lapply(groupTables, function(tab) {
         t(apply(tab, 1, function(x) {
            # if the sum of values is larger than amount of minimal present, then set everything NA
            if(sum(!is.na(x)) < max(length(x) * minPercentPresent, minTotalPresent)) {
               x[1:length(x)] <- NA 
            }
            return(x)
         }))  
      }))
   } 
   
   compileImputatedGroups <- function(imputedGroups, m, pheno, features, annotation) {
      imputationList <- lapply(1:m, function(i) {
         compiledTable <- do.call(cbind, lapply(imputedGroups, function(x) { 
            colN <- colnames(x)
            return(x[,grep(paste("imp", i, "\\.", sep=""), colN)])
         }))
         
         colnames(compiledTable) <- sub("^.+?\\.imp[0-9]+?\\.", "", colnames(compiledTable))

         if (!all(colnames(compiledTable) %in% rownames(pData(pheno)))) {
            stop("Problems when compiling the imputed groups: trimming of colnames failed.")
         }
         
         compiledTable <- compiledTable[,rownames(pData(pheno))]
         return(createExpressionSet(compiledTable, pheno, features, annotation))
      })
   }
   
   
   if(!is(data.input, "ExpressionSet")) {
      stop("data is not an object of class ExpressionSet")
   }
   if (!is.numeric(minPercentPresent)) {
      stop("minPercentPresent has to be numeric")
   }
   if (!is.numeric(minTotalPresent)) {
      stop("minTotalPresent has to be numeric")
   }
   if (!(minPercentPresent >= 0 && minPercentPresent < 1)) {
      stop("minPercentPresent has to be between 0 and 1")
   }
   if (!(minTotalPresent > 0)) {
      stop("minTotalPresent has to be bigger than 0")
   }
   if (!is.numeric(m)) {
      stop("m has to be numeric")
   }
   if(!m > 2) {
      stop("m hast to be greater than 2")
   }
   if(!groupingCol %in% colnames(pData(data.input))) {
      stop("groupingCol has to be a column name of pData(data.input)")
   }
   
   # create groups list
   groups <- as.character(pData(data.input)[, groupingCol])
   names(groups) <- rownames(pData(data.input))
   
   # prepare list for amelia
   groupTables <- extractGroups(data.input, groups)
   
   # false positive filter
   groupTables <- correctFalsePositives(groupTables, minPercentPresent=minPercentPresent, minTotalPresent=minTotalPresent)
   
   imputations <- lapply(groupTables, function(x) { 
      return(amelia(x, m=m, empri=empri, ...))
   })
   
   #impute with amelia
   imputedGroups <- lapply(imputations, function(x) { 
      return(as.data.frame(x$imputations))
   })
   
   # create list of expression sets
   allImputations <- compileImputatedGroups(imputedGroups, m, pheno=phenoData(data.input), features=featureData(data.input), annotation=annotation(data.input))
   
   imp.output <- new("MImputedExpressionSets", data=allImputations, minPercentPresent = minPercentPresent, minTotalPresent=minTotalPresent, groupingCol = groupingCol, numberImputations = m)
   
   returnList <- list()
   returnList[["data"]] <- imp.output
   returnList[["imputation"]] <- imputations
   return(returnList)
})



setGeneric("imputeTimeWithAmelia", function(data.input, ...) standardGeneric("imputeTimeWithAmelia"))
setMethod("imputeTimeWithAmelia", "ExpressionSet", function(data.input, timeCol="time", minTotalPresent=ncol(data.input)/2, m=10, empri=.01 * nrow(data.input), ...) {   
   if (!is.numeric(m)) {
      stop("m has to be numeric")
   }
   if(!m > 2) {
      stop("m hast to be greater than 2")
   }
   if(!timeCol %in% colnames(pData(data.input))) {
      stop("timeCol has to be a column name of pData(data.input)")
   }
   
   time <- pData(data.input)[,timeCol]
   
   dataTable <- t(exprs(data.input)[,names(time)])
   
   # filter out observations with not more than minTotalPresent observations
   # filter out any columns which do not vary
   filter <- apply(dataTable, 2, function(x) {
      return((sum(is.na(x)) <= (length(x) - minTotalPresent)) && length(unique(x)) > 1)
   })
   
   dataTable <- dataTable[,filter]
       
   dataTable <- cbind(as.numeric(as.character(time)), dataTable)
   colnames(dataTable)[1] <- ".time"
   

   # filter 
   imputations <- amelia(dataTable, m=m, empri=empri, ts=".time") 
   
   # step, remove time
   
   #impute with amelia
   imputedGroups <- lapply(imputations, function(x) { 
      return(t(as.data.frame(x$imputations)))
   })
   
   # create list of expression sets
   allImputations <- compileImputatedGroups(imputedGroups, m, pheno=phenoData(data.input), features=featureData(data.input), annotation=annotation(data.input))
   
   imp.output <- new("MImputedExpressionSets", data=allImputations, minPercentPresent = minPercentPresent, minTotalPresent=minTotalPresent, groupingCol = groupingCol, numberImputations = m)
   
   return(c(data=imp.output, imputations=imputations))
})


# -----------------------------------------------------------
# Manipulation
# -----------------------------------------------------------

setGeneric("fillNAs", function(object, value) standardGeneric("fillNAs"))
setMethod("fillNAs", c(object="ExpressionSet", value="numeric"), function(object, value) { 
   tab <- exprs(object)
   tab[is.na(tab)] <- value
   exprs(object)  <- tab
   return(object)
})

setGeneric("calculateFeatureCorrelations",  function(object, ...) standardGeneric("calculateFeatureCorrelations"))
setMethod("calculateFeatureCorrelations", "ExpressionSet", function(object, use="pairwise.complete.obs", method="pearson") {
   # to cope with the standard deviation is zero problem
   filterForCor <- function(input) {
      non.zero.var <- logical()
      for(i in 1:ncol(input)) {
         non.zero.var[i] <- var(input[,i]) > 0
      }
      return(input[,non.zero.var])
   }
   
   
   return(cor(filterForCor(t(as.matrix(exprs(object)))), use=use, method=method))
})

