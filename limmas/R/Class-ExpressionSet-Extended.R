# Extend the Class ExpressionSet, allowing the functionality to be the same


# ------------------------------------------
# public functions
# ------------------------------------------

# Missingness
# function for calculating percentage of NAs per sample
# data is ExpressionSet
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
# shift 
# -----------------------------------------------------------

setGeneric("shiftData", function(data, scalefactor, ...) standardGeneric("shiftData"))
setMethod("shiftData", "ExpressionSet",
          function(data, scalefactor = 1000,
                   FUN = function(x) { exprs(x) <- exprs(x) / scalefactor; return(x) }){
             
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
# Detection limit
# -----------------------------------------------------------


setGeneric(".getGroupData", function(data, group, ...) standardGeneric(".getGroupData"))
setMethod(".getGroupData", "ExpressionSet",
      function(data, group, groupCol="groups", ...) { 
         groupData <- exprs(data.transformed)[,rownames(pData(data)[pData(data.transformed)[,groupCol] == group,])]
         groupData <- groupData[!apply(groupData, 1, function(y) { all(is.na(y)) }),]
         medianExpression <- apply(groupData, 1, function(y) { median(y, na.rm = T)})
         order <- order(medianExpression)
         groupData <- groupData[order,]
         medianExpression <- medianExpression[order]
         naCount <- apply(groupData, 1, function(y) { sum(is.na(y)) })
         na.max <- max(naCount)
         na.min <- 0
         na.vec <- na.min:na.max
         na.count <- unlist(lapply(na.vec, function(z) { sum(naCount==z) }))
         i.max <- max(groupData, na.rm = TRUE)
         i.min <- min(groupData, na.rm = TRUE)
         
         ret <- list()
         ret[["groupData"]] <- groupData
         ret[["medianExpression"]] <- medianExpression
         ret[["naCount"]] <- naCount
         ret[["na.max"]] <- na.max
         ret[["na.min"]] <- na.min
         ret[["na.vec"]] <- na.vec
         ret[["na.count"]] <- na.count
         ret[["i.max"]] <- i.max
         ret[["i.min"]] <- i.min
         return(ret)
})

setGeneric("plotNAs", function(data, group, ...) standardGeneric("plotNAs"))
setMethod("plotNAs", "ExpressionSet",
          function(data, group, groupCol="groups", ...) {  
             g <- .getGroupData(data.transformed, group=group, groupCol=groupCol)
             
         plot(g[["medianExpression"]], g[["naCount"]],
           main = paste(group, " (", ncol(g[["groupData"]]), " samples)", sep=""),
           ylab = "NA count",
           xlab = "median expression",
           type = "p",
           cex  = 0.2)
          })


setGeneric("plotNAhistogram", function(data, group, ...) standardGeneric("plotNAhistogram"))
setMethod("plotNAhistogram", "ExpressionSet",
          function(data, group, groupCol="groups", ...) {     
         g <- .getGroupData(data.transformed, group=group, groupCol=groupCol)
         plots <- lapply(g[["na.vec"]], function(z) { density(g[["medianExpression"]][g[["naCount"]] == as.numeric(z)]) }) 
         d.max <- max(unlist(lapply(plots, function(z) { max(z$y) })))
         colvec <- rainbow(g[["na.max"]] + 1)
         plot(plots[[1]], col=colvec[1], lwd=3,
              main = paste("# NAs for group ", group, sep=""),
              xlim = c(g[["i.min"]], g[["i.max"]]), ylim = c(0, d.max))
         zz <- lapply(1:(g[["na.max"]]), function(z) { lines(plots[[z+1]], col=colvec[z + 1], lwd=3) })
         legend("topright", "# NAs", paste(g[["na.vec"]], " (", g[["na.count"]], " features)", sep=""), fill = colvec,  )
})

# -----------------------------------------------------------
# imputation
# -----------------------------------------------------------

setGeneric("imputeIndependentGroupsWithAmelia", function(data.input, ...) standardGeneric("imputeIndependentGroupsWithAmelia"))
setMethod("imputeIndependentGroupsWithAmelia", "ExpressionSet", function(data.input, minPresent=0.5, groupingCol="groups", m=10,  ...) {   
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
   correctFalsePositives <- function (groupTables, minPresent) { 
      return(lapply(groupTables, function(tab) {
         t(apply(tab, 1, function(x) {
            if(sum(is.na(x)) > length(x) * minPresent) {
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
         
         colnames(compiledTable) <- sub("^.+?\\.imp.+?\\.", "", colnames(compiledTable))
         
         compiledTable <- compiledTable[,rownames(pData(pheno))]
         return(createExpressionSet(compiledTable, pheno, features, annotation))
      })
   }
   
   
   if(!is(data.input, "ExpressionSet")) {
      stop("data is not an object of class ExpressionSet")
   }
   if (!is.numeric(minPresent)) {
      stop("minPresent has to be numeric")
   }
   if (! (minPresent > 0 && minPresent < 1)) {
      stop("minPresent has to be between 0 and 1")
   }
   if (!is.numeric(m)) {
      stop("m has to be numeric")
   }
   if(! m > 2) {
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
   groupTables <- correctFalsePositives(groupTables, minPresent=minPresent)   
   
   imputations <- lapply(groupTables, function(x) {
      return(amelia(x, m=m, ...))
   })
   
   ##impute with amelia
   imputedGroups <- lapply(imputations, function(x) {
      return(as.data.frame(x$imputations))
   })
   
   # create list of expression sets
   allImputations <- compileImputatedGroups(imputedGroups, m, phenoData(data.input), featureData(data.input), annotation(data.input))

   imp.output <- new("MImputedExpressionSets", data=allImputations, minPresent = minPresent,
                     groupingCol = groupingCol, numberImputations = m, originalData = data.input)
   
   returnList <- list()
   returnList[["data"]] <- imp.output
   returnList[["imputation"]] <- imputations
   return(returnList)
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


