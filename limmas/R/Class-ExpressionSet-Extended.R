# Extend the Class ExpressionSet, allowing the functionality to be the same


##' Functions for filtering rows while maintaining data integrety
##'
##' @param filter vector for filtering
##' @importClassesFrom Biobase ExpressionSet
##' @export
setMethod("filterRows", c(object="ExpressionSet", filter="logical"), function(object, filter) {
   object <- object[filter,]
   fData <- fData(object)[filter,]
   annotation(object) <- annotation(object)[filter]
   return(object)
})

##' Functions for filtering columns while maintaining data integrety
##'
##' @param filter vector for filtering
##' @importClassesFrom Biobase ExpressionSet
##' @return filtered ExpressionSet
##' @export
setMethod("filterCols", c(object="ExpressionSet", filter="logical"), function(object, filter) {
   object <- object[,filter]
   return(object)
})


##' function for calculating percentage of NAs per sample
##'
##' @param minIntensity minimal intensity cutoff
##' @return  percentage of missingness
##' @importClassesFrom Biobase ExpressionSet
##' @export
setMethod("checkMissingness", "ExpressionSet", function(data, minIntensity=0){
   intensities <- getIntensities(data, minIntensity=minIntensity)
   return(apply(intensities, 2, function(x) {
      sum(is.na(x))/length(x)*100
   }))
})

##' function for checking the number of complete rows
##'
##' @param minIntensity minimal intensity cutoff
##' @importClassesFrom Biobase ExpressionSet
##' @return sum of complete rows
##' @export
setMethod("checkCompleteRows", "ExpressionSet", function(data, minIntensity=0) {
   intensities <- getIntensities(data, minIntensity=minIntensity)
   return(sum(apply(intensities, 1, function(x) {
      sum(is.na(x)) == 0
   })))
})


##' get intensities with above a certain min intensity cutoff
##'
##' @param minIntensity minimal intensity cutoff
##' @importClassesFrom Biobase ExpressionSet
##' @return intensities
##' @export
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


##' peptide filter
##'
##' remove all probes where the peptide count is smaller than or equal to 1
##' get intensities with above a certain min intensity cutoff
##'
##' @param peptideCutoff minimal peptide cutoff
##' @param peptideColumns columns having
##' @param method function which is used for filtering
##' @importClassesFrom Biobase ExpressionSet
##' @return filtered ExpresionSet
##' @export
setMethod("peptideFilter", "ExpressionSet", function (data, peptideCutoff=1, peptideColumns=c("Peptides"), method=c(any, all)) {
   if(!class(data)[1] == "ExpressionSet") {
      stop("data is not an object of class ExpressionSet")
   }

   filter <- NULL

   if(length(peptideColumns) == 1) {
      filter <- fData(data)[,peptideColumns] > peptideCutoff
   } else {
      filter <- unlist(apply(fData(data)[,peptideColumns], 1, function(x) { method(x) > peptideCutoff }))
   }

   return(filterRows(data, filter))
})


##' reverse filter
##'
##' removes positive hits from reverse databases
##'
##' @param symbol filter for symbol
##' @param reverseColumn name of column containing the reverse column
##' @importClassesFrom Biobase ExpressionSet
##' @return filtered ExpressionSet
##' @export
setMethod("reverseFilter", "ExpressionSet", function (data, symbol = "+", reverseColumn="Reverse") {
   rCol <- fData(data)[, reverseColumn]
   rCol[rCol == symbol] <- NA
   return(filterRows(data, !is.na(rCol)))
})

##' contaminant filter
##'
##' removes features known to be contaminantss
##'
##' @param symbol filter for symbol
##' @param contaminantColumn name of column containing the contaminant column
##' @importClassesFrom Biobase ExpressionSet
##' @return filtered ExpressionSet
##' @export
setMethod("contaminantFilter", "ExpressionSet", function(data, symbol = "+", contaminantColumn="Contaminant") {
   cCol <- fData(data)[, contaminantColumn]
   cCol[cCol == symbol] <- NA

   return(filterRows(data, !is.na(cCol)))
})


##' returns complete cases
##'
##' @importClassesFrom Biobase ExpressionSet
##' @return filtered ExpressionSet
##' @export
setMethod("completeCases", "ExpressionSet", function(object) {
   return(filterRows(object, apply(exprs(object), 1, function(x) !any(is.na(x)))))
})

# -----------------------------------------------------------
# Normalization
# -----------------------------------------------------------

##' normalized data
##'
##' @param minIntensity minimal intensity for NA cutoff
##' @param FUN function for normalization
##' @import affyPLM
##' @return normalized ExpressionSet
##' @export
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

##' transforms data
##'
##' @param FUN function for transformation
##' @importClassesFrom Biobase ExpressionSet
##' @return transformed ExpressionSet
##' @export
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

##' scale data
##'
##' @param scalefactor factor for scaling
##' @param FUN function for scaling
##' @importClassesFrom Biobase ExpressionSet
##' @return scaled ExpressionSet
##' @export
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


#TODO
setMethod("getGroupData", "ExpressionSet",
      function(data, group, groupCol="groups", ...) {
         groupData <- exprs(data.transformed)[,rownames(pData(data)[pData(data.transformed)[,groupCol] == group,])]
         groupData <- groupData[!apply(groupData, 1, function(y) { all(is.na(y)) }),]
         medianExpression <- apply(groupData, 1, function(y) { median(y, na.rm = T)})
         sdExpression <- apply(groupData, 1, function(y) { sd(y, na.rm = T)})
         order     <- order(medianExpression)
         groupData <- groupData[order,]
         medianExpression <- medianExpression[order]
         sdExpression <- sdExpression[order]
         naCount   <- apply(groupData, 1, function(y) { sum(is.na(y)) })
         na.max    <- max(naCount)
         na.min    <- 0
         na.vec    <- na.min:na.max
         na.count  <- unlist(lapply(na.vec, function(z) { sum(naCount==z) }))
         i.max     <- max(groupData, na.rm = TRUE)
         i.min     <- min(groupData, na.rm = TRUE)

         ret <- list()
         ret[["groupData"]]        <- groupData
         ret[["medianExpression"]] <- medianExpression
         ret[["sdExpression"]]     <- sdExpression
         ret[["naCount"]]          <- naCount
         ret[["na.max"]]           <- na.max
         ret[["na.min"]]           <- na.min
         ret[["na.vec"]]           <- na.vec
         ret[["na.count"]]         <- na.count
         ret[["i.max"]]            <- i.max
         ret[["i.min"]]            <- i.min
         return(ret)
})


# plot the NA counts
setMethod("plotMedianVsNAs", "ExpressionSet",
   function(data, group, groupCol="groups", ...) {
      g <- getGroupData(data.transformed, group=group, groupCol=groupCol)

      plot(g[["medianExpression"]], g[["naCount"]],
      main = paste(group, " (", ncol(g[["groupData"]]), " samples)", sep=""),
      ylab = "NA count",
      xlab = "median expression",
      type = "p",
      cex  = 0.2)
   })



setMethod("plotMedianVsSD", "ExpressionSet",
          function(data, group, groupCol="groups", ...) {
             g <- getGroupData(data.transformed, group=group, groupCol=groupCol)

             plot(g[["medianExpression"]], g[["sdExpression"]],
                  main = paste(group, " (", ncol(g[["groupData"]]), " samples)", sep=""),
                  ylab = "standard deviation",
                  xlab = "median expression",
                  type = "p",
                  cex  = 0.2)
          })


setMethod("plotNAsVsSD", "ExpressionSet",
          function(data, group, groupCol="groups", ...) {
             g <- getGroupData(data.transformed, group=group, groupCol=groupCol)

             plot(g[["naCount"]], g[["sdExpression"]],
                  main = paste(group, " (", ncol(g[["groupData"]]), " samples)", sep=""),
                  ylab = "standard deviation",
                  xlab = "NA count",
                  type = "p",
                  cex  = 0.2)
          })



setMethod("plotNAdensity", "ExpressionSet", function(data, group, groupCol="groups", ...) {
   # get preprocessed data for the groups
   g <- getGroupData(data, group=group, groupCol=groupCol)
   # create and store the density plots for printing together
   plots <- lapply(g[["na.vec"]], function(z) { density(g[["medianExpression"]][g[["naCount"]] == as.numeric(z)]) })
   # determine boundries of the plot
   d.max <- max(unlist(lapply(plots, function(z) { max(z$y) })))
   # and the colors
   colvec <- rainbow(g[["na.max"]] + 1)

   # plot all together
   plot(plots[[1]], col=colvec[1], lwd=3,
        main = paste("# NAs for group ", group, sep=""),
        xlim = c(g[["i.min"]], g[["i.max"]]), ylim = c(0, d.max))
   zz <- lapply(1:(g[["na.max"]]), function(z) { lines(plots[[z+1]], col=colvec[z + 1], lwd=3) })
   legend("topright", "# NAs", paste(g[["na.vec"]], " (", g[["na.count"]], " features)", sep=""), fill = colvec,  )
})



# g <- getGroupData(data.transformed, group=group, groupCol=groupCol)
#
# names(g)
#
# windowSize = 0.3
# steps = seq(g[["i.min"]], g[["i.max"]], 0.5)
# perc  = unlist(lapply(steps, function(i) {
#    j <- i + windowSize
#
#    window <- g[["medianExpression"]] >= i & g[["medianExpression"]] <= j
#
#    g[["medianExpression"]][ window ]
#    ret <- sum(g[["naCount"]][ window ]) / (g[["na.max"]] * sum(window))
#    if(is.nan(ret)) {
#       return(0)
#    } else {
#       return(ret)
#    }
# }))
#
# plot(steps, perc, type="b")
#
# plot(cumsum(perc))

# -----------------------------------------------------------
# imputation
# -----------------------------------------------------------

##' impute independent groups with amelia
##'
##' @param scalefactor factor for scaling
##' @param FUN function for scaling
##' @importClassesFrom Biobase ExpressionSet
##' @return scaled ExpressionSet
##' @export
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


# NA values will be replaced by given 'value'
setMethod("fillNAs", c(object="ExpressionSet", value="numeric"), function(object, value) {
   tab <- exprs(object)
   tab[is.na(tab)] <- value
   exprs(object)  <- tab
   return(object)
})

# this function will calculate a correlation matrix for the samples with given function.
# by default it will calculate the correlations using the complete observations.
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


