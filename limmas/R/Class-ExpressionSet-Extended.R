# --------------------------------------------------------
# Extension of the Class ExpressionSet
# Authors: Thomas Schwarzl <schwarzl@embl.de> with help from Elisa D'Arcangelo
# Additional functions for ExpressionSet
# --------------------------------------------------------

##' @name filterRows
##' @title filter rows
##' @description Functions for filtering rows while maintaining data integrety
##' @param object ExpressionSet. input
##' @param filter vector. vector for filtering
##' @importClassesFrom Biobase ExpressionSet
##' @return filtered ExpressionSet
##' @export
setMethod("filterRows", c(object="ExpressionSet", filter="logical"), function(object, filter) {
   object <- object[filter,]
   fData <- fData(object)[filter,]
   annotation(object) <- annotation(object)[filter]
   return(object)
})

##' @name filterCols
##' @title filter columns
##' @description Functions for filtering columns while maintaining data integrety
##' @param object ExpressionSet. input
##' @param filter vector. vector for filtering
##' @importClassesFrom Biobase ExpressionSet
##' @return filtered ExpressionSet
##' @export
setMethod("filterCols", c(object="ExpressionSet", filter="logical"), function(object, filter) {
   object <- object[,filter]
   return(object)
})

##' @name checkMissingness
##' @title check missingness of samples in the data set
##' @description function for calculating percentage of missing values (NAs) per sample
##' @param minIntensity minimal intensity cutoff
##' @return percentage of missingness
##' @importClassesFrom Biobase ExpressionSet
##' @export
setMethod("checkMissingness", "ExpressionSet", function(data, minIntensity = 0) {
   intensities <- getIntensities(data, minIntensity=minIntensity)
   return(apply(intensities, 2, function(x) {
      sum(is.na(x))/length(x)*100
   }))
})

##' @name checkCompleteRows
##' @title checkCompleteRows
##' @description function for checking the number of complete rows
##' @param minIntensity minimal intensity cutoff
##' @importClassesFrom Biobase ExpressionSet
##' @return sum of complete rows
##' @export
setMethod("checkCompleteRows", "ExpressionSet", function(data, minIntensity = 0) {
   intensities <- getIntensities(data, minIntensity=minIntensity)
   return(sum(apply(intensities, 1, function(x) {
      sum(is.na(x)) == 0
   })))
})


##' @name getIntensities
##' @title get intensities
##' @description get intensities with above a certain min intensity cutoff
##' @param minIntensity minimal intensity cutoff
##' @importClassesFrom Biobase ExpressionSet
##' @return intensities
##' @export
setMethod("getIntensities", c(data="ExpressionSet", minIntensity="numeric"), function(data, minIntensity) {
   if(minIntensity < 0) {
      stop("minIntensity cannot be less than 0.")
   }

   intensities <- exprs(data)
   # intensities equal or less than minIntensity are set NA
   intensities[intensities <= minIntensity] <- NA

   return(intensities)
})

# -----------------------------------------------------------
# Filter
# -----------------------------------------------------------


##' @name peptideFilter
##' @title peptide filter
##' @description remove all probes where the peptide count is smaller than or equal to 1
##' get intensities with above a certain min intensity cutoff
##' @param peptideCutoff minimal peptide cutoff
##' @param peptideColumns columns having
##' @param method function which is used for filtering
##' @importClassesFrom Biobase ExpressionSet
##' @return filtered ExpresionSet
##' @export
setMethod("peptideFilter", "ExpressionSet", function (data,
                                                      peptideCutoff = 1,
                                                      peptideColumns = c("Peptides"),
                                                      method = c(any, all)) {
   if(!extends(class(data), "ExpressionSet")) {
      stop("'data' is not an object of class ExpressionSet or derived from an ExpressionSet")
   }

   filter <- NULL

   if(length(peptideColumns) == 1) {
      filter <- fData(data)[,peptideColumns] > peptideCutoff
   } else {
      filter <- unlist(apply(fData(data)[,peptideColumns], 1, function(x) {
                        method(x) > peptideCutoff
                }))
   }

   return(filterRows(data, filter))
})


##' @name reverseFilter
##' @title reverse filter
##' @description removes positive hits from reverse databases
##' @param symbol filter for symbol
##' @param reverseColumn name of column containing the reverse column
##' @importClassesFrom Biobase ExpressionSet
##' @return filtered ExpressionSet
##' @export
setMethod("reverseFilter", "ExpressionSet", function (data, symbol = "+", reverseColumn="Reverse") {
   if(!extends(class(data), "ExpressionSet")) {
      stop("'data' is not an object of class ExpressionSet or derived from an ExpressionSet")
   }

   if(!is.character(symbol)) {
      stop("'symbol' has to be a character symbol")
   }

   if(!is.character(reverseColumn)) {
      stop("'reverseColumn' has to be a character string")
   }

   rCol <- fData(data)[, reverseColumn]
   rCol[rCol == symbol] <- NA
   return(filterRows(data, !is.na(rCol)))
})

##' @name contaminantFilter
##' @title contaminant filter
##' @description removes features known to be contaminants
##' @param symbol filter for symbol
##' @param contaminantColumn name of column containing the contaminant column
##' @importClassesFrom Biobase ExpressionSet
##' @return filtered ExpressionSet
##' @export
setMethod("contaminantFilter", "ExpressionSet", function(data, symbol = "+", contaminantColumn="Contaminant") {
   if(!extends(class(data), "ExpressionSet")) {
      stop("'data' is not an object of class ExpressionSet or derived from an ExpressionSet")
   }

   if(!is.character(symbol)) {
      stop("'symbol' has to be a character symbol")
   }

   if(!is.character(contaminantColumn)) {
      stop("'contaminantColumn' has to be a character string")
   }

   cCol <- fData(data)[, contaminantColumn]
   cCol[cCol == symbol] <- NA

   return(filterRows(data, !is.na(cCol)))
})

##' @name completeCase
##' @title CompleteCases
##' @description returns complete cases
##' @importClassesFrom Biobase ExpressionSet
##' @return filtered ExpressionSet
##' @export
setMethod("completeCases", "ExpressionSet", function(object) {
   return(filterRows(object, apply(exprs(object), 1, function(x) !any(is.na(x)))))
})

# -----------------------------------------------------------
# Normalization
# -----------------------------------------------------------

##' @name normalizeData
##' @title normalized data
##' @description Function to normalize a dataset with a missing values minimum cutoff
##' @param minIntensity minimal intensity for NA cutoff, all values <= minIntensity are set to NA
##' @param FUN function for normalization of data
##' @import affyPLM
##' @return normalized ExpressionSet
##' @export
setMethod("normalizeData", "ExpressionSet", function(data, minIntensity=0, FUN=normalize.ExpressionSet.quantiles, ...) {
   if(!is.function(FUN)) {
      stop("'FUN' has to be a function")
   }

   exprs(data)[exprs(data) <= minIntensity] <- NA
   return(FUN(data, ...))
})

# -----------------------------------------------------------
# Transformation
# -----------------------------------------------------------

##' @name transformData
##' @title transforms data
##' @description Transforms transformation
##' @param FUN function for transformation
##' @importClassesFrom Biobase ExpressionSet
##' @return transformed ExpressionSet
##' @export
setMethod("transformData", "ExpressionSet",
          function(data,
                   FUN = function(x) { exprs(x) <- log(exprs(x)); return(x) }) {

   if(!is.function(FUN)) {
      stop("fun has to be a function")
   }

   return(FUN(data))
})



# -----------------------------------------------------------
# Scaling
# -----------------------------------------------------------

##' @name scaleData
##' @title scale data
##' @description This function scales samples X by given scalefactor: $\fraction{X}{scalefactor}$
##' @param scalefactor numeric. factor for scaling
##' @param FUN function. scaling operation
##' @importClassesFrom Biobase ExpressionSet
##' @return ExpressionSet. Scaled data
##' @export
setMethod("scaleData", "ExpressionSet",
          function(data, scalefactor = 1000,
                   FUN = function(x) { exprs(x) <- exprs(x) / scalefactor; return(x) }) {
             if(!is.numeric(scalefactor)) {
                stop("scalefactor has to be numeric")
             }
             if(!is.function(FUN)) {
                stop("fun has to be a function")
             }

             return(FUN(data))
          })


# -----------------------------------------------------------
# Detection limit
# -----------------------------------------------------------

##' @name getGroupData
##' @title get group data
##' @param group character string which specifies the group
##' @param groupCol name of the column specifying the groups in the pData slot of the ExpressionSet
##' @description This function calculates a number of statistics for NA relationship in the data
##' @importClassesFrom Biobase ExpressionSet
##' @return a list containing following elements:
##' groupData: expression data for the group
##' stats: data frame containing columns
##'        medianExpression: vector of median expression for the group
##'        sdExpression: standard deviation for expression
##'        naCount for expression
##' na.max: maximal NAs in sample
##' na.min: minimal NAs in samkple
##' na.vec: NA vector with all NA counts
##' na.count: NA count
##' i.max: max median expression
##' i.min: min median expression
##' @export
setMethod("getGroupData", "ExpressionSet",
      function(data, group, groupCol="groups", ...) {
         if(!is.character(groupCol)) {
            stop("'groupCol' is not a character string")
         }

         if(!groupCol %in% colnames(pData(data))) {
            stop(paste0("'groupCol' is set to ", groupCol, " which is not a column name in the pData slot of the ExpressionSet given."))
         }

         if(!group %in% pData(data)[,groupCol]) {
            stop(paste0("'group' is set to ", group, " which is not an element the column ", groupCol, " in the pData slot of the ExpressionSet given." ))
         }

         groupData <- exprs(data.transformed)[,rownames(pData(data)[pData(data.transformed)[,groupCol] == group,])]
         groupData <- groupData[!apply(groupData, 1, function(y) { all(is.na(y)) }),]
         medianExpression <- apply(groupData, 1, function(y) { median(y, na.rm = T)})
         sdExpression <- apply(groupData, 1, function(y) { sd(y, na.rm = T)})
         order     <- order(medianExpression)
         groupData <- groupData[order,]
         medianExpression <- medianExpression[order]
         sdExpression <- sdExpression[order]

         naCount   <- apply(groupData, 1, function(y) { sum(is.na(y)) })

         stats     <- data.frame(medianExpression = medianExpression,
                                 sdExpression     = sdExpression,
                                 naCount          = naCount)

         na.max    <- max(naCount)
         na.min    <- 0
         na.vec    <- na.min:na.max
         na.count  <- unlist(lapply(na.vec, function(z) { sum(naCount==z) }))
         i.max     <- max(groupData, na.rm = TRUE)
         i.min     <- min(groupData, na.rm = TRUE)

         ret <- list()
         ret[["groupData"]]        <- groupData
         ret[["stats"]]            <- stats
         ret[["na.max"]]           <- na.max
         ret[["na.min"]]           <- na.min
         ret[["na.vec"]]           <- na.vec
         ret[["na.count"]]         <- na.count
         ret[["i.max"]]            <- i.max
         ret[["i.min"]]            <- i.min
         return(ret)
})

##' @name plotsNA
##' @title plots for ExpressionSets with missing values
##' @param group group
##' @param groupCol group column
##' @param ... additional graphical parameters
##' @description
##' plotMedianVsNAs: plot the NA counts
##' plotMedianVsSD: Plots the relationship of median expression versus standard deviation
##' plotNAsVsSD: Plots the missing value counts versus the standard deviation
##' plotNAdensity: Plot the missing value (NA) density
##'
##' @importClassesFrom Biobase ExpressionSet
##' @export
setMethod("plotMedianVsNAs", "ExpressionSet",
   function(data, group, groupCol="groups", ...) {
      g <- getGroupData(data.transformed, group=group, groupCol=groupCol)

      boxplot(g[["stats"]][,"medianExpression"]~g[["stats"]][,"naCount"],
            main = paste(group, " (", ncol(g[["groupData"]]), " samples)", sep=""),
            xlab = "NA count",
            ylab = "median expression",
            type = "p",
            cex  = 0.2)
   })


##' @rdname plotsNA
##' @importClassesFrom Biobase ExpressionSet
##' @export
setMethod("plotMedianVsSD", "ExpressionSet",
          function(data, group, groupCol="groups", ...) {
             g <- getGroupData(data.transformed, group=group, groupCol=groupCol)

             plot(g[["stats"]][,"medianExpression"], g[["stats"]][,"sdExpression"],
                  main = paste(group, " (", ncol(g[["groupData"]]), " samples)", sep=""),
                  ylab = "standard deviation",
                  xlab = "median expression",
                  type = "p",
                  cex  = 0.2,
                  ...)
          })

##' @rdname plotsNA
##' @importClassesFrom Biobase ExpressionSet
##' @export
setMethod("plotNAsVsSD", "ExpressionSet",
          function(data, group, groupCol="groups", ...) {
             g <- getGroupData(data.transformed, group=group, groupCol=groupCol)

            boxplot(g[["stats"]][,"sdExpression"]~g[["stats"]][,"naCount"],
                    main = paste(group, " (", ncol(g[["groupData"]]), " samples)", sep=""),
                    ylab = "standard deviation",
                    xlab = "NA count",
                    type = "p",
                    cex  = 0.2,
                    ...)
          })


##' @rdname plotsNA
##' @importClassesFrom Biobase ExpressionSet
##' @export
setMethod("plotNAdensity", "ExpressionSet", function(data, group, groupCol="groups", ...) {
   # get preprocessed data for the groups
   g <- getGroupData(data, group=group, groupCol=groupCol)
   # create and store the density plots for printing together
   plots <- lapply(g[["na.vec"]], function(z) {
      density(g[["stats"]][,"medianExpression"][g[["stats"]][,"naCount"] == as.numeric(z)])
   })

   # determine boundries of the plot
   d.max <- max(unlist(lapply(plots, function(z) { max(z$y) })))
   # and the colors
   colvec <- rainbow(g[["na.max"]] + 1)

   # plot all together
   plot(plots[[1]], col=colvec[1], lwd=3,
        main = paste("# NAs for group ", group, sep=""),
        xlim = c(g[["i.min"]], g[["i.max"]]), ylim = c(0, d.max), ...)
   zz <- lapply(1:(g[["na.max"]]), function(z) { lines(plots[[z+1]], col=colvec[z + 1], lwd=3) })
   legend("topright", "# NAs", paste(g[["na.vec"]], " (", g[["na.count"]], " features)", sep=""), fill = colvec,  )
})


##' @name plotChangeMissing
##' @title plots chance of missingness
##' @param group group
##' @param groupCol group column
##' @param windowSize size of window
##' @param stepSize size of the steps to go
##' @description Plot the chance of a value being missing for a certain median expression
##' @return ggplot2 plot
##' @importClassesFrom Biobase ExpressionSet
##' @export
setMethod("plotChangeMissing", "ExpressionSet", function(data, group, groupCol="groups", windowSize = 0.3, stepSize = 0.5) {
   g <- getGroupData(data.transformed, group=group, groupCol=groupCol)

   windowSize = 0.5
   stepSize = 0.1

   steps = seq(g[["i.min"]], g[["i.max"]], stepSize)
   perc  = unlist(lapply(steps, function(i) {
      j <- i + windowSize

      window <- g[["stats"]][,"medianExpression"] >= i & g[["stats"]][,"medianExpression"] <= j

      ret <- sum(g[["stats"]][,"naCount"][ window ]) / (g[["na.max"]] * sum(window))
      if(is.nan(ret)) {
         return(0)
      } else {
         return(ret)
      }
   }))

   X <- data.frame(steps, perc)


   ggplot(X, aes(x = steps, y=perc)) +
      geom_point(size = I(3)) +
      geom_line() +
      theme_minimal() +
      theme(text = element_text(size=20)) +
      xlim(c(g[["i.min"]], g[["i.max"]])) +
      xlab("median expression") +
      ylab("change of missingness") +
      ggtitle(paste0("Chance of missingness, window size: ", windowSize))
})

# -----------------------------------------------------------
# imputation
# -----------------------------------------------------------

##' @name imputeIndependentGroupsWithAmelia
##' @title impute independent groups with amelia
##' @description impute independent groups with amelia
##'
##' @param minPresent
##' @param m number of imputations (default: 10, usually between 3 and 10)
##' @param groupingCol character string specifying the column containing group information
##' @param quiet logical. if TRUE it will suppress the messages from Amelia
##' @param ... parameters will be forwarded to the amelia function from package Amelia
##' check out their specific parameters and options for parallelizing
##' @importClassesFrom Biobase ExpressionSet
##' @return MImputedExpressionSets
##' @export
setMethod("imputeIndependentGroupsWithAmelia", "ExpressionSet", function(data.input,
                                                                         minPresent = 0.5,
                                                                         groupingCol = "groups",
                                                                         m = 10,
                                                                         quiet = F,
                                                                         ...) {
   suppressPackageStartupMessages(require(Amelia))

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
      if(quiet) {
         return(suppressMessages(amelia(x, m=m, ...)))
      } else {
         return(amelia(x, m=m, ...))
      }
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

##' @name fillNAs
##' @title fill missing values (NAs)
##' @param value value to replace NA
##' @description NA values will be replaced by given 'value'
##' @return complete data
##' @export
setMethod("fillNAs", c(object="ExpressionSet", value="numeric"), function(object, value) {
   tab <- exprs(object)
   tab[is.na(tab)] <- value
   exprs(object)  <- tab
   return(object)
})

##' @title calculate feature correlations
##' @name calculateFeatureCorrelations
##' @param use parameter passed on to \code{cor} function
##' @param method method passed on to \code{cor} function
##' @description this function will calculate a correlation matrix for the samples with given function.
##' by default it will calculate the correlations using the complete observations.
##' @return feature
##' @export
setMethod("calculateFeatureCorrelations", "ExpressionSet", function(object,
                                                                    use = "pairwise.complete.obs",
                                                                    method = "pearson") {
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


