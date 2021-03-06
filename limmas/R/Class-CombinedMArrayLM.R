# --------------------------------------------------------
# Class CombinedMArrayLM
# Authors: Thomas Schwarzl <schwarzl@embl.de> with help from Elisa D'Arcangelo
# Combined MArrayLM from multiple imputation data
# --------------------------------------------------------
# Constructor:
##' @rdname MImputedExpressionSets
##' @export
CombinedMArrayLM <- function(ids,
                             coefficients,
                             tstat,
                             p.value,
                             featurelist,
                              ...) {
   return(new(Class                = "CombinedMArrayLM",
              ids            = ids,
              coefficients   = coefficients,
              tstat          = tstat,
              p.value        = p.value,
              featurelist    = featurelist,
              ...))
}
# --------------------------------------------------------

##' @title topTableImpute
##' @name topTableImpute
##' @description This function is a modification of the topTable function in the limma package to work with the imputated dataset
##' @param fit is a CombinedMArrayLM object
##' @param sort.by statistic to sort by; choices are p.value ("p" or "P"), t-value ("t" or "T"), "logFC"
##' @param coef coefficient of interest from the fitted model
##' @param adjust.method method to use for adjusting p-values for multiple testing (see help page of p.adjust() for list of options)
##' @param lfc minimum absolute log2-fold-change required
##' @param p.value cutoff value for adjust p-values. Only proteins with lower p-values are listed
##' @param featurelist data.frame or character vector containing protein information, e.g. protein names
##' @param number maximum number of proteins to list
##' @return top result table
##' @export
setMethod("topTableImpute", "CombinedMArrayLM", function(fit,
                                                         sort.by = "p",
                                                         coef    = 1,
                                                         adjust.method = "BH",
                                                         lfc     = 0,
                                                         p.value = 1,
                                                         featurelist = NULL,
                                                         number  = 10) {

   if (!is.null(featurelist) && is.null(dim(featurelist))) {
      featurelist <- data.frame(ID = featurelist, stringsAsFactors = FALSE)
   }

   M        <- as.matrix(fit@coefficients)[,coef]
   tstat    <- as.matrix(fit@tstat)[,coef]
   P.value  <- as.matrix(fit@p.value)[,coef]
   rownum   <- 1:length(M)

   sort.by <- match.arg(sort.by, c("logFC", "M", "P", "p", "T", "t"))

   if(sort.by == "M")
      sort.by <- "logFC"
   if(sort.by == "T")
      sort.by <- "t"
   if(sort.by == "p")
      sort.by <- "P"

   # Adjust the raw p-values for multiple testing using method specified in adjust.method (e.g. BH, FDR, etc.)
   adj.P.value <- p.adjust(P.value, method=adjust.method)

   # If a specific p.value threshold or lfc threshold is specified, find the significantly DE proteins according
   # to those thresholds.
   if(p.value < 1 | lfc > 0) {
      sig <- (adj.P.value < p.value) & (abs(M) > lfc)
      if(any(is.na(sig)))
         # Check for NA values
         sig[is.na(sig)] <- FALSE
      if(all(!sig))
         # If no significantly DE proteins, return an empty dataframe
         return(data.frame())

      # Extract significant proteins (name, logFC, tstat, p.val, p.val.adj)
      featurelist <- featurelist[sig, , drop=F]
      M           <- M[sig]
      tstat       <- tstat[sig]
      P.value     <- P.value[sig]
      adj.P.value <- adj.P.value[sig]
      rownum      <- rownum[sig]
   }

   # Sort the proteins according to the value specified in sort.by (logFC, p.value or t.value)
   ord <- switch(sort.by,
                 logFC = order(abs(M), decreasing=T),
                 P = order(P.value, decreasing=F),
                 t = order(abs(tstat), decreasing=T))

   # Check that the number of top proteins requested is not more than the number of proteins in the dataset.
   if(length(M) < number)
      number <- length(M)

   top <- ord[1:number]

   if(is.null(featurelist))
      top_results_table <- data.frame(logFC = M[top])
   else
      top_results_table <- data.frame(featurelist[top,,drop=F], logFC=M[top], stringsAsFactors=F)

   # Create the sorted dataframe for outputting
   top_results_table <- data.frame(top_results_table,
                                   t = tstat[top],
                                   P.value = P.value[top],
                                   adj.P.value = adj.P.value[top])

   rownames(top_results_table) <- as.character(rownum)[top]
   return(top_results_table)
})



##' @title get significant features
##' @name getSignificantFeatures
##' @description This function is a modification of the topTable function in the limma package to work with the imputated dataset
##' @param fit is a CombinedMArrayLM object
##' @param sort.by statistic to sort by; choices are p.value ("p" or "P"), t-value ("t" or "T"), "logFC"
##' @param coef coefficient of interest from the fitted model
##' @param adjust.method method to use for adjusting p-values for multiple testing (see help page of p.adjust() for list of options)
##' @param lfc minimum absolute log2-fold-change required
##' @param p.value cutoff value for adjust p-values. Only proteins with lower p-values are listed
##' @param featurelist data.frame or character vector containing protein information, e.g. protein names
##' @param number maximum number of proteins to list
##' @param logFCcutoff function which defines the log fold change cutoff. NA if no cutoff is applied
##' @param adj.P.value adjusted p-value cutoff. NA is none
##' @return top result table
##' @export
setMethod("getSignificantFeatures", "CombinedMArrayLM", function(fit,
                                                                 p.value      = 1,
                                                                 adj.P.value  = 0.05,
                                                                 logFCcutoff  = NA,
                                                                 ...) {
   tt <- topTableImpute(fit, number=nrow(fit@coefficients), p.value=p.value, ...)

   if(is.function(logFCcutoff)) {
         tt <- tt[logFCcutoff(tt[,"logFC"]), ]
   }

   if(!is.na(adj.P.value)) {
      if(is.numeric(adj.P.value)) {
         tt <- tt[tt[,"adj.P.value"] < adj.P.value,]
      } else {
         stop("adj.P.value has to be numeric")
      }
   }

   return(tt)
})


##' @title writeSignificantFeatures
##' @name writeSignificantFeatures
##' @description writes signifiant features
##' @param fit is a CombinedMArrayLM object
##' @param file file name
##' @export
setMethod("writeSignificantFeatures", c(fit="CombinedMArrayLM", file="character"), function(fit,
                                                                                            file,
                                                                                            ...) {
   tt <- getSignificantFeatures(fit,  ...)
   write.table(tt, file=file, sep="\t", quote=F, row.names=F)
})

##' @title getFeatureFilter
##' @name getFeatureFilter
##' @description get feature filters
##' @param fit is a CombinedMArrayLM object
##' @param data data
##' @param p.value p.value index
##' @param adjust multiple hypothesis testing method
##' @param onlyPositive get features with positive logFC
##' @param onlyNegative get features with negative logFC
##' @param mode and / or logic
##' @return filter
##' @export
setMethod("getFeatureFilter", c(fit="CombinedMArrayLM", data="MImputedExpressionSets"), function(fit,
                                                                                                 data,
                                                                                                 p.value=0.05,
                                                                                                 adjust="BH",
                                                                                                 onlyPositive=T,
                                                                                                 onlyNegative=F,
                                                                                                 mode="or",
                                                                                                 ...) {
   ncoef <- ncol(fit@coefficients)

   tt <- lapply(1:ncoef, function(x) {
      return(sort(as.numeric(rownames(getSignificantFeatures(fit, coef=x, adjust=adjust, p.value=p.value, onlyPositive=onlyPositive, onlyNegative=onlyNegative)))))#, ...))
   })

   if(mode == "and") {
      differ <- tt[[1]]
      if(length(tt) > 1) {
         for(i in 2:length(tt)) {
            differ <- intersect(differ, tt[[i]])
         }
      }
   } else if (mode == "or") {
      differ <- unique(sort(unlist(tt)))
   } else {
      stop("incorrect mode")
   }

   return(differ)
})

##' @title filterFeatures
##' @name filterFeatures
##' @description filter featues
##' @param fit is a CombinedMArrayLM object
##' @param data data
##' @return filtered data
##' @export
setMethod("filterFeatures", c(fit="CombinedMArrayLM", data="MImputedExpressionSets"), function(fit,
                                                                                               data,
                                                                                               ...) {
   return(data[getFeatureFilter(fit, data, ...),])
})
