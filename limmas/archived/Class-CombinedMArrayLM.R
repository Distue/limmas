# --------------------------------------------------------
# Class MImputedExpressionSets
# holds ExpressionSets of multiple imputed data
# --------------------------------------------------------
setClass("CombinedMArrayLM",
         representation=representation(
            ids            = "character",
            coefficients   = "matrix",
            tstat          = "matrix",
            p.value        = "matrix",
            featurelist    = "data.frame")
         )
         


# assertions
setValidity("CombinedMArrayLM", function(object) {
   msg <- NULL
   
   if (is.null(msg)) TRUE else msg
})            


# This function is a modification of the topTable function in the limma package.
# fit is a fitted object coming from combine.MI.fits()
# sort.by = statistic to sort by; choices are p.value ("p" or "P"), t-value ("t" or "T"), "logFC"
# coef = coefficient of interest from the fitted model
# adjust.method = method to use for adjusting p-values for multiple testing 
# (see help page of p.adjust() for list of options)
# lfc = minimum absolute log2-fold-change required
# p.value = cutoff value for adjust p-values. Only proteins with lower p-values are listed
# featurelist = dataframe or character vector containing protein information, e.g. protein names
# number = maximum number of proteins to list

setGeneric("topTableImpute", function(fit, ...) standardGeneric("topTableImpute"))
setMethod("topTableImpute", "CombinedMArrayLM", function(fit, sort.by="p", coef=1, adjust.method="BH", lfc=0, p.value=1, featurelist=NULL, number=10) {

   if (!is.null(featurelist) && is.null(dim(featurelist))) {
      featurelist <- data.frame(ID = featurelist, stringsAsFactors = FALSE)
   }
   
   M <- as.matrix(fit@coefficients)[,coef]
   tstat <- as.matrix(fit@tstat)[,coef]
   P.value <- as.matrix(fit@p.value)[,coef]
   rownum <- 1:length(M)
   
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
   if(p.value < 1 | lfc > 0){
      sig <- (adj.P.value < p.value) & (abs(M) > lfc)
      if(any(is.na(sig))) 
         sig[is.na(sig)] <- FALSE     # Check for NA values
      if(all(!sig))  
         return(data.frame())            # If no significantly DE proteins, return an empty dataframe
      
      featurelist <- featurelist[sig, , drop=F]       # Extract significant proteins (name, logFC, tstat, p.val, p.val.adj)
      M <- M[sig]
      tstat <- tstat[sig]
      P.value <- P.value[sig]
      adj.P.value <- adj.P.value[sig]
      rownum <- rownum[sig]
   }
   
   # Sort the proteins according to the value specified in sort.by (logFC, p.value or t.value)                     
   ord <- switch(sort.by, logFC = order(abs(M), decreasing=T), P = order(P.value, decreasing=F), 
                 t = order(abs(tstat), decreasing=T))
   
   # Check that the number of top proteins requested is not more than the number of proteins in the dataset.                     
   if(length(M) < number) number <- length(M)
   
   top <- ord[1:number]
   if(is.null(featurelist))
      top_results_table <- data.frame(logFC = M[top])
   else
      top_results_table <- data.frame(featurelist[top,,drop=F], logFC=M[top], stringsAsFactors=F)
   
   # Create the sorted dataframe for outputting
   top_results_table <- data.frame(top_results_table, t = tstat[top], P.value = P.value[top], adj.P.value = adj.P.value[top])
   rownames(top_results_table) <- as.character(rownum)[top]
   return(top_results_table)
})


setGeneric("getSignificantGenes", function(fit, ...) standardGeneric("getSignificantGenes"))
setMethod("getSignificantGenes", "CombinedMArrayLM", function(fit, p.value=0.05, ...) {
   tt <- topTableImpute(fit, number=nrow(fit@coefficients), p.value=p.value, ...)
})
   
          
       
             