# --------------------------------------------------------
# Class MImputedMArrayLM
# Authors: Thomas Schwarzl <schwarzl@embl.de> with help from Elisa D'Arcangelo
# MArrayLM of multiple imputed data
# --------------------------------------------------------
# Constructor:
##' @rdname MImputedExpressionSets
##' @export
MImputedMArrayLM <- function(data,
                             ...) {
   return(new(Class  = "MImputedMArrayLM",
              data   = data,
              ...))
}
# --------------------------------------------------------


##' @title getData
##' @name getData
##' @description get the data
##' @return data
##' @export
setMethod("getData", "MImputedMArrayLM", function(object) {
   return(object@data)
})

##' @title get the number of imputations
##' @name getNumberImputations
##' @description returns the number of imputations.
##' @export
setMethod("getNumberImputations", "MImputedMArrayLM", function(object) {
   return(length(object@data))
})



# setGeneric("data", function(object) standardGeneric("data"))

# setGeneric("data<-", function(object, value) standardGeneric("data<-"))


# setMethod(data, "MImputedMArrayLM", function(object) {
#    slot(object, "data")
# })


# setReplaceMethod("data", "MImputedMArrayLM", function(object, value){
#    slot(object, "data") <- value
#    validObject(object)
#    return(object)
# })


##' @name contrastFit
##' @title contrast fit
##' @description TODO
##' @param contrasts contrasts
##' @importMethodsFrom limma eBayes makeContrasts
##' @export
setMethod("contrastFit", "MImputedMArrayLM", function(fit, contrasts) {
   contrastMatrix <- makeContrasts(contrasts=contrasts, levels=fit@data[[1]]$design)

   #fit each probeset to model
   fit@data <- lapply(fit@data, function(x) {
      return(eBayes(contrasts.fit(x, contrastMatrix)))
   })

   return(fit)
})

##' @name contrastLimit
##' @title contrastLimit
##' @description TODO
##' @param contrasts contrasts
##' @importMethodsFrom limma eBayes
##' @export
setMethod("contrastLimit", "MImputedMArrayLM", function(fit, contrasts) {
   # for each contrast, get the detection limits
   fit@data <- lapply(fit@data, function(x) {
      return(eBayes(contrasts.fit(x, contrastMatrix)))
   })

   return(fit)
})

##' @name vsFit
##' @title vsFit
##' @description TODO
##' @param design design matrix
##' @param vs vs TODO
##' @return return RODO
##' @export
setMethod("vsFit", c(fit    = "MImputedMArrayLM",
                     design = "matrix",
                     vs     = "character"), function(fit, design, vs) {
   if (!vs %in% colnames(design)) {
      stop("'vs' is not a colname of design")
   }
   return(vsFit(fit, design, which(colnames(design)==vs)))
})

##' @rdname vsFit
##' @export
setMethod("vsFit", c(fit    = "MImputedMArrayLM",
                     design = "matrix",
                     vs     = "numeric"), function(fit, design, vs) {
   #create contrast
   selected.group <- colnames(design)[vs]
   other.groups   <- colnames(design)[-vs]
   contrast       <- paste(selected.group, "=", selected.group, " - (", paste(other.groups, collapse="+"), ")/", length(other.groups), sep="")

   return(limmaContrastFit(data.list, design, contrast))
})

##' @name checkMissingness
##' @title check missingness
##' @description checks the missingness
##' @return TODO
##' @export
setMethod("checkMissingness", "MImputedMArrayLM", function(data) {
   checkMissingness(eset(data, 1))
})

##' @name combineFits
##' @title combineFits
##' @param fit fit
##' @description combines fits
##' @return TODO
##' @export
setMethod("combineFits", "MImputedMArrayLM", function(fit) {
   efit.list <- fit@data
   featurelist <- efit.list[[1]]$genes
   m <- length(efit.list) # Number of imputed datasets.
   nproteins <- nrow(efit.list[[1]]$coeff)  # Number of proteins.
   coef.all <- lapply(efit.list, coef)  # Extract the model coefficients for each imputed dataset.

   # Set up matrices for storing results.
   t.val <- p.value <- p.value.adj <- coefficients <- matrix(NA, nrow=nproteins, ncol=ncol(coef.all[[1]]))

   # Apply to each gene separately.
   for(j in 1:nproteins){

      # Need to get the average of the coefficients and variance-covariance of the coefficients
      # for the m imputed datasets.

      # Extract the m variance-covariance matrices.
      variances <- lapply(efit.list, function(x) x$cov.coefficients*x$s2.post[j])

      # Extract the m vectors of model coefficients.
      coef <- lapply(coef.all, function(x) x[j,])

      # Average of the variance-covariance matrices.
      sum.variance <- Reduce('+', variances)
      vbar <- sum.variance/m

      # Average of the coefficient vectors.
      sum.coefs <- Reduce('+', coef)
      cbar <- sum.coefs/m

      # Calculate the between imputation variance-covariance.
      evar <- var(do.call("rbind", coef))

      # Calculate the total variance for the imputed data.
      variance.tot <- vbar+evar*(1/m+1)

      # Calculate the df for the imputed data.
      r <- (1 + 1/m)*evar/vbar
      df <- (m-1)*(1+1/r)^2
      if(is.matrix(df))
         df <- diag(df)
      if(is.matrix(r))
         r <- diag(r)

      # Calculate the t.values = coef/SE(coef).
      temp <- data.frame(coefficients=cbar, SE=sqrt(diag(variance.tot)))
      t.val[j,] <- temp[,1]/temp[,2]
      coefficients[j,] <- cbar
      # Calculate the corresponding p.values using Rubin's df for imputed
      # data.
      p.value[j,] <- 2*pt(abs(t.val[j,]), df, lower.tail=FALSE)
   }

   ids <- rownames(efit.list[[1]])
   out <- CombinedMArrayLM(ids = ids,
                           coefficients=coefficients,
                           tstat=t.val,
                           p.value=p.value,
                           featurelist=featurelist)
   return(out)
})


