##' @rdname topTableImpute
##' @export
setGeneric("topTableImpute", function(fit, ...) standardGeneric("topTableImpute"))

##' @rdname getSignificantFeatures
##' @export
setGeneric("getSignificantFeatures", function(fit, ...) standardGeneric("getSignificantFeatures"))

##' @rdname writeSignificantFeatures
##' @export
setGeneric("writeSignificantFeatures", function(fit, file, ...) standardGeneric("writeSignificantFeatures"))

##' @rdname getFeatureFilter
##' @export
setGeneric("getFeatureFilter", function(fit, data, ...) standardGeneric("getFeatureFilter"))

##' @rdname filterFeatures
##' @export
setGeneric("filterFeatures", function(fit, data, ...) standardGeneric("filterFeatures"))

##' @rdname originalNamesCol
##' @export
setGeneric("originalNamesCol", function(object) standardGeneric("originalNamesCol"))

##' @rdname originalNamesCol
##' @export
setGeneric("originalNamesCol<-", function(object, value) standardGeneric("originalNamesCol<-"))

##' @rdname sampleNamesCol
##' @export
setGeneric("sampleNamesCol", function(object) standardGeneric("sampleNamesCol"))

##' @rdname sampleNamesCol
##' @export
setGeneric("sampleNamesCol<-", function(object, value) standardGeneric("sampleNamesCol<-"))

##' @rdname getOriginalNames
##' @export
setGeneric("getOriginalNames", function(object) standardGeneric("getOriginalNames"))

##' @rdname getSampleNames
##' @export
setGeneric("getSampleNames", function(object) standardGeneric("getSampleNames"))

##' @rdname getAnnotatedDataFrame
##' @export
setGeneric("getAnnotatedDataFrame", function(object) standardGeneric("getAnnotatedDataFrame"))

##' @rdname getData
##' @export
setGeneric("getData", function(object) standardGeneric("getData"))

##' @rdname getNumberImputations
##' @export
setGeneric("getNumberImputations", function(object) standardGeneric("getNumberImputations"))

##' @rdname contrastFit
##' @export
setGeneric("contrastFit", function(fit, contrasts) standardGeneric("contrastFit"))

##' @rdname contrastLimit
##' @export
   setGeneric("contrastLimit", function(fit, contrasts) standardGeneric("contrastLimit"))

##' @rdname vsFit
##' @export
setGeneric("vsFit", function(fit, design, vs) standardGeneric("vsFit"))

##' @rdname combineFits
##' @export
setGeneric("combineFits", function(fit) standardGeneric("combineFits"))

##' @rdname getOriginalData
##' @export
setGeneric("getOriginalData", function(object) standardGeneric("getOriginalData"))

##' @rdname minPresent
##' @export
setGeneric("minPresent", function(object) standardGeneric("minPresent"))

##' @rdname minPresent
##' @export
setGeneric("minPresent<-", function(object, value) standardGeneric("minPresent<-"))

##' @rdname groupingCol
##' @export
setGeneric("groupingCol", function(object) standardGeneric("groupingCol"))

##' @rdname groupingCol
##' @export
setGeneric("groupingCol<-", function(object, value) standardGeneric("groupingCol<-"))

##' @rdname groupingCol
##' @export
setGeneric("getGroupingCol", function(object) standardGeneric("getGroupingCol"))

##' @rdname numberImputations
##' @export
setGeneric("numberImputations", function(object) standardGeneric("numberImputations"))

##' @rdname numberImputations
##' @export
setGeneric("numberImputations<-", function(object, value) standardGeneric("numberImputations<-"))



##' @rdname eset
##' @export
setGeneric("eset", function(object, imputation) standardGeneric("eset"))

##' @rdname intensities
##' @export
setGeneric("intensities", function(object, imputation) standardGeneric("intensities"))

##' @rdname limmasFit
##' @export
setGeneric("limmasFit", function(object, design) standardGeneric("limmasFit"))

##' @rdname completeCases
##' @export
setGeneric("completeCases", function(object) standardGeneric("completeCases"))

##' @rdname filterRows
##' @export
setGeneric("filterRows", function(object, filter) standardGeneric("filterRows"))

##' @rdname filterCols
##' @export
setGeneric("filterCols", function(object, filter) standardGeneric("filterCols"))

##' @rdname fillNAsWithValues
##' @export
setGeneric("fillNAsWithValues", function(object, value) standardGeneric("fillNAsWithValues"))

##' @rdname plotExpression
##' @export
setGeneric("plotExpression", function(object, ID, ...) standardGeneric("plotExpression"))

##' @rdname getAverageExpression
##' @export
setGeneric("getAverageExpression", function(object) standardGeneric("getAverageExpression"))

##' @rdname calcGroupEstimations
##' @export
setGeneric("calcGroupEstimations", function(object) standardGeneric("calcGroupEstimations"))


##' @rdname estimateLimits
##' @export
setGeneric("estimateLimits", function(object,  design, contrasts) standardGeneric("estimateLimits"))

##' @rdname checkMissingness
##' @export
setGeneric("checkMissingness", function(data, ...) standardGeneric("checkMissingness"))

##' @rdname filterCols
##' @export
setGeneric("filterCols", function(object, filter) standardGeneric("filterCols"))

##' @rdname filterRows
##' @export
setGeneric("filterRows", function(object, filter) standardGeneric("filterRows"))

##' @rdname reverseFilter
##' @export
setGeneric("reverseFilter", function(data, ...) standardGeneric("reverseFilter"))

##' @rdname contaminantFilter
##' @export
setGeneric("contaminantFilter", function(data, ...) standardGeneric("contaminantFilter"))

##' @rdname completeCases
##' @export
setGeneric("completeCases", function(object) standardGeneric("completeCases"))

##' @rdname normalizeData
##' @export
setGeneric("normalizeData", function(data, ...) standardGeneric("normalizeData"))

##' @rdname transformData
##' @export
setGeneric("transformData", function(data, ...) standardGeneric("transformData"))

##' @rdname scaleData
##' @export
setGeneric("scaleData", function(data, scalefactor, ...) standardGeneric("scaleData"))

##' @rdname getGroupData
##' @export
setGeneric("getGroupData", function(data, group, ...) standardGeneric("getGroupData"))

##' @rdname plotMedianVsNAs
##' @export
setGeneric("plotMedianVsNAs", function(data, group, ...) standardGeneric("plotMedianVsNAs"))

##' @rdname plotMedianVsSD
##' @export
setGeneric("plotMedianVsSD", function(data, group, ...) standardGeneric("plotMedianVsSD"))

##' @rdname plotNAsVsSD
##' @export
setGeneric("plotNAsVsSD", function(data, group, ...) standardGeneric("plotNAsVsSD"))

##' @rdname fillNAs
##' @export
setGeneric("plotNAdensity", function(data, group, ...) standardGeneric("plotNAdensity"))

##' @rdname fillNAs
##' @export
setGeneric("imputeIndependentGroupsWithAmelia", function(data.input, ...) standardGeneric("imputeIndependentGroupsWithAmelia"))

##' @rdname fillNAs
##' @export
setGeneric("fillNAs", function(object, value) standardGeneric("fillNAs"))

##' @rdname calculateFeatureCorrelations
##' @export
setGeneric("calculateFeatureCorrelations",  function(object, ...) standardGeneric("calculateFeatureCorrelations"))

##' @rdname checkCompleteRows
##' @export
setGeneric("checkCompleteRows", function(data, ...) standardGeneric("checkCompleteRows"))

##' @rdname getIntensities
##' @export
setGeneric("getIntensities", function(data, minIntensity) standardGeneric("getIntensities"))

##' @rdname peptideFilter
##' @export
setGeneric("peptideFilter", function(data, ...) standardGeneric("peptideFilter"))



