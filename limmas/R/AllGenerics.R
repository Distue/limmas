##' topTableImpute method generics
##'
##' @docType methods
##' @name topTableImpute
##' @rdname topTableImpute-methods
##' @title topTableImpute method
##' @param ... Additional argument list
##' @return top result table
##' @export
##' @author Thomas Schwarzl
if(!isGeneric("topTableImpute"))
   setGeneric("topTableImpute", function(fit, ...) standardGeneric("topTableImpute"))

##' getSignificantFeatures method generics
##'
##' @docType methods
##' @name getSignificantFeatures
##' @rdname getSignificantFeatures-methods
##' @title getSignificantFeatures method
##' @param ... Additional argument list
##' @return results
##' @export
##' @author Thomas Schwarzl
if(!isGeneric("getSignificantFeatures"))
   setGeneric("getSignificantFeatures", function(fit, ...) standardGeneric("getSignificantFeatures"))


##' writeSignificantFeatures method generics
##'
##' @docType methods
##' @name writeSignificantFeatures
##' @rdname writeSignificantFeatures-methods
##' @title writeSignificantFeatures method
##' @param file file name
##' @param ... Additional argument list
##' @return results
##' @export
##' @author Thomas Schwarzl
if(!isGeneric("writeSignificantFeatures"))
   setGeneric("writeSignificantFeatures", function(fit, file, ...) standardGeneric("writeSignificantFeatures"))


##' getFeatureFilter method generics
##'
##' @docType methods
##' @name getFeatureFilter
##' @rdname getFeatureFilter-methods
##' @title getFeatureFilter method
##' @param data the data
##' @param ... Additional argument list
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("getFeatureFilter"))
   setGeneric("getFeatureFilter", function(fit, data, ...) standardGeneric("getFeatureFilter"))  


##' filterFeatures method generics
##'
##' @docType methods
##' @name filterFeatures
##' @rdname filterFeatures-methods
##' @title filterFeatures method
##' @param data the data
##' @param ... Additional argument list
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("filterFeatures"))
   setGeneric("filterFeatures", function(fit, data, ...) standardGeneric("filterFeatures"))   

##' originalNamesCol method generics
##'
##' @docType methods
##' @name originalNamesCol
##' @rdname originalNamesCol-methods
##' @title originalNamesCol method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("originalNamesCol"))
   setGeneric("originalNamesCol", function(object) standardGeneric("originalNamesCol"))

##' originalNamesCol<- method generics
##'
##' @docType methods
##' @name originalNamesCol<-
##' @rdname originalNamesCol<--methods
##' @title originalNamesCol<- method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("originalNamesCol<-"))
   setGeneric("originalNamesCol<-", function(object, value) standardGeneric("originalNamesCol<-"))

##' sampleNamesCol method generics
##'
##' @docType methods
##' @name sampleNamesCol
##' @rdname sampleNamesCol-methods
##' @title sampleNamesCol method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("sampleNamesCol"))
   setGeneric("sampleNamesCol", function(object) standardGeneric("sampleNamesCol"))

##' sampleNamesCol<- method generics
##'
##' @docType methods
##' @name sampleNamesCol<-
##' @rdname sampleNamesCol<--methods
##' @title sampleNamesCol<- method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("sampleNamesCol<-"))
   setGeneric("sampleNamesCol<-", function(object, value) standardGeneric("sampleNamesCol<-"))

##' getOriginalNames method generics
##'
##' @docType methods
##' @name getOriginalNames
##' @rdname getOriginalNames-methods
##' @title getOriginalNames method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("getOriginalNames"))
   setGeneric("getOriginalNames", function(object) standardGeneric("getOriginalNames"))

##' getSampleNames method generics
##'
##' @docType methods
##' @name getSampleNames
##' @rdname getSampleNames-methods
##' @title getSampleNames method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("getSampleNames"))
   setGeneric("getSampleNames", function(object) standardGeneric("getSampleNames"))

##' getAnnotatedDataFrame method generics
##'
##' @docType methods
##' @name getAnnotatedDataFrame
##' @rdname getAnnotatedDataFrame-methods
##' @title getAnnotatedDataFrame method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("getAnnotatedDataFrame"))
   setGeneric("getAnnotatedDataFrame", function(object) standardGeneric("getAnnotatedDataFrame"))

##' getData method generics
##'
##' @docType methods
##' @name getData
##' @rdname getData-methods
##' @title getData method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("getData"))
   setGeneric("getData", function(object) standardGeneric("getData"))

##' getNumberImputations method generics
##'
##' @docType methods
##' @name getNumberImputations
##' @rdname getNumberImputations-methods
##' @title getNumberImputations method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("getNumberImputations"))
   setGeneric("getNumberImputations", function(object) standardGeneric("getNumberImputations"))

##' contrastFit method generics
##'
##' @docType methods
##' @name contrastFit
##' @rdname contrastFit-methods
##' @title contrastFit method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("contrastFit"))
   setGeneric("contrastFit", function(fit, contrasts) standardGeneric("contrastFit"))

##' contrastLimit method generics
##'
##' @docType methods
##' @name contrastLimit
##' @rdname contrastLimit-methods
##' @title contrastLimit method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("contrastLimit"))
   setGeneric("contrastLimit", function(fit, contrasts) standardGeneric("contrastLimit"))

##' vsFit method generics
##'
##' @docType methods
##' @name vsFit
##' @rdname vsFit-methods
##' @title vsFit method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("vsFit"))
   setGeneric("vsFit", function(fit, design, vs) standardGeneric("vsFit"))

##' combineFits method generics
##'
##' @docType methods
##' @name combineFits
##' @rdname combineFits-methods
##' @title combineFits method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("combineFits"))
   setGeneric("combineFits", function(fit) standardGeneric("combineFits"))

##' getOriginalData method generics
##'
##' @docType methods
##' @name getOriginalData
##' @rdname getOriginalData-methods
##' @title getOriginalData method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("getOriginalData"))
   setGeneric("getOriginalData", function(object) standardGeneric("getOriginalData"))

##' minPresent method generics
##'
##' @docType methods
##' @name minPresent
##' @rdname minPresent-methods
##' @title minPresent method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("minPresent"))
   setGeneric("minPresent", function(object) standardGeneric("minPresent"))

##' minPresent<- method generics
##'
##' @docType methods
##' @name minPresent<-
##' @rdname minPresent<--methods
##' @title minPresent<- method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("minPresent<-"))
   setGeneric("minPresent<-", function(object, value) standardGeneric("minPresent<-"))

##' data method generics
##'
##' @docType methods
##' @name data
##' @rdname data-methods
##' @title data method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("data"))
   setGeneric("data", function(object) standardGeneric("data"))

##' data<- method generics
##'
##' @docType methods
##' @name data<-
##' @rdname data<--methods
##' @title data<- method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("data<-"))
   setGeneric("data<-", function(object, value) standardGeneric("data<-"))

##' groupingCol method generics
##'
##' @docType methods
##' @name groupingCol
##' @rdname groupingCol-methods
##' @title groupingCol method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("groupingCol"))
   setGeneric("groupingCol", function(object) standardGeneric("groupingCol"))

##' groupingCol<- method generics
##'
##' @docType methods
##' @name groupingCol<-
##' @rdname groupingCol<--methods
##' @title groupingCol<- method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("groupingCol<-"))
   setGeneric("groupingCol<-", function(object, value) standardGeneric("groupingCol<-"))

##' getGroupingCol method generics
##'
##' @docType methods
##' @name getGroupingCol
##' @rdname getGroupingCol-methods
##' @title getGroupingCol method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("getGroupingCol"))
   setGeneric("getGroupingCol", function(object) standardGeneric("getGroupingCol"))

##' numberImputations method generics
##'
##' @docType methods
##' @name numberImputations
##' @rdname numberImputations-methods
##' @title numberImputations method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("numberImputations"))
   setGeneric("numberImputations", function(object) standardGeneric("numberImputations"))

##' numberImputations<- method generics
##'
##' @docType methods
##' @name numberImputations<-
##' @rdname numberImputations<--methods
##' @title numberImputations<- method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("numberImputations<-"))
   setGeneric("numberImputations<-", function(object, value) standardGeneric("numberImputations<-"))

##' pData method generics
##'
##' @docType methods
##' @name pData
##' @rdname pData-methods
##' @title pData method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("pData"))
   setGeneric("pData", function(object) standardGeneric("pData"))

##' fData method generics
##'
##' @docType methods
##' @name fData
##' @rdname fData-methods
##' @title fData method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("fData"))
   setGeneric("fData", function(object) standardGeneric("fData"))

##' annotation method generics
##'
##' @docType methods
##' @name annotation
##' @rdname annotation-methods
##' @title annotation method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("annotation"))
   setGeneric("annotation", function(object) standardGeneric("annotation"))

##' eset method generics
##'
##' @docType methods
##' @name eset
##' @rdname eset-methods
##' @title eset method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("eset"))
   setGeneric("eset", function(object, imputation) standardGeneric("eset"))

##' intensities method generics
##'
##' @docType methods
##' @name intensities
##' @rdname intensities-methods
##' @title intensities method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("intensities"))
   setGeneric("intensities", function(object, imputation) standardGeneric("intensities"))

##' limmasFit method generics
##'
##' @docType methods
##' @name limmasFit
##' @rdname limmasFit-methods
##' @title limmasFit method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("limmasFit"))
   setGeneric("limmasFit", function(object, design) standardGeneric("limmasFit"))

##' completeCases method generics
##'
##' @docType methods
##' @name completeCases
##' @rdname completeCases-methods
##' @title completeCases method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("completeCases"))
   setGeneric("completeCases", function(object) standardGeneric("completeCases"))

##' filterRows method generics
##'
##' @docType methods
##' @name filterRows
##' @rdname filterRows-methods
##' @title filterRows method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("filterRows"))
   setGeneric("filterRows", function(object, filter) standardGeneric("filterRows"))

##' filterCols method generics
##'
##' @docType methods
##' @name filterCols
##' @rdname filterCols-methods
##' @title filterCols method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("filterCols"))
   setGeneric("filterCols", function(object, filter) standardGeneric("filterCols"))

##' fillNAsWithValues method generics
##'
##' @docType methods
##' @name fillNAsWithValues
##' @rdname fillNAsWithValues-methods
##' @title fillNAsWithValues method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("fillNAsWithValues"))
   setGeneric("fillNAsWithValues", function(object, value) standardGeneric("fillNAsWithValues"))

##' plotExpression method generics
##'
##' @docType plotExpression methods
##' @name plotExpression
##' @rdname plotExpression-methods
##' @title plotExpression method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("plotExpression"))
   setGeneric("plotExpression", function(object, ID, ...) standardGeneric("plotExpression"))

##' getAverageExpression method generics
##'
##' @docType methods
##' @name getAverageExpression
##' @rdname getAverageExpression-methods
##' @title getAverageExpression method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("getAverageExpression"))
   setGeneric("getAverageExpression", function(object) standardGeneric("getAverageExpression"))

##' calcGroupEstimations method generics
##'
##' @docType methods
##' @name calcGroupEstimations
##' @rdname calcGroupEstimations-methods
##' @title calcGroupEstimations method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("calcGroupEstimations"))
   setGeneric("calcGroupEstimations", function(object) standardGeneric("calcGroupEstimations"))

##' estimateLimits method generics
##'
##' @docType methods
##' @name estimateLimits
##' @rdname estimateLimits-methods
##' @title estimateLimits method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("estimateLimits"))
   setGeneric("estimateLimits", function(object,  design, contrasts) standardGeneric("estimateLimits"))

##' checkMissingness method generics
##'
##' @docType methods
##' @name checkMissingness
##' @rdname checkMissingness-methods
##' @title  method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("checkMissingness"))
   setGeneric("checkMissingness", function(data, ...) standardGeneric("checkMissingness"))

##' filterCols method generics
##'
##' @docType methods
##' @name filterCols
##' @rdname filterCols-methods
##' @title filterCols method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("filterCols"))
   setGeneric("filterCols", function(object, filter) standardGeneric("filterCols"))

##' filterRows method generics
##'
##' @docType methods
##' @name filterRows
##' @rdname filterRows-methods
##' @title filterRows method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("filterRows"))
   setGeneric("filterRows", function(object, filter) standardGeneric("filterRows"))

##' reverseFilter method generics
##'
##' @docType methods
##' @name reverseFilter
##' @rdname reverseFilter-methods
##' @title reverseFilter method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("reverseFilter"))
   setGeneric("reverseFilter", function(data, ...) standardGeneric("reverseFilter"))

##' contaminantFilter method generics
##'
##' @docType methods
##' @name contaminantFilter
##' @rdname contaminantFilter-methods
##' @title contaminantFilter method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("contaminantFilter"))
   setGeneric("contaminantFilter", function(data, ...) standardGeneric("contaminantFilter"))

##' completeCases method generics
##'
##' @docType methods
##' @name completeCases
##' @rdname completeCases-methods
##' @title completeCases method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("completeCases"))
   setGeneric("completeCases", function(object) standardGeneric("completeCases"))

##' normalizeData method generics
##'
##' @docType methods
##' @name normalizeData
##' @rdname normalizeData-methods
##' @title normalizeData method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("normalizeData"))
   setGeneric("normalizeData", function(data, ...) standardGeneric("normalizeData"))

##' transformData method generics
##'
##' @docType methods
##' @name transformData
##' @rdname transformData-methods
##' @title transformData method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("transformData"))
   setGeneric("transformData", function(data, ...) standardGeneric("transformData"))

##' scaleData method generics
##'
##' @docType methods
##' @name scaleData
##' @rdname scaleData-methods
##' @title scaleData method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("scaleData"))
   setGeneric("scaleData", function(data, scalefactor, ...) standardGeneric("scaleData"))

##' getGroupData method generics
##'
##' @docType methods
##' @name getGroupData
##' @rdname getGroupData-methods
##' @title getGroupData method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("getGroupData"))
   setGeneric("getGroupData", function(data, group, ...) standardGeneric("getGroupData"))

##' plotMedianVsNAs method generics
##'
##' @docType methods
##' @name plotMedianVsNAs
##' @rdname plotMedianVsNAs-methods
##' @title plotMedianVsNAs method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("plotMedianVsNAs"))
   setGeneric("plotMedianVsNAs", function(data, group, ...) standardGeneric("plotMedianVsNAs"))

##' plotMedianVsSD method generics
##'
##' @docType methods
##' @name plotMedianVsSD
##' @rdname plotMedianVsSD-methods
##' @title plotMedianVsSD method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("plotMedianVsSD"))
   setGeneric("plotMedianVsSD", function(data, group, ...) standardGeneric("plotMedianVsSD"))

##' plotNAsVsSD method generics
##'
##' @docType methods
##' @name plotNAsVsSD
##' @rdname plotNAsVsSD-methods
##' @title plotNAsVsSD method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("plotNAsVsSD"))
   setGeneric("plotNAsVsSD", function(data, group, ...) standardGeneric("plotNAsVsSD"))

##' fillNAs method generics
##'
##' @docType methods
##' @name fillNAs
##' @rdname fillNAs-methods
##' @title  method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("plotNAdensity"))
   setGeneric("plotNAdensity", function(data, group, ...) standardGeneric("plotNAdensity"))

##' fillNAs method generics
##'
##' @docType methods
##' @name fillNAs
##' @rdname fillNAs-methods
##' @title fillNAs method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("imputeIndependentGroupsWithAmelia"))
   setGeneric("imputeIndependentGroupsWithAmelia", function(data.input, ...) standardGeneric("imputeIndependentGroupsWithAmelia"))

##' fillNAs method generics
##'
##' @docType methods
##' @name fillNAs
##' @rdname fillNAs-methods
##' @title fillNAs method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("fillNAs"))
   setGeneric("fillNAs", function(object, value) standardGeneric("fillNAs"))

##' calculateFeatureCorrelations method generics
##'
##' @docType methods
##' @name calculateFeatureCorrelations
##' @rdname calculateFeatureCorrelations-methods
##' @title  calculateFeatureCorrelations method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("calculateFeatureCorrelations"))
   setGeneric("calculateFeatureCorrelations",  function(object, ...) standardGeneric("calculateFeatureCorrelations"))

##' checkCompleteRows method generics
##'
##' @docType methods
##' @name checkCompleteRows
##' @rdname checkCompleteRows-methods
##' @title  checkCompleteRows method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("checkCompleteRows"))
   setGeneric("checkCompleteRows", function(data, ...) standardGeneric("checkCompleteRows"))


##' getIntensities method generics
##'
##' @docType methods
##' @name getIntensities
##' @rdname getIntensities-methods
##' @title  getIntensities method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("getIntensities"))
   setGeneric("getIntensities", function(data, minIntensity) standardGeneric("getIntensities"))

##' peptideFilter method generics
##'
##' @docType methods
##' @name peptideFilter
##' @rdname peptideFilter-methods
##' @title  peptideFilter method
##' @return results
##' @export
##' @author Thomas Schwarzl 
if(!isGeneric("peptideFilter"))
   setGeneric("peptideFilter", function(data, ...) standardGeneric("peptideFilter"))



