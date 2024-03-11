#' Determine Gender Using massiR Package
#'
#' @param eset An ExpressionSet object containing expression data and phenotype data.
#' @return An ExpressionSet object with an added column in phenoData called PredictedGender.
#' @importFrom massiR massi_y massi_y_plot massi_cluster massi_cluster_plot
#' @examples
#' # Assuming 'eset' is your ExpressionSet object:
#' esetWithGender <- determineGenderWithMassiR(eset)
#' @export
determineGenderWithMassiR <- function(eset) {
  if (!requireNamespace("massiR", quietly = TRUE)) {
    stop("Package 'massiR' is required but not installed.")
  }

  # For testing
  eset <- GSE53552_rma

  data(y.probes)
  # This package currently supports:
  # "illumina_humanwg_6_v1" "illumina_humanwg_6_v2" "illumina_humanwg_6_v1" "illumina_humanht_12"   "affy_hugene_1_0_st_v1" "affy_hg_u133_plus_2"
  if(eset@annotation == "hgu133plus2") {
    yChromosomeTestProbes <- y.probes$affy_hg_u133_plus_2
  }

  # Extract expression data from the ExpressionSet
  exprData <- Biobase::exprs(eset)

  # Use massiR to predict gender
  massiYOutput <- massiR::massi_y(expression_data = eset, y_probes = yChromosomeTestProbes)

  # Plot the data to take a quick look:
  massiR::massi_y_plot(massiYOutput)

  # Select the threshold: The threshold can be determined by quantiles of probe
  # variance (CV): 1=All probes, 2=Upper 75%, 3=Upper 50%, 4=Upper 25%.
  massiSelectOutput <- massiR::massi_select(expression_data = eset, y_probes = yChromosomeTestProbes, threshold=3)

  # Finally predict the gender based on the selected probes:
  resultsOfGenderPrediction <- massiR::massi_cluster(massiSelectOutput)

  massiR::massi_cluster_plot(massiSelectOutput, resultsOfGenderPrediction)


  # Add PredictedGender to the phenotype data of the ExpressionSet
  phenoData <- pData(eset)
  phenoData$PredictedGender <- resultsOfGenderPrediction$massi.results$sex

  # Recommended code from vignette:
  # Get the sex for each sample from the massi_cluster results
  #esetSampleResults <- data.frame(resultsOfGenderPrediction[[2]])
  #sexData <- data.frame(esetSampleResults[c("ID", "sex")])
  # Extract the order of samples in the ExpressionSet and match with results
  #esetNames <- colnames(exprs(eset))
  # match the sample order in massiR results to the same as the ExpressionSet object
  #exData <- sexData[match(esetNames, sexData$ID),]
  # create an annotatedDataFrame to add to ExpressionSet
  #newPhenoData <- new("AnnotatedDataFrame", data = sexData)
  # add the annotatedDataFrame to the Expressionset as phenoData
  #phenoData(eset) <- newPhenoData

  # Update the ExpressionSet with the new phenotype data
  pData(eset) <- phenoData

  return(eset)
}
