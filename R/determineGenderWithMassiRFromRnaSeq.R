#' Determine Gender Using the massiR Package
#'
#' This function predicts the gender of samples in an ExpressionSet object using
#' Y chromosome probes. It accounts for different platforms by mapping annotations
#' and assumes RNA-seq if the platform is not recognized.
#'
#' @param eset An ExpressionSet object containing expression data and phenotype data.
#' @return An updated ExpressionSet object with an added column in phenoData named
#'         PredictedGender, which contains the gender prediction for each sample.
#' @importFrom massiR massi_y massi_y_plot massi_cluster massi_cluster_plot
#' @importFrom biomaRt useMart getBM
#' @importFrom Biobase ExpressionSet exprs pData
#' @examples
#' # Assuming 'eset' is an ExpressionSet object with expression and phenotype data:
#' esetWithGender <- determineGenderWithMassiR(eset)
#' @export
determineGenderWithMassiRFromRnaSeq <- function(eset) {
  if (!requireNamespace("massiR", quietly = TRUE)) {
    stop("Package 'massiR' is required but not installed.")
  }

  # Find the location of the XIST gene (ENSG00000229807)
  xistGene <- which(rownames(eset@assayData$exprs) == "ENSG00000229807")
  mean_Xist <- mean(eset@assayData$exprs[xistGene,])

  #Some samples are all female and this bypasses the whole thing
  if(mean_Xist < 10) {
    # Add PredictedGender to the phenotype data of the ExpressionSet
    phenoData <- pData(eset)
    listPredictedGender <- rep("female", length(colnames(Biobase::exprs(eset))))
    phenoData$PredictedGender <- listPredictedGender

    # Update the ExpressionSet with the new phenotype data
    pData(eset) <- phenoData

    return(eset)
  }

  # Load y chromosome probes data
  data(y.probes)

  # Determine Y chromosome test probes based on eset annotation
  determineYChromosomeTestProbes <- function(eset, y.probes) {
    platformMappings <- list(
      hgu133plus2 = "affy_hg_u133_plus_2",
      illuminaHumanwg6V1 = "illumina_humanwg_6_v1",
      illuminaHumanwg6V2 = "illumina_humanwg_6_v2",
      illuminaHumanht12 = "illumina_humanht_12",
      affyHugene10StV1 = "affy_hugene_1_0_st_v1"
    )

    platform <- platformMappings[[eset@annotation]]

    if (!is.null(platform) && !is.null(y.probes[[platform]])) {
      yChromosomeTestProbes <- y.probes[[platform]]
    } else {
      print("Assuming RNA-seq or platform not recognized.")
      yChromosomeTestProbes <- determineRnaSeqy.probes()
    }

    return(yChromosomeTestProbes)
  }

  determineRnaSeqy.probes <- function() {
    mart <- biomaRt::useMart('ensembl', dataset = "hsapiens_gene_ensembl")
    geneAttributes <- biomaRt::getBM(attributes = c("ensembl_gene_id", "chromosome_name"),
                                     filters = 'chromosome_name', values = 'Y', mart = mart)
    yChromosomeProbes <- data.frame(row.names = geneAttributes$ensembl_gene_id)
    return(yChromosomeProbes)
  }

  yChromosomeTestProbes <- determineYChromosomeTestProbes(eset, y.probes)

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
