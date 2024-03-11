#' Match Metadata to Principal Components
#'
#' This function matches the metadata (covariates) in an ExpressionSet object to the principal components of expression data. It uses the degCovariates function from the DEGreport R library.
#'
#' @param eset An ExpressionSet object containing expression data and metadata.
#' @param metadataCols A character vector specifying the column names in the metadata to be analyzed.
#' @return A dataframe listing the principal component ID, the percent of variation explained by the PC, the column name of the metadata, and the confidence that this covariate is explained by the principal component.
#' @importFrom DEGreport degCovariates
#' @importFrom Biobase exprs pData
#' @importFrom stats prcomp
#' @examples
#' # Assuming 'eset' is your ExpressionSet:
#' result <- matchMetadataToPCs(eset, metadataCols = c("covariate1", "covariate2"))
#' @export
matchMetadataToPCs <- function(eset, metadataCols) {
  if (!requireNamespace("DEGreport", quietly = TRUE)) {
    stop("DEGreport package is not installed. Please install it to use this function.")
  }

  # For testing
  #GSE53552_rma_genes_with_metadata@phenoData@data$batch <- batch_info$Batch
  #eset <- GSE53552_rma_genes_with_metadata
  #exprData <- Biobase::exprs(GSE53552_rma)
  #metadataCols = c("source_name_ch1 organism_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2", "characteristics_ch1.3", "characteristics_ch1.4", "batch", "ARM")


  # Extract expression data
  exprData <- Biobase::exprs(eset)

  pcaResult <- stats::prcomp(t(exprData), scale. = TRUE)
  pcaResultSummary <- summary(pcaResult)
  percentVarianceExplained <- pcaResultSummary$importance[2,] * 100

    # Extract metadata (covariates) specified by metadataCols
  phenoData <- as.data.frame(Biobase::pData(eset))
  metadataColsIndices <- which(colnames(phenoData) %in% metadataCols)
  phenoData <- phenoData[, metadataColsIndices]

  # Ensure that row names of phenoData match the column names of exprData to align with PCA
  if (!all(rownames(phenoData) == colnames(exprData))) {
    stop("Mismatch between sample names in expression data and phenotype data.")
  }

  # Use degCovariates to find association between PCs and covariates
  degResults <- DEGreport::degCovariates(exprData, phenoData, minPC = length(phenoData))

  # Format the results to include PC ID, percent variation explained, metadata column name, and confidence
  formattedResults <- data.frame(
    PCID = degResults$significants$PC,
    PercentVariationExplained = percentVarianceExplained[1:length(degResults$significants$PC)],
    MetadataCol = degResults$significants$term,
    Confidence = degResults$significants$p.value
  )

  # Order results by PCID and remove those with no cofidence (i.e. NAs)
  formattedResults <- formattedResults[order(formattedResults$PercentVariationExplained, decreasing = T), ]
  formattedResults <- formattedResults[!is.na(formattedResults$Confidence),]

  return(formattedResults)
}
