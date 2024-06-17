#' Create ARM Column from Metadata in ExpressionSet
#'
#' Takes an `ExpressionSet` object and creates a new metadata column `ARM` by combining and cleaning
#' information from `characteristics_ch1.1`, `characteristics_ch1.2`, and `characteristics_ch1.3`.
#'
#' @param eset An `ExpressionSet` object containing metadata to be processed.
#' @return An updated `ExpressionSet` object with a new `ARM` column in the metadata.
#' @importFrom Biobase phenoData pData
#' @examples
#' # Assuming `eset` is your ExpressionSet object:
#' # esetUpdated <- createArmsFromMetaDataGSE145513(eset)
#' @export
createArmsFromMetaDataGSE145513 <- function(eset) {
  requireNamespace("Biobase", quietly = TRUE)

  # Validate input
  if (!inherits(eset, "ExpressionSet")) {
    stop("Input must be an ExpressionSet object.")
  }

  # Access the phenotype data
  phenoData <- Biobase::pData(eset)

  phenoData$characteristics_ch14 <- paste0(phenoData$`characteristics_ch1.2`, "_",
                                           phenoData$`characteristics_ch1.3`)

  # Ensure required columns exist
  requiredColumns <- c("characteristics_ch14")
  if (!all(requiredColumns %in% colnames(phenoData))) {
    stop("Not all required metadata columns are present in the ExpressionSet object.")
  }

  phenoData$ARM <- phenoData$characteristics_ch14

  # Cleans up the special characters that complicate things letter when it comes to graphing
  phenoData$ARM <- gsub("agent: ", "", phenoData$ARM)
  phenoData$ARM <- gsub("\\+", "_pos", phenoData$ARM)
  phenoData$ARM <- gsub("-", "_neg", phenoData$ARM)
  phenoData$ARM <- gsub(" ", "_", phenoData$ARM)
  phenoData$ARM <- gsub("\\(", "_", phenoData$ARM)
  phenoData$ARM <- gsub("\\)", "", phenoData$ARM)

  # Update the phenoData in the eset object
  pData(eset) <- phenoData

  return(eset)
}
