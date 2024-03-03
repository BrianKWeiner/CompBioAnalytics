#' Clean Up Metadata in ExpressionSet Object Based on Condition
#'
#' Conditionally updates metadata columns within an ExpressionSet object. For rows where
#' "characteristics_ch1.4" is missing, it shifts data from "characteristics_ch1.3" to "characteristics_ch1.4",
#' and "characteristics_ch1.2" to "characteristics_ch1.3", then sets "characteristics_ch1.2" to "treatment: NA".
#' Rows with data in "characteristics_ch1.4" are left unaffected.
#'
#' @param eset An `ExpressionSet` object containing metadata that needs conditional cleaning.
#' @return An updated `ExpressionSet` object with conditionally cleaned metadata.
#' @importFrom Biobase phenoData pData
#' @examples
#' # Assuming `eset` is your ExpressionSet object:
#' # esetCleaned <- cleanUpMetadataConditionally(eset)
#' @export
cleanUpMetadataForGSE53552 <- function(eset) {
  requireNamespace("Biobase", quietly = TRUE)

  # Validate input
  if (!inherits(eset, "ExpressionSet")) {
    stop("Input must be an ExpressionSet object.")
  }

  # Access the phenotype data
  phenoData <- pData(eset)

  # Ensure required columns exist
  requiredColumns <- c("characteristics_ch1.2", "characteristics_ch1.3", "characteristics_ch1.4")
  if (!all(requiredColumns %in% colnames(phenoData))) {
    stop("Not all required metadata columns are present in the ExpressionSet object.")
  }

  # Identify rows where "characteristics_ch1.4" is NA or empty
  rowsToUpdate <- which(is.na(phenoData$characteristics_ch1.4) | phenoData$characteristics_ch1.4 == "")

  # Conditional update based on rows identified
  if (length(rowsToUpdate) > 0) {
    phenoData$characteristics_ch1.4[rowsToUpdate] <- phenoData$characteristics_ch1.3[rowsToUpdate]
    phenoData$characteristics_ch1.3[rowsToUpdate] <- phenoData$characteristics_ch1.2[rowsToUpdate]
    phenoData$characteristics_ch1.2[rowsToUpdate] <- "treatment: NA"
  }

  # Update the phenoData in the eset object
  pData(eset) <- phenoData

  return(eset)
}
