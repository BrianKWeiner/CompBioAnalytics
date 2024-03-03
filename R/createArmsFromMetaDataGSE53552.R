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
#' # esetUpdated <- createArmsFromMetaDataGSE53552(eset)
#' @export
createArmsFromMetaDataGSE53552 <- function(eset) {
  requireNamespace("Biobase", quietly = TRUE)

  # Validate input
  if (!inherits(eset, "ExpressionSet")) {
    stop("Input must be an ExpressionSet object.")
  }

  # Access the phenotype data
  phenoData <- pData(eset)

  # Ensure required columns exist
  requiredColumns <- c("characteristics_ch1.1", "characteristics_ch1.2", "characteristics_ch1.3")
  if (!all(requiredColumns %in% colnames(phenoData))) {
    stop("Not all required metadata columns are present in the ExpressionSet object.")
  }

  # Strip text before and including the colon, then combine the columns
  cleanedCols <- lapply(phenoData[, requiredColumns], function(column) {
    sapply(column, function(value) {
      parts <- strsplit(as.character(value), ":\\s*")[[1]]
      if (length(parts) > 1) tail(parts, n = 1) else value
    })
  })

  # Combine cleaned columns into a new ARM column
  phenoData$ARM <- apply(do.call(cbind, cleanedCols), 1, paste, collapse = "_")

  # Now replace all the ARM data to get rid of spaces and hyphens so this does
  # not interfear with limma an contrasts.


  phenoData$ARM <- gsub("[-().; ]+", "_", phenoData$ARM)
  phenoData$ARM <- gsub("\\s+", "_", phenoData$ARM)


  # Update the phenoData in the eset object
  pData(eset) <- phenoData

  return(eset)
}
