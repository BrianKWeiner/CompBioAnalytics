#' Update Sample Names in AffyBatch Object
#'
#' Truncates sample names in an `AffyBatch` object, typically to remove file extensions and paths,
#' leaving only the unique identifier part of the name (e.g., "GSM1296024").
#'
#' @param affyBatch An `AffyBatch` object as obtained from `affy::ReadAffy()`.
#' @return The `AffyBatch` object with updated sample names.
#' @importFrom affy $.AffyBatch
#' @examples
#' # Assuming `affyBatch` is your AffyBatch object:
#' # affyBatchUpdated <- updateSampleNames(affyBatch)
#' @export
updateSampleNames <- function(affyBatch) {
  requireNamespace("affy", quietly = TRUE)

  # Validate input
  if (!inherits(affyBatch, "AffyBatch")) {
    stop("The input must be an AffyBatch object.")
  }

  # Extract current sample names
  currentNames <- sampleNames(affyBatch)

  # Update sample names by truncating to the first occurrence of "_"
  updatedNames <- sapply(currentNames, function(name) {
    parts <- unlist(strsplit(name, "_", fixed = TRUE))
    return(parts[1])
  })

  # Apply updated sample names to the AffyBatch object
  sampleNames(affyBatch) <- updatedNames

  return(affyBatch)
}
