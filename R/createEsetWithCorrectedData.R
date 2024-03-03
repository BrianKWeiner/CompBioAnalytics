#' Combine Corrected Data with Existing ExpressionSet
#'
#' This function creates a new `ExpressionSet` object with batch-corrected expression data,
#' while retaining and updating the metadata (experimentData, phenoData, featureData, annotation,
#' and protocolData) from an existing `ExpressionSet` object.
#'
#' @param correctedData A matrix of batch-corrected expression data.
#' @param affyBatch An `AffyBatch` object from which initial metadata is derived.
#' @param existingEset An existing `ExpressionSet` object to use for updating metadata in the new `ExpressionSet`.
#' @return A new `ExpressionSet` object containing the corrected expression data along with updated metadata.
#' @importFrom Biobase ExpressionSet experimentData phenoData featureData annotation
#' @examples
#' # Assuming `correctedData` is your matrix of batch-corrected expression data,
#' # `affyBatch` is your original AffyBatch object,
#' # and `existingEset` is an existing ExpressionSet object:
#' # newEset <- createEsetWithCorrectedData(correctedData, affyBatch, existingEset)
#' @export
createEsetWithCorrectedData <- function(correctedData, affyBatch, existingEset) {
  requireNamespace("Biobase", quietly = TRUE)

  # Validate inputs
  if (!inherits(affyBatch, "AffyBatch")) {
    stop("affyBatch must be an AffyBatch object")
  }
  if (!is.matrix(correctedData)) {
    stop("correctedData must be a matrix")
  }
  if (!inherits(existingEset, "ExpressionSet")) {
    stop("existingEset must be an ExpressionSet object")
  }

  # Ensure the sample names in correctedData match those in existingEset
  if (!all(colnames(correctedData) %in% sampleNames(existingEset))) {
    stop("Sample names in correctedData must match those in existingEset")
  }

  # Create a new ExpressionSet object with the corrected data
  newEset <- new("ExpressionSet", exprs = correctedData)

  # Update metadata from existingEset
  experimentData(newEset) <- experimentData(existingEset)
  phenoData(newEset) <- phenoData(existingEset)
  featureData(newEset) <- featureData(existingEset)
  annotation(newEset) <- annotation(existingEset)
  protocolData(newEset) <- protocolData(existingEset)

  # Align sample names between the newEset and the existingEset
  sampleNames(newEset) <- sampleNames(existingEset)[match(colnames(correctedData), sampleNames(existingEset))]

  return(newEset)
}
