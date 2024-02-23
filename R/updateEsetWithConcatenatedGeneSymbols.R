#' Update ExpressionSet with Concatenated Gene Symbols
#'
#' Replaces Affymetrix probe IDs in an `ExpressionSet` object with a concatenation of the original probe ID and
#' corresponding HGNC gene symbol, separated by an underscore. This helps preserve the uniqueness of probe IDs
#' while providing gene symbol information.
#'
#' @param eset An `ExpressionSet` object containing Affymetrix probe IDs.
#' @return An updated `ExpressionSet` object with feature names replaced by concatenated probe ID and
#'         gene symbol strings. If a gene symbol is not available for a probe ID, the original probe ID is retained.
#' @importFrom Biobase featureNames
#' @examples
#' # Assuming `eset` is your ExpressionSet object and `mapProbesToSymbols` function is defined and available:
#' # esetUpdated <- updateEsetWithConcatenatedGeneSymbols(eset)
#' @export
updateEsetWithConcatenatedGeneSymbols <- function(eset) {
  # Assuming mapProbesToSymbols function is already defined and available
  geneSymbols <- mapProbesToSymbols(eset)

  # Retrieve the original probe IDs from the eset
  originalProbeIds <- featureNames(eset)

  # Initialize a vector to hold the new feature names, preserving the order
  newFeatureNames <- vector("character", length = length(originalProbeIds))

  # Update feature names with concatenated probe ID and gene symbol
  for (i in seq_along(originalProbeIds)) {
    probeId <- originalProbeIds[i]
    symbol <- ifelse(probeId %in% names(geneSymbols), geneSymbols[probeId], "")
    newFeatureNames[i] <- ifelse(nchar(symbol) > 0, paste(probeId, symbol, sep = "_"), probeId)
  }

  # Update the featureNames of the eset
  featureNames(eset) <- newFeatureNames

  return(eset)
}
