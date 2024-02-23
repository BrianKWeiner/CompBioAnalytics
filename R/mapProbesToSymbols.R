#' Map Affymetrix Probe IDs to HGNC Gene Symbols
#'
#' This function maps Affymetrix probe IDs to their corresponding HGNC gene symbols using
#' the specified Bioconductor annotation package. It handles multiple mappings by choosing
#' one gene symbol per probe ID, aiming for consistency.
#'
#' @param eset An `ExpressionSet` object containing Affymetrix probe IDs.
#' @return A named vector with Affymetrix probe IDs as names and corresponding HGNC gene symbols as values.
#'         If a probe ID maps to multiple gene symbols, one symbol is arbitrarily selected.
#' @importFrom AnnotationDbi select
#' @importFrom hgu133plus2.db hgu133plus2SYMBOL
#' @examples
#' # Assuming `eset` is your ExpressionSet object:
#' # geneSymbols <- mapProbesToSymbols(eset)
#' @export
mapProbesToSymbols <- function(eset) {
  requireNamespace("AnnotationDbi", quietly = TRUE)
  requireNamespace("hgu133plus2.db", quietly = TRUE)

  # Extract probe IDs from the ExpressionSet object
  probeIds <- featureNames(eset)

  # Map probe IDs to gene symbols
  mappings <- AnnotationDbi::select(hgu133plus2.db,
                                    keys = probeIds,
                                    columns = "SYMBOL",
                                    keytype = "PROBEID")

  # Create a named vector: probe IDs as names, gene symbols as values
  symbols <- setNames(mappings$SYMBOL, mappings$PROBEID)

  return(symbols)
}
