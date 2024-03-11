#' Estimate Cell Type Proportions Using ImSig
#'
#' @param eset An ExpressionSet object with expression data.
#' @param correlationThreshold The correlation threshold (r) for feature selection in ImSig.
#' @return A dataframe with cell type proportions per sample.
#' @importFrom imsig imsig
#' @examples
#' # Assuming 'eset' is your ExpressionSet object:
#' cellProportions <- estimateCellTypeProportions(eset, correlationThreshold = 0.7)
#' @export
estimateCellTypeProportions <- function(eset, correlationThreshold = 0.7) {
  if (!requireNamespace("imsig", quietly = TRUE)) {
    stop("Package 'imsig' is required but not installed.")
  }

  # For testing
  #eset <- GSE53552_corrected_eset_unique_genes
  #correlationThreshold <- 0.7

  # Convert the ExpressionSet to a suitable format for ImSig
  exprData <- as.data.frame(Biobase::exprs(eset))
  rownames(exprData) <- rownames(Biobase::exprs(eset))

  # Ensure genes are HGNC symbols and no duplicates or missing values

  # Run ImSig to estimate cell type proportions
  cellProportions <- imsig::imsig(exp = exprData, r = correlationThreshold)

  # Change NANs to 0s
  cellProportions[is.na(cellProportions)] <- 0

  # If additional processing of cellProportions is needed, do it here
  # Assuming cellProportions is directly in the desired format based on the package's documentation

  return(cellProportions)
}
