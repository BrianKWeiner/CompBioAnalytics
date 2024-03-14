#' Apply Limma-Voom Transformation to ExpressionSet
#'
#' Transforms count data in an ExpressionSet object using the Limma-Voom method,
#' preparing it for linear modeling. This transformation is useful for RNA-seq data.
#' The function returns a new ExpressionSet object with expression data transformed using Limma-Voom.
#'
#' @param eset An ExpressionSet object containing raw count data.
#' @return A new ExpressionSet object with expression data transformed using Limma-Voom.
#' @examples
#' # Assuming 'eset' is your ExpressionSet object with RNA-seq count data:
#' eset_voom <- applyLimmaVoomTransformation(eset)
#'
#' @export
applyLimmaVoomTransformation <- function(eset) {
  # Ensure necessary packages are installed
  if (!requireNamespace("Biobase", quietly = TRUE)) stop("Package 'Biobase' is required but not installed.")
  if (!requireNamespace("limma", quietly = TRUE)) stop("Package 'limma' is required but not installed.")

  # For testing
  #eset <- GSE145513_eset

  # Extract count data from the ExpressionSet
  counts <- Biobase::exprs(eset)

  # Apply Limma-Voom transformation
  # Design matrix, assuming a basic case. Adjust as needed for your study.
  design <- model.matrix(~ 1, data = Biobase::pData(eset))
  v <- limma::voom(counts, design, plot = TRUE)

  # Create a new ExpressionSet with the transformed data
  eset_transformed <- new("ExpressionSet",
                          assayData = list(exprs = v$E),
                          phenoData = eset@phenoData,
                          featureData = eset@featureData,
                          experimentData = eset@experimentData,
                          annotation = Biobase::annotation(eset))

  return(eset_transformed)
}
