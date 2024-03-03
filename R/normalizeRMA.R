#' Perform RMA Normalization on AffyBatch
#'
#' This function applies Robust Multi-array Average (RMA) normalization to an AffyBatch object,
#' which involves background adjustment, quantile normalization, and summarization. This should
#' be performed before any batch correction methods.
#'
#' @param affyBatch An `AffyBatch` object to be normalized.
#' @return An `ExpressionSet` object containing the RMA-normalized expression values.
#' @importFrom affy rma
#' @examples
#' # Assuming `affyBatch` is your AffyBatch object:
#' # esetRMA <- normalizeRMA(affyBatch)
#' @export
normalizeRMA <- function(affyBatch) {
  requireNamespace("affy", quietly = TRUE)

  # Validate input
  if (!inherits(affyBatch, "AffyBatch")) {
    stop("The input must be an AffyBatch object.")
  }

  # Apply RMA normalization
  esetRMA <- affy::rma(affyBatch)

  return(esetRMA)
}
