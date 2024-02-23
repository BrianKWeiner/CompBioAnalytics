#' Transform Expression Data in ExpressionSet
#'
#' This function modifies the expression matrix within an ExpressionSet object by performing
#' the following steps:
#' 1. Replaces any values less than 0 with 0.
#' 2. Adds 1 to every value in the matrix.
#' 3. Applies a log2 transformation to all values.
#' The modified ExpressionSet object is returned.
#'
#' @param eset An `ExpressionSet` object.
#' @return A modified `ExpressionSet` object with the expression data transformed.
#' @importFrom Biobase exprs
#' @examples
#' # Assuming `eset` is your ExpressionSet object:
#' # esetTransformed <- transformExpressionData(eset)
#' @export
transformExpressionData <- function(eset) {
  # Access the expression data
  exprMatrix <- exprs(eset)

  # Step 1: Replace values < 0 with 0
  exprMatrix[exprMatrix < 0] <- 0

  # Step 2: Add 1 to every value
  exprMatrix <- exprMatrix + 1

  # Step 3: Apply log2 transformation
  exprMatrix <- log2(exprMatrix)

  # Update the expression data in the eset
  exprs(eset) <- exprMatrix

  # Return the modified eset
  return(eset)
}
