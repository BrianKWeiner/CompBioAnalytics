#' Perform ComBat Batch Correction
#'
#' Applies ComBat batch correction to an AffyBatch object's expression data,
#' using the batch information provided. This function relies on the `sva` package.
#'
#' @param affyData An `AffyBatch` object.
#' @param batchInfo A data frame with batch information for each sample, typically
#'        produced by the `categorizeCelDataIntoBatches` function. Must contain columns
#'        'FileName' and 'Batch'.
#' @return A matrix of batch-corrected expression data.
#' @importFrom sva ComBat
#' @importFrom Biobase exprs pData
#' @examples
#' # Assuming `affyData` is your AffyBatch object and `batchInfo` contains the batch data:
#' # correctedData <- performComBatBatchCorrection(affyData, batchInfo)
#' @export
performComBatBatchCorrection <- function(affyData, batchInfo) {
  requireNamespace("sva", quietly = TRUE)
  requireNamespace("Biobase", quietly = TRUE)

  # Extract expression data from AffyBatch object
  exprData <- exprs(affyData)

  # Ensure the order of batchInfo matches the order of samples in affyData
  sampleNames <- colnames(exprData)
  batchInfo <- batchInfo[match(sampleNames, batchInfo$FileName), ]

  # Perform ComBat correction
  batch <- factor(batchInfo$Batch)
  correctedData <- sva::ComBat(dat = exprData,
                               batch = batch,
                               mod = NULL,
                               par.prior = TRUE,
                               prior.plots = FALSE)

  return(correctedData)
}
