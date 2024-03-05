#' Perform ComBat Batch Correction, Visualize PCA, and Return Corrected Data
#'
#' @param eset An ExpressionSet object used for analysis.
#' @param batch A vector specifying the batch for each sample.
#' @param mod A model matrix specifying covariates of interest (other than batch).
#' @param outputPath Optional path for saving the PCA plot.
#' @return A matrix of batch-corrected expression data.
#' @importFrom sva ComBat
#' @importFrom ggplot2 ggplot geom_point aes labs ggsave
#' @importFrom stats prcomp
#' @examples
#' # eset is your ExpressionSet, batch is your batch vector
#' correctedData <- performComBatBatchCorrectionPCA(eset, batch, outputPath="PCA_Batch_Correction.png")
#' @export
performComBatBatchCorrectionPCA <- function(eset, batch, mod=NULL, outputPath="PCA_Batch_Correction.png") {
  if (!requireNamespace("sva", quietly = TRUE) || !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Required packages 'sva' and 'ggplot2' are not installed.")
  }


  #For testing
  #eset <- GSE53552_rma_genes
  #batch <- batch_info
  #mod <- NULL
  #outputPath <- "PCA_Batch_Correction.png"

  exprData <- Biobase::exprs(eset)

  # Ensure the order of batchInfo matches the order of samples in affyData
  sampleNames <- colnames(exprData)
  batchData <- batch[match(sampleNames, batch$FileName), ]

  # Factorize the batch data
  batchData <- factor(batchData$Batch)

  # PCA before batch correction
  pcaBefore <- prcomp(t(exprData))

  # Apply ComBat batch correction
  correctedData <- ComBat(dat=exprData, batch=batchData, mod=mod, par.prior=TRUE, prior.plots=FALSE)

  # PCA after batch correction
  pcaAfter <- prcomp(t(correctedData))

  # Prepare PCA data for plotting
  pcaBeforeDF <- data.frame(PC1 = pcaBefore$x[,1], PC2 = pcaBefore$x[,2], Batch = batch, Stage = 'Before')
  pcaAfterDF <- data.frame(PC1 = pcaAfter$x[,1], PC2 = pcaAfter$x[,2], Batch = batch, Stage = 'After')
  pcaDF <- rbind(pcaBeforeDF, pcaAfterDF)
  colnames(pcaDF) <- c("PC1", "PC2", "FileName", "OriginalDateTime", "FormattedDateTime", "Batch", "Stage")
  pcaDF$Batch <- factor(pcaDF$Batch)
  pcaDF$Stage <- factor(pcaDF$Stage, ordered = T, levels = c("Before", "After"))


  # Plot
  p <- ggplot(pcaDF, aes(x = PC1, y = PC2, color = Batch)) +
    geom_point() +
    facet_wrap(~ Stage, ncol = 2) +
    labs(title = "PCA Before and After Batch Correction", x = "PC1", y = "PC2") +
    theme_minimal() +
    scale_color_viridis_d()

  # Display the plot
  print(p)

  # Save the plot
  if (!is.null(outputPath)) {
    ggsave(outputPath, plot = p, device = "png", width = 12, height = 6)
  }

  # Return the corrected expression data
  return(correctedData)
}
