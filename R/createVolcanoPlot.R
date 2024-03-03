#' Create Volcano Plot from Differential Expression Results
#'
#' Generates a volcano plot for visualizing differential expression analysis results,
#' highlighting genes based on log2 fold change and adjusted p-value thresholds. Allows
#' for customization of significant gene labeling, fold change, and p-value thresholds.
#'
#' @param topTableData A dataframe containing differential expression analysis results.
#' Must include columns for log fold change (`logFC`) and adjusted p-value (`adj.P.Val`).
#' @param labelGenes Logical, indicating whether to label significant genes on the plot.
#' Default is `TRUE`.
#' @param log2FCThreshold Numeric, specifying the log2 fold change threshold for
#' highlighting significant genes. Default is `1`.
#' @param pValThreshold Numeric, specifying the adjusted p-value threshold for
#' highlighting significant genes. Default is `0.05`.
#' @return An EnhancedVolcano plot object.
#' @importFrom EnhancedVolcano EnhancedVolcano
#' @examples
#' # Assuming `topTableData` is your dataframe from limma's topTable with required columns:
#' createVolcanoPlot(topTableData, labelGenes = TRUE, log2FCThreshold = 1, pValThreshold = 0.05)
#'
#' @export
createVolcanoPlot <- function(topTableData, labelGenes = TRUE, log2FCThreshold = 1, pValThreshold = 0.05) {
  requireNamespace("EnhancedVolcano", quietly = TRUE)

  if (!"data.frame" %in% class(topTableData)) {
    stop("Input must be a dataframe output from limma's topTable function.")
  }

  requiredCols <- c("logFC", "adj.P.Val")
  if (!all(requiredCols %in% names(topTableData))) {
    stop(paste("Input dataframe must include the following columns:", paste(requiredCols, collapse = ", ")))
  }

  # Simplify and reformat gene names if labeling is enabled
  formattedGeneNames <- if (labelGenes) {
    sapply(rownames(topTableData), function(name) {
      parts <- strsplit(name, "_", fixed = TRUE)[[1]]
      if (length(parts) > 1) {
        paste(tail(parts, 1), paste(parts[-length(parts)], collapse = "_"), sep = "\n")
      } else {
        name
      }
    })
  } else {
    rownames(topTableData)
  }

  # Determine significance based on log2FC and adjusted p-value thresholds
  significant <- with(topTableData, logFC > log2FCThreshold & adj.P.Val < pValThreshold | logFC < -log2FCThreshold & adj.P.Val < pValThreshold)
  sigGenes <- if (labelGenes) formattedGeneNames[significant] else character(0)

  EnhancedVolcano::EnhancedVolcano(
    topTableData,
    lab = formattedGeneNames,
    x = 'logFC',
    y = 'adj.P.Val',
    xlim = c(-3, 3),
    title = 'Differential Expression',
    xlab = 'Log2 Fold Change',
    ylab = '-Log10 Adjusted P-value',
    pCutoff = pValThreshold,
    FCcutoff = log2FCThreshold,
    pointSize = 2.0,
    labSize = 3.0,
    col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
    colAlpha = 1,
    legendLabels = c('Not Significant', 'Log2FC', 'Adjusted P-value', 'Log2FC & Adjusted P-value'),
    legendPosition = 'right',
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    selectLab = sigGenes # Only label significant genes if labeling is enabled
  )
}
