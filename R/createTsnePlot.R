#' Produce t-SNE plots from ExpressionSet object
#'
#' This function creates t-SNE plots for visualizing high-dimensional expression data
#' contained in an ExpressionSet object, color-coded by specified phenotype columns.
#' Each phenotype column generates a separate plot, and the function returns a faceted
#' plot combining these individual plots.
#'
#' @param eset An ExpressionSet object containing expression data in the assayData slot
#' and phenotype data in the phenoData slot.
#' @param phenotypeCols A character vector specifying which columns from the phenotype
#' data to use for color-coding the t-SNE plots.
#' @param perplexity Perplexity parameter for the t-SNE algorithm (default is 30).
#' @return A ggplot object displaying a faceted plot of t-SNE visualizations, each facet
#' corresponding to a different phenotype variable for coloring.
#' @importFrom Rtsne Rtsne
#' @importFrom ggplot2 ggplot geom_point aes facet_wrap
#' @importFrom cowplot plot_grid
#' @examples
#' # eset is your ExpressionSet object with expression and phenotype data
#' # phenotypeCols might include columns like "Treatment" or "CellType"
#' tsneFacetedPlot <- createTsnePlot(eset, c("Treatment", "CellType"))
#'
#' @export
createTsnePlot <- function(eset, phenotypeCols, perplexity = 30) {
  require(Rtsne)
  require(ggplot2)
  require(cowplot)

  # For testing
  #eset <- GSE53552_corrected_eset
  #phenotypeCols <- c("source_name_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2", "characteristics_ch1.3", "characteristics_ch1.4", "batch", "PredictedGender")

  perplexity <- 30

  exprsData <- Biobase::exprs(eset)
  phenoData <- Biobase::pData(eset)

  # Perform t-SNE
  tsneRes <- Rtsne(t(exprsData), perplexity = perplexity, is_distance = FALSE)

  # Prepare data for plotting
  plotData <- data.frame(tsneRes$Y, phenoData)

  # Ensure all phenotypeCols data are treated as factors
  plotData[phenotypeCols] <- lapply(plotData[phenotypeCols], factor)

  plotList <- list()
  for (col in phenotypeCols) {
    if (!col %in% names(phenoData)) {
      stop(paste("Column", col, "not found in phenotype data."))
    }

    plotList[[col]] <- ggplot(plotData, aes_string(x = "X1", y = "X2", color = col)) +
      geom_point() +
      labs(x = "t-SNE 1", y = "t-SNE 2") +
      theme_minimal() +
      theme(legend.position = "none")
  }

  # Determine layout for combined plot
  nCols <- min(length(phenotypeCols), 4)  # Use up to 4 columns

  # Combine plots
  combinedPlot <- plot_grid(plotlist = plotList, ncol = nCols, align = 'v', labels = names(plotList))
  combinedPlot

  return(combinedPlot)
}
