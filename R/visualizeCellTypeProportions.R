#' Visualize Cell Type Proportions
#'
#' @param cellTypeDF A dataframe containing cell type proportions per sample,
#' expected to be the output of `estimateCellTypeProportions`.
#' @return A ggplot object showing faceted bar graphs for each cell class.
#' @importFrom ggplot2 ggplot geom_bar aes labs facet_wrap
#' @examples
#' # Assuming 'cellTypeDF' is your dataframe from `estimateCellTypeProportions`:
#' visualizeCellTypeProportions(cellTypeDF)
#' @export
visualizeCellTypeProportions <- function(cellTypeDF) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed.")
  }

  # For testing
  #cellTypeDF <- GSE53552_proportion_cell_types

  cellTypeDF$SampleID <- rownames(cellTypeDF)

  # Convert wide format to long format for ggplot2
  longDF <- reshape2::melt(cellTypeDF, id.vars = "SampleID", variable.name = "CellType", value.name = "Proportion")

  # Create the plot
  p <- ggplot(longDF, aes(x = SampleID, y = Proportion, fill = CellType)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Cell Type Proportions Across Samples", x = "Sample ID", y = "Proportion") +
    facet_wrap(~ CellType, scales = "free_y") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) # Rotate x labels for readability

  return(p)
}
