#' Create Tile Plot for ANOVA Results
#'
#' Visualizes the significance of ANOVA results between different arms for various cell types
#' using a tile plot. Each tile's color intensity represents the negative log (base 10) of
#' the adjusted p-value, highlighting the significance of differences observed.
#'
#' @param finalResults Dataframe containing ANOVA results, including columns for ARM1, ARM2,
#' CellType, and neg.log.10.p.adj, which is the negative logarithm (base 10) of the
#' adjusted p-value. This dataframe is expected to be the output from `analyzeCellTypeANOVA`.
#'
#' @return A ggplot object displaying a tile plot, which can be further customized or directly
#' rendered in R graphical devices.
#'
#' @import ggplot2
#' @examples
#' # Assuming 'finalResults' contains the output from `analyzeCellTypeANOVA`:
#' plot <- createTilePlotCellTypeAnovaResults(finalResults)
#' print(plot)
#'
#' @export
createTilePlotCellTypeAnovaResults <- function(finalResults) {
  if(nrow(finalResults) == 0) {
    stop("finalResults is empty. Tile plot cannot be created.")
  }

  finalResults$ArmComparison <- paste(finalResults$Arm1, finalResults$Arm2, sep = " : ")

  p <- ggplot(finalResults, aes(x = CellType, y = ArmComparison, fill = neg.log.10.p.adj)) +
    geom_tile(color = "black", size = 0.5) +  # Add black borders with size 1.5 around tiles
    scale_fill_gradient(low = "blue", high = "red", name = "-log10(p.adj)") +
    labs(title = "Tile Plot of Cell Type Significance Across ARM Comparisons",
         x = "Cell Type", y = "ARM Comparison") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank(),  # Remove minor grid lines
          panel.background = element_rect(color = "black", fill = NA))  # Add border to empty (white) tiles

  return(p)
}
