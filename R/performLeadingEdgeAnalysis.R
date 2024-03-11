#' Perform Leading Edge Analysis from fgsea Results
#'
#' Visualizes the most frequent genes in the leading edges of upregulated or downregulated gene sets from fgsea results using a tile plot.
#'
#' @param fgseaRes Data frame containing fgsea results.
#' @param direction Character specifying the direction of regulation ("up" or "down").
#' @param maxGeneSize Integer specifying the maximum number of leading edge gene to display. Default is 20.
#' @param maxGeneSets Integer specifying the maximum number of gene sets to display. Default is 25.
#' @importFrom ggplot2 ggplot geom_tile aes labs theme
#' @examples
#' # Assuming 'fgseaRes' is your dataframe containing fgsea results:
#' performLeadingEdgeAnalysis(fgseaRes, direction = "up", maxGeneSize = 20, maxGeneSets = 25)
#' @export
performLeadingEdgeAnalysis <- function(fgseaRes, direction = "up", maxGeneSize = 20, maxGeneSets = 25) {
  require(ggplot2)

  # For testing
  #fgseaRes <- fgseaRes
  #direction <- "up"
  #maxGeneSize <- 20
  #maxGeneSets <- 25

  graphTitleFileName <- "Top_X_Regulated_Leading_Edge_Genes.png"

  # Filter fgsea results based on direction
  if (direction == "up") {
    fgseaSubset <- fgseaRes[fgseaRes$NES > 0, ]
    graphTitleFileName <- "Top_Up_Regulated_Leading_Edge_Genes.png"
  } else if (direction == "down") {
    fgseaSubset <- fgseaRes[fgseaRes$NES < 0, ]
    graphTitleFileName <- "Top_Down_Regulated_Leading_Edge_Genes.png"
  } else {
    stop("Invalid direction. Please choose 'up' or 'down'.")
  }

  # Sort and select top gene sets based on NES
  fgseaSubset <- fgseaSubset[order(fgseaSubset$NES, decreasing = TRUE), ]

  # Extract leading edge genes
  leadingEdges <- fgseaSubset$leadingEdge
  names(leadingEdges) <- fgseaSubset$pathway

  # Count gene frequencies
  geneFrequency <- table(unlist(leadingEdges))

  # Select top 20 most frequent genes
  topGenes <- names(sort(geneFrequency, decreasing = TRUE))[1:20]

  # Initialize a dataframe to hold the presence (1) or absence (0) of top genes in each pathway's leading edge
  plotData <- matrix(0,
                     nrow = length(fgseaSubset$pathway),
                     ncol = length(topGenes),
                     dimnames = list(fgseaSubset$pathway, topGenes))

  # Fill the matrix with presence/absence data
  for (i in seq_along(leadingEdges)) {
    for (gene in leadingEdges[[i]]) {
      if (gene %in% topGenes) {
        plotData[rownames(plotData) == names(leadingEdges)[i], gene] <- 1
      }
    }
  }

  # After converting the matrix to a dataframe but before melting it for ggplot
  plotData <- as.data.frame(plotData)
  plotData$Pathway <- rownames(plotData)

  # Calculate column sums for sorting
  columnSums <- rowSums(plotData[, -ncol(plotData)]) # Exclude the Pathway column
  plotData$Sum <- columnSums

  # Order the dataframe by column sums in descending order
  plotData <- plotData[order(-plotData$Sum), ]

  # Discard the pathways that do not have any of the top leading edge genes
  pathwaysToDropIndices <- which(plotData$Sum == 0)
  plotData <- plotData[-pathwaysToDropIndices,]

  # If remaining sets of genes is greater than maxGeneSets then trim it
  pathwaysThatRemain <- dim(plotData)[1][1]
  if (pathwaysThatRemain > maxGeneSets) {  plotData <- plotData[1:maxGeneSets,]}

  # Convert the Pathway column into a factor with levels set according to the order in the dataframe
  plotData$Pathway <- factor(plotData$Pathway, levels = plotData$Pathway, ordered = TRUE)

  # Drop the Sum column before melting
  plotData$Sum <- NULL

  # Melt the dataframe for ggplot
  plotDataMelted <- reshape2::melt(plotData, id.vars = "Pathway")

  # Create the tile plot with the ordered pathways
  leadingEdgeGene_plot <- ggplot(plotDataMelted, aes(x = Pathway, y = variable, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(title = "Leading Edge Analysis", x = "Gene Sets", y = paste0("Top ",maxGeneSize," Genes")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          legend.position = "none")

  print(leadingEdgeGene_plot)

  ggsave(graphTitleFileName, plot = leadingEdgeGene_plot, device = "png", width = 14, height = 6.2)

  return(leadingEdgeGene_plot)


}
