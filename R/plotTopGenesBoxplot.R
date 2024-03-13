#' Plot Top Up- and Down-Regulated Genes from Limma Analysis
#'
#' This function generates two sets of box plots for the top 20 up- and down-regulated genes based on differential
#' expression analysis results from Limma. It showcases log2 expression values across specified comparison
#' groups. Genes are categorized as up- or down-regulated based on their T-statistics, and plots are annotated
#' with adjusted p-values.
#'
#' @param topTableData A dataframe output from Limma's topTable function, containing columns for logFC,
#'        adjusted p-values (`adj.P.Val`), T-statistics (`t`), and gene symbols.
#' @param eset An `ExpressionSet` object used in the Limma analysis, containing the expression data and
#'        phenotype data of samples.
#' @param group1 A string specifying the first comparison group or condition in the Limma contrast.
#' @param group2 A string specifying the second comparison group or condition in the Limma contrast.
#'
#' @return A list containing two ggplot objects: one for the top up-regulated genes and another for the top
#'         down-regulated genes, each saved as a PNG file.
#'
#' @importFrom ggplot2 ggplot geom_boxplot facet_wrap labs theme ggtitle ggsave
#' @importFrom Biobase exprs pData
#' @importFrom reshape2 melt
#'
#' @examples
#' topTableData <- limma::topTable(fit, coef="conditiontreated", number=Inf)
#' eset <- yourExpressionSet
#' plotTopGenesBoxplot(topTableData, eset, "Condition1", "Condition2")
#' @export

plotTopGenesBoxplot <- function(topTableData, eset, group1, group2) {
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("Biobase", quietly = TRUE)

  #Testing
  #topTableData <- topTable
  #eset <- GSE145513_eset_unique_genes
  #group1 <- unique_arms_GSE145513[1]
  #group2 <- unique_arms_GSE145513[2]

  # Split the data into up- and down-regulated genes based on logFC
  upRegulated <- topTableData[topTableData$logFC > 0, ]
  downRegulated <- topTableData[topTableData$logFC < 0, ]

  # Rank each group by ascending p-value and select the top 10
  topUp <- head(upRegulated[order(upRegulated$P.Val, decreasing = FALSE),], 20)
  topDown <- head(downRegulated[order(downRegulated$P.Val, decreasing = FALSE),], 20)

  # Combine the top 20 up- and down-regulated genes
  topGenes <- rbind(topUp, topDown)

  # Extract expression data for selected genes
  exprData <- Biobase::exprs(eset)[rownames(eset) %in% rownames(topGenes), ]

  # Extract expression data for top genes
  exprData <- exprs(eset)[rownames(exprs(eset)) %in% rownames(topGenes), ]

  # Rearrange the exprDATA so it is in the same order as topGenes
  exprData <- exprData[match(rownames(topGenes), rownames(exprData)),]

  exprData <- data.frame(exprData)

  # Add on gene symbol, p-values, and adusted p-values
  exprData$GeneSymbol <- rownames(exprData)
  exprData$P.value <- topGenes$P.Value
  exprData$adj.P.Val <- topGenes$adj.P.Val
  exprData$t.stat <- topGenes$t
  exprDataLong <- reshape2::melt(exprData, id = c("GeneSymbol", "P.value", "adj.P.Val", "t.stat"))

  colnames(exprDataLong) <- c("GeneSymbol", "P.value", "adj.P.Val", "t.stat", "Sample", "Log2Exp")

  # Add ARM information
  exprDataLong$ARM <- pData(eset)[exprDataLong$Sample, "ARM"]

  # Filter to remove the arms not part of the original contrast
  exprDataLongFilteredArm <- exprDataLong[exprDataLong$ARM %in% c(group1, group2),]

  exprDataLong <- exprDataLongFilteredArm

  group1Clean <- gsub("_", " ", group1)
  group1Clean <- gsub("(.{15,}?)\\s", "\\1\n", group1Clean)

  group2Clean <- gsub("_", " ", group2)
  group2Clean <- gsub("(.{15,}?)\\s", "\\1\n", group2Clean)

  # Clean ARM information: replace underscores with spaces, and break long strings
  exprDataLong$ARM <- gsub("_", " ", exprDataLong$ARM)
  exprDataLong$ARM <- gsub("(.{15,}?)\\s", "\\1\n", exprDataLong$ARM)

  exprDataLong$ARM <- factor(exprDataLong$ARM, levels = c(group1Clean, group2Clean), ordered = T)

  # Format p-values and  adjusted p-values for better legibility
  exprDataLong$P.value <- ifelse(nchar(as.character(exprDataLong$P.value)) > 6,
                                   sprintf("%.2e", exprDataLong$P.value),
                                   exprDataLong$P.value)
  exprDataLong$adj.P.Val <- ifelse(nchar(as.character(exprDataLong$adj.P.Val)) > 6,
                                   sprintf("%.2e", exprDataLong$adj.P.Val),
                                   exprDataLong$adj.P.Val)


  format_p_value <- function(adj.p.value) {
    ifelse(adj.p.value < 0.05,
           sprintf("adj. p: %s", sub("e", " x 10^", format(adj.p.value, scientific = TRUE))),
           sprintf("adj. p: %s", format(adj.p.value, scientific = FALSE)))
  }

  # Split the data based on the sign of the T-statistic
  upRegulated <- exprDataLong[exprDataLong$t.stat > 0, ]
  downRegulated <- exprDataLong[exprDataLong$t.stat < 0, ]

  # Function to create plots
  create_plot <- function(data, title, subtitle) {
    data <- upRegulated
    # Use gene symbol as facet labels and include formatted adjusted p-value
    data$FacetLabel <- paste(data$GeneSymbol, "\n(", lapply(data$adj.P.Val, format_p_value), ")", sep="")

    p <- ggplot(data, aes(x = ARM, y = Log2Exp, fill = ARM)) +
      ggtitle(title, subtitle = subtitle) +
      geom_boxplot() +
      facet_wrap(~ FacetLabel, scales = "free", ncol = 5) +
      labs(y = "Log2 Expression", x = "") +
      theme(strip.text.x = element_text(size = 10),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none")
    return(p)
  }

  contrast <- paste0(group1, " vs ", group2)
  contrast_no_space <- paste0(group1, "_vs_", group2)

  # Create and save plots for up-regulated genes
  pUp <- create_plot(upRegulated, "Top Up-Regulated Genes (by p-value):", contrast)
  ggsave(paste0("Top_Up_Regulated_Genes.",contrast_no_space,".png"), plot = pUp, device = "png", width = 14, height = 6.2)

  # Create and save plots for down-regulated genes
  pDown <- create_plot(downRegulated, "Top Down-Regulated Genes (by p-value):", contrast)
  ggsave(paste0("Top_Down_Regulated_Genes.",contrast_no_space,".png"), plot = pUp, device = "png", width = 14, height = 6.2)


  # Return the plots
  return(list(UpRegulated = pUp, DownRegulated = pDown))
}
