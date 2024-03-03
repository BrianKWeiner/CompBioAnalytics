#' Create Enrichment Plots for Top and Bottom Ranked Pathways
#'
#' Generates enrichment plots for the top 10 pathways with the highest and lowest NES
#' from the results of a GSEA analysis performed using the fgsea package.
#'
#' @param fgseaRes Results from fgsea analysis.
#' @param topTableData A dataframe containing differential expression analysis results.
#' Must include columns for log fold change (`logFC`) and adjusted p-value (`adj.P.Val`).
#' @importFrom fgsea plotEnrichment
#' @importFrom ggplot2 ggsave
#' @examples
#' # Assuming 'fgseaRes' contains the results from your fgsea analysis:
#' # createEnrichmentPlots(fgseaRes, topTableData)
#' @export
createEnrichmentPlots <- function(fgseaRes, topTableData) {
  requireNamespace("fgsea", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)

  #for testing
  #fgseaRes <- GSEAresults
  #topTableData <- topTable
  #rankingMethod <- "t_stat"

  # Sort fgsea results by NES
  fgseaResSorted <- fgseaRes[order(fgseaRes$NES), ]

  # Select top 10 pathways with highest and lowest NES
  topPathways <- tail(fgseaResSorted, 10)
  bottomPathways <- head(fgseaResSorted, 10)

  # Concatenate top and bottom pathways
  selectedPathways <- rbind(bottomPathways, topPathways)

  # Prepare the ranked list based on the specified ranking method
  rankedList <- switch(rankingMethod,
                       t_stat = setNames(topTableData$t, rownames(topTableData)),
                       fc = setNames(topTableData$logFC, rownames(topTableData)),
                       signed_p = {
                         signFC <- sign(topTableData$logFC)
                         pVals <- topTableData$P.Value
                         setNames(signFC * (1/pVals), rownames(topTableData))
                       },
                       stop("Invalid ranking method. Choose 't_stat', 'fc', or 'signed_p'.")
                       )

  exampleRanks <- rankedList


  msigdbGeneSets <- msigdbr::msigdbr(species = "Homo sapiens")
  geneSets <- split(msigdbGeneSets$gene_symbol, msigdbGeneSets$gs_name)

  # Generate enrichment plots
  plots <- lapply(selectedPathways$pathway, function(pathway) {

    cleanTitle <- gsub("_", " ", pathway)

    pd <- plotEnrichmentData(
      pathway = geneSets[[pathway]],
      stats = unlist(rankedList)
    )

    # Enhance it to make it look nicer
    with(pd,
         ggplot(data=curve) +
           ggtitle(cleanTitle) +
           geom_line(aes(x=rank, y=ES), color="green") +
           geom_ribbon(data=stats,
                       mapping=aes(x=rank, ymin=0,
                                   ymax=stat/maxAbsStat*(spreadES/4)),
                       fill="grey") +
           geom_segment(data=ticks,
                        mapping=aes(x=rank, y=-spreadES/16,
                                    xend=rank, yend=spreadES/16),
                        size=0.2) +
           geom_hline(yintercept=posES, colour="red", linetype="dashed") +
           geom_hline(yintercept=negES, colour="red", linetype="dashed") +
           geom_hline(yintercept=0, colour="black") +
           theme(
             panel.background = element_blank(),
             panel.grid.major=element_line(color="grey92")
           ) +
           labs(x="rank", y="enrichment score"))
    return(pd)
  })

  # Optionally, save plots to files
  lapply(seq_along(plots), function(i) {
    ggplot2::ggsave(filename = paste0("enrichment_plot_", i, ".png"),
                    plot = plots[[i]],
                    device = "png")
  })

  return(plots)
}
