#' Perform Pre-Ranked GSEA Using fgsea with MSigDB Category Selection
#'
#' Performs GSEA on limma topTable results using fgsea, allowing users to specify the gene set category from MSigDB.
#'
#' @param topTableData Dataframe from topTable from limma DEG analysis that contains p-value, fold change, and t-stat info
#' @param rankingMethod Method to rank genes: 't_stat', 'fc', or 'signed_p'.
#' @param msigdbCategory Category of MSigDB gene sets to use, e.g., 'H' for Hallmark.
#' @param minSize Minimum size of the gene set to be tested.
#' @param maxSize Maximum size of the gene set to be tested.
#' @param nperm Number of permutations for estimating p-values.
#' @return A list containing the fgsea results for the specified MSigDB category.
#' @importFrom fgsea fgsea
#' @importFrom msigdbr msigdbr
#' @examples
#' # Assuming 'topTableData' is the result of :
#' # GSEAresults <- performGSEA(topTableData, rankingMethod = 't_stat', msigdbCategory = 'H',
#' #                            minSize = 10, maxSize = 500, nperm = 10000)
#' @export
performGSEA <- function(fgseaRes, rankingMethod, msigdbCategory, minSize = 10, maxSize = 500, nperm = 10000) {
  requireNamespace("fgsea", quietly = TRUE)
  requireNamespace("msigdbr", quietly = TRUE)

  # Validate input
  if (!"data.frame" %in% class(topTableData)) {
    stop("Input must be a dataframe output from limma's topTable function.")
  }

  #For testing
  #topTableData <- topTable
  #msigdbCategory = 'H'
  #rankingMethod = 't_stat'
  #minSize = 10
  #maxSize = 500
  #nperm = 10000

  # Set seed for reproducibility
  set.seed(12345678)

  # Fetch gene sets based on the specified MSigDB category
  msigdbGeneSets <- msigdbr::msigdbr(species = "Homo sapiens", category = msigdbCategory)
  geneSets <- split(msigdbGeneSets$gene_symbol, msigdbGeneSets$gs_name)


  # Prepare the ranked list based on the specified ranking method
  rankedList <- switch(rankingMethod,
                       t_stat = setNames(topTableData$t, rownames(topTableData)),
                       fc = setNames(topTableData$logFC, rownames(topTableData)),
                       signed_p = {
                         signFC <- sign(topTableData$logFC)
                         pVals <- topTableData$P.Value
                         setNames(signFC * (1/pVals), rownames(topTableData))
                       },
                       stop("Invalid ranking method. Choose 't_stat', 'fc', or 'signed_p'."))


  # Perform fgsea
  fgseaRes <- fgsea::fgseaSimple(pathways = geneSets,
                           stats = unlist(rankedList),
                           minSize = minSize,
                           maxSize = maxSize,
                           nperm = nperm)

  return(fgseaRes)
}
