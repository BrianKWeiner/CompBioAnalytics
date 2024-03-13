#' Update ExpressionSet with Unique Gene Names Based on Highest Median Expression
#'
#' Updates an ExpressionSet object to ensure unique gene names in assayData$exprs,
#' choosing the probe with the highest median expression for genes with multiple mappings.
#' Allows specifying the column name where gene symbols are stored.
#'
#' @param eset An ExpressionSet object with potentially non-unique gene names.
#' @param geneSymbolColumn The name of the column in featureData where gene symbols are stored.
#'        Default is "Gene Symbol". Can be set to another column name such as "external_gene_name".
#' @return An updated ExpressionSet object with unique gene names and probes selected based on
#'         the highest median expression.
#' @importFrom Biobase ExpressionSet exprs featureData
#' @importFrom dplyr group_by summarise filter slice ungroup select
#' @export
updateEsetWithUniqueGenesWithColumnSpecified <- function(eset, geneSymbolColumn = "Gene Symbol") {
  requireNamespace("Biobase", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)

  #for testing
  #eset <- GSE145513_eset_voom_transformed
  #geneSymbolColumn = "external_gene_name"

  # Extract expression data
  exprData <- exprs(eset)

  # Dynamically extract gene symbols based on the specified column
  geneSymbols <- sapply(eset@featureData@data[[geneSymbolColumn]], function(x) strsplit(x, " /// ")[[1]][1])


  # Calculate median expression for each probe
  medianExprs <- apply(exprData, 1, median)

  # Combine gene symbols with their corresponding median expression values into a dataframe
  geneInfo <- data.frame(geneSymbol = geneSymbols, medianExpr = medianExprs, stringsAsFactors = FALSE)
  rownames(geneInfo) <- rownames(exprData)

  geneSymbols <- geneSymbols[!is.na(geneSymbols)]

  geneInfo <- geneInfo[!is.na(geneInfo$geneSymbol),]

  # Initialize a list to store the selected probe for each unique gene
  selectedProbes <- vector("list", length = length(unique(geneSymbols)))

  # Iterate over unique gene symbols to find the probe with the highest median expression for each
  for (gene in unique(geneSymbols)) {
    probesForGene <- which(geneSymbols == gene)
    if (length(probesForGene) > 1) {
      # If multiple probes map to the same gene, select the one with the highest median expression
      highestExprProbe <- probesForGene[which.max(medianExprs[probesForGene])]
      selectedProbes[[gene]] <- rownames(geneInfo)[highestExprProbe]
    } else if (length(probesForGene) == 1) {
      # If only one probe maps to the gene, select it directly
      selectedProbes[[gene]] <- rownames(geneInfo)[probesForGene]
    }
  }


  # Get indices to filter the expression and featureData
  selectedIndices <- unique(unlist(selectedProbes))
  matchedIndices <- match(selectedProbesVec, rownames(eset@assayData$exprs))

  # Update expression and featureData in eset based on selected indices
  updatedExprs <- exprData[selectedIndices, ]
  updatedFeatureData <- eset@featureData@data[matchedIndices, , drop = FALSE]
  rownames(updatedExprs) <- names(unlist(selectedProbes))

  # Ensure unique row names by appending numeric suffix to duplicates
  rownames(updatedExprs) <- make.unique(rownames(updatedExprs))

  # Update rownames of the featureData to ensure consistency when building the
  # new eSet object.

  rownames(updatedFeatureData) <- rownames(updatedExprs)

  # Create a new ExpressionSet with updated data
  newEset <- Biobase::ExpressionSet(assayData = updatedExprs,
                                    phenoData = eset@phenoData,
                                    featureData = AnnotatedDataFrame(updatedFeatureData),
                                    experimentData = eset@experimentData,
                                    annotation = eset@annotation)

  return(newEset)
}
