#' Plot Heatmap of Selected Genes with Metadata Annotations
#'
#' Generates a heatmap for a subset of genes from an ExpressionSet object, using sample metadata
#' for column annotations. This function allows for detailed visualization of expression data alongside
#' relevant experimental context.
#'
#' @param eset An ExpressionSet object.
#' @param genes A character vector of gene names to be displayed in the heatmap.
#' @param metadataCols A vector indicating which columns of the metadata to use for column annotation.
#'        Can be a numeric vector specifying column indices or a character vector specifying column names.
#' @param scale A character string indicating if the data should be scaled; options are "none", "row", or "column".
#' @param cluster_rows Logical indicating whether to cluster rows.
#' @param cluster_cols Logical indicating whether to cluster columns.
#' @param show_rownames Logical indicating whether to show row names.
#' @param show_colnames Logical indicating whether to show column names.
#' @param ... Additional arguments passed to `pheatmap::pheatmap()`.
#'
#' @return A heatmap plot.
#' @importFrom Biobase exprs pData
#' @importFrom pheatmap pheatmap
#' @examples
#' # Assuming `eset` is your ExpressionSet and `selectedGenes` is a vector of gene names:
#' # plotHeatmap(eset, genes = c("Gene1", "Gene2"), metadataCols = c("Condition", "Treatment"))
#' @export
plotExpressionHeatmap <- function(eset,
                                  genes,
                                  metadataCols,
                                  scale = "none",
                                  cluster_rows = TRUE,
                                  cluster_cols = TRUE,
                                  show_rownames = TRUE,
                                  show_colnames = TRUE,
                                  color_limits = c(-3, 3),
                                  ...) {
  # Ensure necessary packages are installed
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("The 'pheatmap' package is required but not installed.")
  }

  eset <- GSE53552_corrected_eset
  genes <- topVariableGenes
  metadataCols <- columns_to_select[1:3]
  scale <- "row"
  cluster_rows <- TRUE
  cluster_cols <- TRUE
  annotation_row <- NULL
  annotation_col <- NULL
  show_rownames <- TRUE
  show_colnames <- TRUE
  color_limits <- c(-3, 3)

  # Extract expression data
  exprData <- exprs(eset)

  # Subset expression data for selected genes
  if (!all(genes %in% rownames(exprData))) {
    stop("Some specified genes are not found in the ExpressionSet object.")
  }
  exprDataSubset <- exprData[genes, ]

  # Extract metadata for column annotations
  metaData <- pData(eset)
  if (!is.null(metadataCols)) {
    if (is.numeric(metadataCols)) {
      metaData <- metaData[, metadataCols]
    } else if (is.character(metadataCols)) {
      metaData <- metaData[, metadataCols, drop = FALSE]
    }
  }

  # Prepare the annotation as a data frame where each column is a variable to annotate the columns of the heatmap
   annotation_metadata <- as.data.frame((metaData))

  # Custom color scale from royal blue to fire brick red
  color_map <- colorRampPalette(c("royalblue", "white", "firebrick"))(100)


  # pheatmap
  pheatmap::pheatmap(mat = exprDataSubset,
                     scale = scale,
                     annotation_col = annotation_metadata,
                     cluster_rows = cluster_rows,
                     cluster_cols = cluster_cols,
                     show_rownames = show_rownames,
                     show_colnames = show_colnames,
                     color = color_map,
                     breaks = seq(from = min(color_limits), to = max(color_limits), length.out = 101)
  )


}
