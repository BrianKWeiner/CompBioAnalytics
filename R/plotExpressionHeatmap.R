#' Plot Heatmap of Omics Data
#'
#' This function generates a heatmap of omics data (e.g., gene expression, protein abundances) using the pheatmap package.
#' It offers various customization options including scaling, clustering, and annotations.
#'
#' @param data A numeric matrix or data frame with rows as features (genes, proteins) and columns as samples or conditions.
#' @param scale A character string indicating if the data should be scaled; options are "none", "row", or "column".
#' @param cluster_rows Logical indicating whether to cluster rows.
#' @param cluster_cols Logical indicating whether to cluster columns.
#' @param annotation_row A data frame of annotations to be added to the rows, or NULL if not used.
#' @param annotation_col A data frame of annotations to be added to the columns, or NULL if not used.
#' @param show_rownames Logical indicating whether to show row names.
#' @param show_colnames Logical indicating whether to show column names.
#' @param ... Additional arguments passed to `pheatmap::pheatmap()`.
#'
#' @examples
#' # Assuming `expr_matrix` is a matrix of gene expression data
#' plotExpressionHeatmap(expr_matrix, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE)
#'
#' @return A heatmap plot.
#' @importFrom pheatmap pheatmap
#' @export
plotExpressionHeatmap <- function(data,
                                  scale = "none",
                                  cluster_rows = TRUE,
                                  cluster_cols = TRUE,
                                  annotation_row = NULL,
                                  annotation_col = NULL,
                                  show_rownames = TRUE,
                                  show_colnames = TRUE,
                                  color_limits = c(-3, 3), ...) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("The 'pheatmap' package is required but not installed.")
  }

  # Custom color scale from royal blue to fire brick red
  color_map <- colorRampPalette(c("royalblue", "white", "firebrick"))(100)

  # Truncate data to be within specified color limits
  data_truncated <- pmin(pmax(data, min(color_limits)), max(color_limits))

  # Prepare arguments for pheatmap
  args <- list(
    mat = data_truncated,
    scale = scale,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    annotation_row = annotation_row,
    annotation_col = annotation_col,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    color = color_map,
    breaks = seq(from = min(color_limits),
                 to = max(color_limits),
                 length.out = 101),
    ...
  )

  # Generate heatmap
  do.call(pheatmap::pheatmap, args)
}
