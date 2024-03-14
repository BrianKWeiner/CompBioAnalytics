#' Filter Low Expressed Genes from an ExpressionSet
#'
#' This function filters out genes from an ExpressionSet object based on
#' RNA-seq count values. It retains only genes with a specified minimum number
#' of samples having counts greater than a specified cutoff.
#'
#' @param eset An ExpressionSet object containing RNA-seq count values.
#' @param minSamples The minimum number of samples in which the gene must have
#'        counts greater than the cutoff value to be retained. Default is 2.
#' @param countCutoff The count value above which counts are considered
#'        significant for the purpose of filtering. Default is 10.
#'
#' @return A new ExpressionSet object with filtered expression data and
#'         corresponding featureData.
#'
#' @importFrom Biobase ExpressionSet exprs
#' @examples
#' # Assuming 'eset' is your ExpressionSet object with RNA-seq count data:
#' filteredEset <- filterLowExpressedGenesRnaSeq(eset, minSamples = 2, countCutoff = 10)
#'
#' @export
filterLowExpressedGenesRnaSeq <- function(eset, minSamples = 2, countCutoff = 10) {
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Package 'Biobase' is needed but not installed.")
  }

  # For testing
  #eset <- GSE145513_eset
  #minSamples <- 2
  #countCutoff <- 10

  # Extract expression data from the ExpressionSet
  exprData <- Biobase::exprs(eset)

  # Apply filtering criterion
  keepGenes <- apply(exprData, 1, function(geneCounts) {
    sum(geneCounts > countCutoff) >= minSamples
  })

  # Subset the expression data and featureData
  filteredExprData <- exprData[keepGenes, ]
  filteredFeatureData <- if (!is.null(eset@featureData)) {
    eset@featureData@data[keepGenes, , drop = FALSE]
  } else {
    data.frame()
  }

  featureData <- new("AnnotatedDataFrame", data = filteredFeatureData)
  phenoData <- eset@phenoData
  annotationData <- eset@annotation

  # Construct and return new ExpressionSet
  newEset <- Biobase::ExpressionSet(assayData = filteredExprData,
                                    phenoData = phenoData,
                                    featureData = featureData,
                                    annotation = annotationData)



  return(newEset)
}
