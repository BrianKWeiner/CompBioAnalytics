#' Perform Differential Expression Analysis with Explicit Contrast
#'
#' Conducts differential expression analysis comparing two specified groups within the ARM metadata column
#' of an ExpressionSet object. Utilizes the limma package with an explicit contrast setup (group1 - group2).
#'
#' @param eset An `ExpressionSet` object containing expression data and metadata.
#' @param group1 The name of the first group in the comparison.
#' @param group2 The name of the second group in the comparison.
#' @return A data frame containing the differential expression analysis results for the specified contrast.
#' @importFrom Biobase exprs pData
#' @importFrom limma lmFit makeContrasts contrasts.fit eBayes topTable
#' @examples
#' # Assuming `eset` is your ExpressionSet object and ARM contains groups 'Group1' and 'Group2':
#' # results <- performDEAnalysisWithContrast(eset, 'Group1', 'Group2')
#' @export
performDGEAnalysisByLimma <- function(eset, group1, group2) {
  requireNamespace("Biobase", quietly = TRUE)
  requireNamespace("limma", quietly = TRUE)

  # Validate input
  if (!inherits(eset, "ExpressionSet")) {
    stop("Input must be an ExpressionSet object.")
  }

  #eset <- GSE53552_with_arms
  #group1 <- unique_arms_GSE53552[5]
  #group2 <- unique_arms_GSE53552[6]

  # Access expression data and ARM data
  exprData <- exprs(eset)
  armData <- pData(eset)$ARM

  # Ensure specified groups exist in ARM
  if (!(group1 %in% armData) || !(group2 %in% armData)) {
    stop("Specified groups not found in ARM metadata.")
  }

  # Create a factor vector for the comparison
  groupFactor <- factor(armData, levels = c(group1, group2), ordered = FALSE)

  groupFactor <- factor(armData, levels = unique(eset@phenoData@data$ARM), ordered = FALSE)



  # Design matrix for limma analysis
  design <- model.matrix(~ groupFactor)

  colnames(design) <- gsub("groupFactor", "", colnames(design))
  colnames(design) <- gsub("\\(Intercept\\)", "Intercept", colnames(design))
  colnames(design)

  # Specify the contrast explicitly (group1 - group2)
  cmd <- paste0("contrastMatrix <- limma::makeContrasts(", group1, "-", group2, ", levels = design)")
  eval(parse(text = cmd))

  # Fit linear model using limma
  fit <- lmFit(exprData, design)

  # Set contrasts
  fit <- limma::contrasts.fit(fit, contrasts = contrastMatrix)

  # Apply empirical Bayes smoothing
  fit <- limma::eBayes(fit)

  # Extract the top table of results for the specified contrast
  resultsTable <- limma::topTable(fit,
                                  coef = colnames(contrastMatrix)[1],
                                  number = nrow(exprData),
                                  sort.by = "p")

  return(resultsTable)
}
