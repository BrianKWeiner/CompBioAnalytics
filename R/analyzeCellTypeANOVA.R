#' Perform ANOVA on Cell Type Proportions and Post-hoc Comparisons
#'
#' This function analyzes cell type proportions across different sample arms,
#' performs ANOVA and post-hoc Tukey tests to identify statistically significant
#' differences, and then summarizes the findings. It returns a list containing
#' both the detailed ANOVA results and a summary of cell types that showed
#' significant differences between arms.
#'
#' @param cellTypeDF A dataframe of cell type proportions per sample.
#' @param phenoData A dataframe from the @phenoData slot of an ExpressionSet,
#' containing sample arms.
#' @param pValueThreshold The threshold for significance in post-hoc tests
#' (default is 0.05).
#'
#' @return A list with two elements:
#' \itemize{
#'   \item{\code{finalResults}:}{A dataframe with Arm1, Arm2, CellType, p.value,
#'   and corrected p.value for each significant comparison found by ANOVA and
#'   post-hoc analysis.}
#'   \item{\code{summaryFinalResults}:}{A dataframe summarizing the significant
#'   cell types between each unique pair of Arm1 and Arm2. It contains three
#'   columns: Arm1, Arm2, and CellTypes, where CellTypes is a concatenated string
#'   of all the significant CellTypes for that arm pair.}
#' }
#'
#' @importFrom stats aov TukeyHSD p.adjust
#' @importFrom dplyr select filter group_by summarise ungroup
#' @examples
#' # Assuming 'cellTypeDF' is your cell type proportions and 'phenoData' contains the arms:
#' result <- analyzeCellTypeANOVA(cellTypeDF, phenoData, 0.05)
#' # Access the finalResults with result$finalResults
#' # Access the summaryFinalResults with result$summaryFinalResults
#' @export
analyzeCellTypeANOVA <- function(cellTypeDF, phenoData, pValueThreshold = 0.05) {
  results <- list()

  #For testing
  #cellTypeDF <- GSE53552_proportion_cell_types
  #phenoData <- GSE53552_with_arms@phenoData@data
  #pValueThreshold = 0.05
  #phenoData$SampleID <- phenoData$geo_accession

  if(!"SampleID" %in% names(phenoData)) phenoData$SampleID <- phenoData$geo_accession

  rearrangeColumns <- function(df) df[, c(ncol(df), 1:(ncol(df)-1))]
  phenoData <- rearrangeColumns(phenoData)
  cellTypeDF <- rearrangeColumns(cellTypeDF)

  colnames(cellTypeDF) <- gsub("[^[:alpha:]]", "_", colnames(cellTypeDF))


  # Ensure cellTypeDF has sample IDs that match phenoData
  cellTypeDF$SampleID <- rownames(cellTypeDF)
  mergedData <- merge(cellTypeDF, phenoData, by = "SampleID")

  # Ensure that cellTypeDF has different values in the column else discard.
  filterColumns <- function(cellTypeDF) {
    # Filter columns where all values are the same or all values are zeros
    cellTypeDF <- cellTypeDF[, sapply(cellTypeDF, function(x) length(unique(x)) > 1 & !all(x == 0))]
    return(cellTypeDF)
  }
  cellTypeDF <- filterColumns(cellTypeDF)


  # Loop through each cell type
  for (cellType in colnames(cellTypeDF)[-1]) {  # Assuming first column is SampleID
    if (cellType == "SampleID") next
    formula <- as.formula(paste(cellType, "~ ARM"))

    # Try-catch block to handle errors in ANOVA
    tryCatch({
      anovaResult <- aov(formula, data = mergedData)

      # Ensure ANOVA was successfully computed
      if (is.null(anovaResult)) {
        next
      }

      # Post-hoc test with error handling for empty or invalid ANOVA results
      posthoc <- tryCatch({
        TukeyHSD(anovaResult)
      }, error = function(e) {
        return(NULL)  # Return NULL if TukeyHSD fails
      })

      if (is.null(posthoc)) {
        next  # Skip to the next iteration if posthoc test fails
      }

      # Extracting significant comparisons
      comparisons <- posthoc[[1]][, 4] < pValueThreshold
      if (any(comparisons)) {
        sigComparisons <- as.data.frame(posthoc[[1]][comparisons, ])
        if (dim(sigComparisons)[2] == 1) {
          sigComparisons <- as.data.frame(t(posthoc[[1]][comparisons, ]))
        }
        rownames(sigComparisons) <- rownames(posthoc$ARM)[comparisons]

        # Prepare results
        sigResults <- data.frame(Arm1 = sub("-.*", "", rownames(sigComparisons)),
                                 Arm2 = sub(".*-", "", rownames(sigComparisons)),
                                 CellType = cellType,
                                 p.adj = sigComparisons[, 4],
                                 neg.log.10.p.adj = -log10(sigComparisons[, 4]) )

        results[[cellType]] <- sigResults
      }
    }, error = function(e) {
      warning(sprintf("ANOVA failed for cell type '%s': %s", cellType, e$message))
    })
  }

  # Combine all results into a single dataframe
  if (length(results) == 0) {
    return(data.frame(Arm1 = character(), Arm2 = character(), CellType = character(), p.adj = numeric(), neg.log.10.p.adj = numeric()))
  }

  finalResults <- do.call(rbind, results)

  # Order the finalResults dataframe by p.adj in ascending order
  finalResults <- finalResults[order(finalResults$p.adj), ]


  # Assuming finalResults is your dataframe with columns Arm1, Arm2, and CellType
  # Generate summaryFinalResults
  if(nrow(finalResults) > 1) {
    summaryFinalResults <- finalResults %>%
      group_by(Arm1, Arm2) %>%
      summarise(CellTypes = toString(unique(CellType))) %>%
      ungroup()
    summaryFinalResults <- as.data.frame(summaryFinalResults)
  } else if (nrow(finalResults) == 1) {
    # If only one row in finalResults, use that row to create summaryFinalResults
    summaryFinalResults <- finalResults %>%
      select(Arm1, Arm2, CellType) %>%
      mutate(CellTypes = toString(CellType)) %>%
      select(-CellType)
    summaryFinalResults <- as.data.frame(summaryFinalResults)
  } else {
    # If no rows in finalResults, return an empty dataframe with the expected structure
    summaryFinalResults <- data.frame(Arm1 = character(0), Arm2 = character(0), CellTypes = character(0))
  }

  # Return both finalResults and summaryFinalResults as a list
  return(list(finalResults = finalResults, summaryFinalResults = summaryFinalResults))
}
