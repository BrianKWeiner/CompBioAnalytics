#' Generate Report on Phenotypic Data in ExpressionSet
#'
#' Analyzes the phenotypic data (`phenoData`) of an `ExpressionSet` object,
#' providing counts of each unique value for each column. For columns with all unique values,
#' only the first three values are listed followed by '...'. Generates a text report
#' summarizing these counts, useful for understanding the composition of phenotypic variables.
#'
#' @param eset An `ExpressionSet` object.
#' @return A character vector, each element a line of the report, detailing
#'         the counts of unique values for each phenotypic data column.
#' @importFrom Biobase phenoData
#' @examples
#' # Assuming `eset` is your ExpressionSet object:
#' report <- generatePhenoDataReport(eset)
#' cat(report, sep = "\n")
#' @export
generatePhenoDataReport <- function(eset) {
  requireNamespace("Biobase", quietly = TRUE)

  # Access the phenoData data frame
  phenoData <- pData(eset)

  # Initialize an empty vector to store the report lines
  report <- character()

  # Iterate through each column in phenoData
  for (colName in names(phenoData)) {
    # Calculate the frequency of each unique value in the current column
    valueCounts <- table(phenoData[[colName]])

    # Check if all values are unique
    allUnique <- length(valueCounts) == nrow(phenoData)

    # Create report lines for the current column
    colReport <- paste("Column", colName, "has", length(valueCounts), "unique values:")

    # Adjust the report for columns with all unique values
    if (allUnique && length(valueCounts) > 3) {
      valueReports <- c(paste(" -", names(valueCounts)[1:3], ":", valueCounts[1:3]), " - ...")
    } else {
      valueReports <- paste(" -", names(valueCounts), ":", valueCounts)
    }

    # Combine the column report with its value reports
    report <- c(report, colReport, valueReports)
  }

  return(report)
}
