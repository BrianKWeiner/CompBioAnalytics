#' Read Affymetrix CEL Files from Directory
#'
#' Reads all Affymetrix .CEL files located in a specified directory and returns an `AffyBatch` object
#' containing the raw data from these files. This function requires the `affy` package.
#'
#' @param dirName The path to the directory containing the CEL files.
#' @return An `AffyBatch` object containing the raw data from the CEL files.
#' @importFrom affy ReadAffy
#' @examples
#' # Assuming the CEL files are located in "/Users/briankweiner/R_code/TestData/cel_files/":
#' # affyData <- readAffyCELFiles("/Users/briankweiner/R_code/TestData/cel_files/")
#' @export
readAffyCELFiles <- function(dirName) {
  if (!requireNamespace("affy", quietly = TRUE)) {
    stop("The 'affy' package is required but not installed.")
  }

  # Set the working directory to the directory containing CEL files
  oldDir <- setwd(dirName)
  on.exit(setwd(oldDir)) # Ensure the original working directory is restored

  # Read the CEL files and create an AffyBatch object
  affyData <- affy::ReadAffy()

  return(affyData)
}
