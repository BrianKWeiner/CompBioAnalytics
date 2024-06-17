#' Create DGEList from featureCounts Files
#'
#' This function reads multiple featureCounts output files from a specified directory,
#' extracts essential columns (Gene ID and counts), and assembles them into a DGEList object
#' for use with the edgeR package.
#'
#' @param directoryPath A character string specifying the path to the directory containing
#' featureCounts output files.
#' @return A DGEList object suitable for downstream differential gene expression analysis.
#' @importFrom readr read_tsv
#' @importFrom edgeR DGEList
#' @examples
#' # Assuming your featureCounts output files are in the specified directory:
#' dgelist <- createDGEListFromFeatureCounts("/Users/briankweiner/R_code/TestData/featurecount_files")
#' @export
createDGEListFromFeatureCounts <- function(directoryPath) {
  require(readr)
  require(edgeR)

  #For testing:
  #directoryPath <- "/Users/briankweiner/R_code/TestData/featurecount_files"

  # List all featureCounts.txt files in the directory
  filePaths <- list.files(directoryPath, pattern = "\\.counts\\.txt$", full.names = TRUE)

  # Function to read and process a single file
  processFile <- function(filePath) {
    # Read the file, skipping the first line and only keeping the first and last column
    data <- read_tsv(filePath, skip = 1, col_select = c(Geneid, ends_with("out.bam")), col_types = cols(.default = "c"))

    # Rename the last column to the basename of the file (sample name)
    colnames(data)[2] <- tools::file_path_sans_ext(basename(filePath))
    data
  }

  # Apply the function to each file and combine the results
  countDataList <- lapply(filePaths, processFile)
  combinedCounts <- Reduce(function(x, y) merge(x, y, by = "Geneid"), countDataList)

  # Convert to matrix and create DGEList
  # First convert all columns except the first (gene IDs) to numeric
  combinedCounts[, -1] <- sapply(combinedCounts[, -1], as.numeric)
  countMatrix <- as.matrix(combinedCounts[-1])
  rownames(countMatrix) <- combinedCounts$Geneid

  # Check if any count data cannot be converted to numeric
  # Excludes the first column which contains the Gene IDs
  if(any(sapply(combinedCounts[, -1], function(x) any(is.na(as.numeric(x)))))) {
    stop("Some count data could not be converted to numeric. Please check input files.")
  }


  # Create DGEList
  dgeList <- DGEList(counts = countMatrix)

  return(dgeList)
}
