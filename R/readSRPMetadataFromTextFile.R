#' Read Metadata from Text File
#'
#' Parses a metadata text file into a DataFrame. The first line of the file, which is usually
#' a comment, is skipped. The second line contains column headers. The DataFrame is indexed
#' by the "Run" column.
#'
#' @param filePath The path to the text file containing the metadata.
#' @return A DataFrame with metadata, indexed by the "Run" column.
#' @importFrom readr read_delim
#' @examples
#' metadata <- readMetadataFromTextFile("/Users/briankweiner/R_code/TestData/featurecount_files/SRP501793_SraRunTable.txt")
#' @export
readSRPMetadataFromTextFile <- function(filePath) {

  # For testing:
  #filePath <- "/Users/briankweiner/R_code/TestData/featurecount_files/SRP501793_SraRunTable.txt"
  # Read the file, skipping the first line and using the second line as header
  metadata <- readr::read_delim(filePath, delim = ",", col_names = TRUE, show_col_types = FALSE)

  # Set the row names to the "Run" column
  rownames(metadata) <- metadata$Run

  return(metadata)
}
