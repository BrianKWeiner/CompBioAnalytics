#' Read GEO Series Matrix File and Return an ExpressionSet
#'
#' This function reads a GEO series matrix file, extracting gene expression data along with sample metadata,
#' and constructs an ExpressionSet object. The ExpressionSet object contains the expression matrix and
#' phenotypic data for the samples.
#'
#' @param file_path The path to the GEO series matrix file.
#' @return An `ExpressionSet` object containing the expression data and sample metadata.
#' @examples
#' geo_eset <- readGEO("/Users/briankweiner/R_code/TestData/geo_files/GSE53552_series_matrix.txt")
#' Biobase::exprs(geo_eset) # To access the expression data
#' Biobase::pData(geo_eset) # To access the sample metadata
#' @export
#' @importFrom GEOquery getGEO
#' @importFrom Biobase ExpressionSet exprs pData
readGEO <- function(file_path) {
  # Ensure the GEOquery and Biobase packages are available
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("The 'GEOquery' package is required but not installed. Please install 'GEOquery' to use readGEO().")
  }
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("The 'Biobase' package is required but not installed. Please install 'Biobase' to use readGEO().")
  }

  checkFileName <- function(file_path) {
    # Extract the file name from the file path
    file_name <- basename(file_path)

    # Check if the file name ends with 'matrix.txt'
    if (!grepl("matrix\\.txt$", file_name)) {
      stop("The file name does not end with 'matrix.txt'. Please provide a valid file.")
    }

    # If the file name is valid, you can proceed
    warn("File name is valid. Proceeding...\n")
  }

  # Load the GEO series matrix file using GEOquery
  geo_series_eset <- GEOquery::getGEO(filename = file_path)

  return(geo_series_eset)
}



