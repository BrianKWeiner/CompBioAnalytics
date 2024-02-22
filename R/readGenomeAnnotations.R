#' Read Genome Annotations from GFF File
#'
#' This function reads genomic annotations from a GFF (General Feature Format) file and returns a data frame with the annotations.
#' The data frame will include columns for sequence name, source, feature type, start position, end position, score, strand, phase, and attributes.
#'
#' @param file_path The path to the GFF file.
#' @param features Optional vector of feature types to filter (e.g., "gene", "exon"). If NULL, all features are returned.
#' @return A data frame with columns for sequence, source, feature, start, end, score, strand, phase, and attributes.
#' @examples
#' gff_annotations <- readGenomeAnnotations("/Users/briankweiner/R_code/TestData/gff_files/GCF_000005845.2_ASM584v2_genomic.gff")
#' head(gff_annotations)
#' @export
#' @importFrom readr read_delim
readGenomeAnnotations <- function(file_path, features = NULL) {
  # Ensure the readr package is available for reading the file
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("The 'readr' package is required but not installed. Please install 'readr' to use readGenomeAnnotations().")
  }

  # Read the GFF file
  gff_data <- readr::read_delim(file_path,
                                delim = "\t",
                                comment = "#",
                                col_names = c("seqname",
                                              "source",
                                              "feature",
                                              "start",
                                              "end",
                                              "score",
                                              "strand",
                                              "phase",
                                              "attributes"),
                                na = ".",
                                quote = "",
                                col_types = readr::cols(.default = readr::col_character()))

  # Filter by feature type if specified
  if (!is.null(features)) {
    gff_data <- gff_data[gff_data$feature %in% features, ]
  }

  return(gff_data)
}
