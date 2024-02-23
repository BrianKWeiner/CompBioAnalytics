#' Read Variants from a VCF File
#'
#' This function reads variant call format (VCF) files and extracts the variant information into a data frame.
#' The resulting data frame includes columns for the chromosome, position, ID, reference allele, alternate allele(s),
#' quality, filter, and information (INFO) fields. Additional columns for genotype information may be included if present.
#'
#' @param file_path The path to the VCF file.
#' @param include_genotypes A logical value indicating whether to include genotype information for each sample. Defaults to FALSE.
#' @return A data frame with variant information, including columns for chromosome, position, ID, reference allele,
#'         alternate allele(s), quality, filter, and INFO fields. If `include_genotypes` is TRUE, additional columns
#'         for genotype information will be included.
#' @examples
#' vcf_data <- readVCF("/Users/briankweiner/R_code/TestData/vcf_files/GCA_000001215.4_current_ids.sample.vcf)
#' head(vcf_data)
#' @export
#' @importFrom VariantAnnotation readVcf
#' @importFrom VariantAnnotation scanVcfHeader
readVCF <- function(file_path, include_genotypes = FALSE) {
  # Ensure the VariantAnnotation package is available for working with VCF files
  if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
    stop("The 'VariantAnnotation' package is required but not installed. Please install 'VariantAnnotation' to use readVCF().")
  }

  # Read the VCF file header to get sample names if genotypes are to be included
  vcf_header <- VariantAnnotation::scanVcfHeader(file_path)
  sample_names <- colnames(VariantAnnotation::geno(vcf_header))

  # Read the VCF file
  vcf <- VariantAnnotation::readVcf(file_path, genome = "hg19")

  # Extract key variant information into a data frame
  variants_df <- data.frame(
    chromosome = as.character(GenomeInfoDb::seqnames(vcf@rowRanges)),
    position = start(vcf@rowRanges),
    id = names(vcf@rowRanges),
    ref = as.character(vcf@fixed$REF),
    alt = paste(sapply(vcf@fixed$ALT, function(x) paste(x, collapse = ","))),
    qual = vcf@fixed$QUAL,
    filter = sapply(vcf@fixed$FILTER, function(x) paste(x, collapse = ";"))
  )

  # Include genotype information if requested
  if (include_genotypes) {
    geno_data <- do.call(cbind, lapply(sample_names, function(sample) {
      vcf@colData$Samples
    }))
    variants_df <- cbind(variants_df, geno_data)
  }

  return(variants_df)
}
