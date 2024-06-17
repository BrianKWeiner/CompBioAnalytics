#' Construct ExpressionSet from GEO SOFT and Count Data
#'
#' Parses a SOFT formatted family file and a corresponding count data file to
#' construct an ExpressionSet object with metadata and expression data.
#'
#' @param softGSEID The ID for a GSE file from NCBI's GEO.
#' @param countDataPath The path to the GSEXXXXXX_count.txt.gz file containing sample counts.
#' @return An ExpressionSet object populated with metadata and expression data.
#' @importFrom GEOquery getGEOfile
#' @importFrom biomaRt useMart getBM
#' @importFrom Biobase ExpressionSet
#' @importFrom data.table fread
#' @examples
#' eset <- createExpressionSetFromGEORnaSeq("GSE145513", "/Users/briankweiner/R_code/TestData/geo_files/GSE145513_count.txt.gz")
#'
#' @export
createExpressionSetFromGEORnaSeq <- function(softGSEID, countDataPath) {
  if (!requireNamespace("Biobase", quietly = TRUE)) stop("Package 'Biobase' is needed but not installed.")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' is needed but not installed.")

  # For testing:
  #softGSEID <- "GSE145513"
  #countDataPath <- "/Users/briankweiner/R_code/TestData/geo_files/GSE145513_count.txt.gz"

  geoFile = GEOquery::getGEO(softGSEID, GSEMatrix = TRUE, getGPL = FALSE)

  # Step 1: Read and parse the SOFT file for metadata
  processedMetadata <- pData(phenoData(geoFile[[1]]))
  processedMetadata$GSM_ID <- processedMetadata$geo_accession


  # Step 2: Load count data
  exprData <- data.table::fread(countDataPath)

  # Assuming the first column contains gene identifiers and subsequent columns contain counts
  rownames <- exprData[[1]]
  exprData <- as.matrix(exprData[,-1, with = FALSE])
  rownames(exprData) <- rownames

  updateExprDataColnames <- function(exprData, descriptions_df) {
    # Create a named vector to map descriptions to GSM IDs
    # Ensure that 'Description' is unique for each 'GSM_ID'; otherwise, handle duplicates as needed
    description_to_gsm <- setNames(descriptions_df$GSM_ID, descriptions_df$description)

    # Match the column names of exprData with the keys (descriptions) of the mapping vector
    matched_gsm_ids <- description_to_gsm[ colnames(exprData) ]

    # Update colnames of exprData with matched GSM IDs where available
    # This assumes every column name in exprData has a direct match in descriptions_df$Description
    # If not all descriptions are matched, handle unmatched cases as needed
    new_colnames <- ifelse(is.na(matched_gsm_ids), colnames(exprData), matched_gsm_ids)

    colnames(exprData) <- new_colnames


    return(exprData)
  }

  newExprData <- updateExprDataColnames(exprData, processedMetadata)
  newExprData <- as.data.frame(newExprData)

  fetchGeneInfo <- function(gene_ids) {
    # Connect to the Ensembl BioMart database
    ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")

    # Specify the attributes (information to retrieve) for the query
    attributes <- c("ensembl_gene_id", "external_gene_name", "entrezgene_id", "gene_biotype")

    # Query BioMart to retrieve information for the specified gene IDs
    gene_info <- biomaRt::getBM(attributes = attributes,
                                filters = "ensembl_gene_id",
                                values = gene_ids,
                                mart = ensembl)

    # Some post-query processing if necessary, e.g., removing duplicates
    # gene_info <- unique(gene_info)

    return(gene_info)
  }

  # Example usage
  gene_ids <- rownames(exprData)  # Assuming exprData is available
  gene_info_df <- fetchGeneInfo(gene_ids)

  # Check the first few rows of the fetched data
  head(gene_info_df)

  collapseGeneInfo <- function(gene_info_df) {
    library(dplyr)

    # Group by ensembl_gene_id and summarize each column by concatenating unique values
    gene_info_collapsed <- gene_info_df %>%
      group_by(ensembl_gene_id) %>%
      summarise(across(.cols = c(external_gene_name, entrezgene_id, gene_biotype),
                       .fns = ~paste(unique(.x), collapse = " /// "),
                       .names = "{.col}"))

    return(gene_info_collapsed)
  }

  # Assuming gene_info_df is the dataframe returned by fetchGeneInfo
  gene_info_df_collapsed <- collapseGeneInfo(gene_info_df)

  # Convert row names of exprData to a dataframe
  exprData_genes <- data.frame(ensembl_gene_id = rownames(exprData))

  # Perform a left join to match exprData genes with collapsed gene info
  final_gene_info <- left_join(exprData_genes, gene_info_df_collapsed, by = "ensembl_gene_id")

  cleanNAValues <- function(df) {
    df <- df %>%
      mutate(across(everything(), ~na_if(.x, "<NA>"))) %>% # Replace "<NA>" strings with true NA
      mutate(across(everything(), ~na_if(.x, "NA"))) %>% # Replace "NA" strings with true NA
      replace(is.na(.), "") # Ensure all NA values are blank

    return(df)
  }

  # Apply the function to your final gene information dataframe
  final_gene_info_cleaned <- cleanNAValues(final_gene_info)
  rownames(final_gene_info_cleaned) <- final_gene_info_cleaned$ensembl_gene_id

  newExprDataMatrix <- as.data.frame(newExprData)
  rownames(newExprDataMatrix) <- rownames(newExprData)
  colnames(newExprDataMatrix) <- colnames(newExprData)
  newExprDataMatrix <- as.matrix(newExprDataMatrix)

  # Step 3: Construct ExpressionSet
  # Placeholder: Construct phenodata and featureData based on your parsed metadata
  phenoData <- new("AnnotatedDataFrame", data = processedMetadata)
  featureData <- new("AnnotatedDataFrame", data = final_gene_info_cleaned)

  newExprDataMatrix <- newExprDataMatrix[, rownames(phenoData)]

  colnames(newExprDataMatrix)
  rownames(phenoData)

  eset <- Biobase::ExpressionSet(assayData = newExprDataMatrix,
                                 phenoData = phenoData,
                                 featureData = featureData,
                                 annotation = processedMetadata$platform_id[1])

  return(eset)
}
