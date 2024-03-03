#' Categorize CEL Files into Processing Batches Based on Scan Dates
#'
#' Takes an `AffyBatch` object and categorizes CEL files into processing batches based on the
#' precise date and time they were run, attempting to accommodate variable date formats.
#'
#' @param affyData An `AffyBatch` object.
#' @return A data frame with columns 'FileName', "OriginalDateTime", 'FormattedDateTime',
#' and 'Batch', listing each CEL file, its scan date and time, and its assigned processing
#' batch number.
#' @importFrom Biobase sampleNames
#' @examples
#' # Assuming `affyData` is your AffyBatch object:
#' # batches <- categorizeCelDataIntoBatches(affyData)
#' @export
categorizeCelDataIntoBatches <- function(affyData) {
  requireNamespace("Biobase", quietly = TRUE)

  # Extract the scan dates from the AffyBatch object
  scanDates <- affyData@protocolData@data$ScanDate

  #scanDates <- GSE53552_cel_data@protocolData@data$ScanDate

  # Define possible date formats
  # This may need to be modified in the future if there are more weird date formats
  dateFormats <- c("%m/%d/%y %H:%M:%S",      # "11/04/08 15:33:30"
                   "%Y-%m-%dT%H:%M:%SZ",     # "2009-06-09T23:21:23Z"
                   "%d/%m/%Y %H:%M:%S")      # "04/12/11 15:33:30"

  # Initialize vector to store parsed dates
  parsedDates <- as.POSIXct(rep(NA, length(scanDates)))

  # Track dates that couldn't be parsed
  unparsedDates <- character(0)

  # Attempt to parse each date using as.POSIXct with multiple formats
  for (i in seq_along(scanDates)) {
    for (format in dateFormats) {
      parsedDate <- as.POSIXct(scanDates[i], format = format)
      if (!is.na(parsedDate)) {
        parsedDates[i] <- parsedDate
        break
      }
    }
    # If no format worked, record the unparsed date
    if (is.na(parsedDates[i])) {
      unparsedDates <- c(unparsedDates, scanDates[i])
    }
  }

  # Warn about unparsed dates
  if (length(unparsedDates) > 0) {
    warning("The following date(s) could not be resolved: ", paste(unparsedDates, collapse = ", "),
            "\nPlease modify the dateFormats variable in the categorizeCelDataIntoBatches.R code.")
  }

  # Sort scan dates to ensure chronological order
  orderedDates <- sort(parsedDates, na.last = TRUE)

  # Initialize batches
  batches <- numeric(length(parsedDates))
  batchNumber <- 1

  #Make the the first instance of the sorted dates to be batch 1
  batches[!is.na(orderedDates)][1] <- batchNumber
  batches[which(orderedDates[1] == parsedDates)] <- batchNumber

  # Assign batches based on 24-hour difference
  for (i in 2:length(orderedDates)) {
    if (!is.na(orderedDates[i]) &&
        difftime(orderedDates[i], orderedDates[i-1], units = "hours") <= 24) {
      batches[which(orderedDates[i] == parsedDates)] <- batchNumber
    } else if (!is.na(orderedDates[i])) {
      batchNumber <- batchNumber + 1
      batches[which(orderedDates[i] == parsedDates)] <- batchNumber
    }
  }

  # Create a data frame to return
  batchInfo <- data.frame(
    FileName = sampleNames(affyData),
    OriginalDateTime = affyData@protocolData@data$ScanDate,
    FormattedDateTime = parsedDates,
    Batch = batches,
    stringsAsFactors = FALSE
  )

  return(batchInfo)
}
