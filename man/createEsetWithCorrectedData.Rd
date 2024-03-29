% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createEsetWithCorrectedData.R
\name{createEsetWithCorrectedData}
\alias{createEsetWithCorrectedData}
\title{Combine Corrected Data with Existing ExpressionSet}
\usage{
createEsetWithCorrectedData(correctedData, affyBatch, existingEset)
}
\arguments{
\item{correctedData}{A matrix of batch-corrected expression data.}

\item{affyBatch}{An \code{AffyBatch} object from which initial metadata is derived.}

\item{existingEset}{An existing \code{ExpressionSet} object to use for updating metadata in the new \code{ExpressionSet}.}
}
\value{
A new \code{ExpressionSet} object containing the corrected expression data along with updated metadata.
}
\description{
This function creates a new \code{ExpressionSet} object with batch-corrected expression data,
while retaining and updating the metadata (experimentData, phenoData, featureData, annotation,
and protocolData) from an existing \code{ExpressionSet} object.
}
\examples{
# Assuming `correctedData` is your matrix of batch-corrected expression data,
# `affyBatch` is your original AffyBatch object,
# and `existingEset` is an existing ExpressionSet object:
# newEset <- createEsetWithCorrectedData(correctedData, affyBatch, existingEset)
}
