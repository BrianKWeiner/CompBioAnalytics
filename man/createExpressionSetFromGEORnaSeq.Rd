% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createExpressionSetFromGEORnaSeq.R
\name{createExpressionSetFromGEORnaSeq}
\alias{createExpressionSetFromGEORnaSeq}
\title{Construct ExpressionSet from GEO SOFT and Count Data}
\usage{
createExpressionSetFromGEORnaSeq(softFilePath, countDataPath)
}
\arguments{
\item{countDataPath}{The path to the GSEXXXXXX_count.txt.gz file containing sample counts.}

\item{softGSEID}{The ID for a GSE file from NCBI's GEO.}
}
\value{
An ExpressionSet object populated with metadata and expression data.
}
\description{
Parses a SOFT formatted family file and a corresponding count data file to
construct an ExpressionSet object with metadata and expression data.
}
\examples{
eset <- createExpressionSetFromGEORnaSeq("GSE145513", "/Users/briankweiner/R_code/TestData/geo_files/GSE145513_count.txt.gz")

}
