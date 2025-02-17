% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qc.R
\name{baseRatioDF}
\alias{baseRatioDF}
\title{Base content data frame from JSON file}
\usage{
baseRatioDF(json, which = "before")
}
\arguments{
\item{json}{A JSON file with QC information from \code{\link{qcFastq}}}
}
\value{
A data frame with FASTQ file information
}
\description{
Base content data frame from JSON file
}
\examples{
# Input
fq1 <- system.file("extdata", "Fox3_Std_small.fq.gz", package = "Rfastp")
outDir <- tempdir()

# Analysis
qc_res <- qc_fastq(file = fq1, outDir = outDir)
df <- baseRatioDF(json = qc_res$json)

# Output
df
}
\author{
Deepak Tanwar (tanward@ethz.ch)
}
