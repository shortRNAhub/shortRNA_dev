% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annoprep2.R
\name{gettRNA}
\alias{gettRNA}
\title{Obtain tRNA sequences from the GttRNA database}
\usage{
gettRNA(sp = "mm10", mt = TRUE, addCCA = TRUE)
}
\arguments{
\item{sp}{Species for which tRNA sequences to be obtained. Currently
available for mm10, mm39, hg19, and hg38}

\item{mt}{Whether to also download and include mitochondrial tRNAs from
mitotRNA database. Please see also `getMttRNA`}

\item{addCCA}{Whether to add CCA modifications to sequences or not.}
}
\description{
Obtain tRNA sequences from the GttRNA database
}
