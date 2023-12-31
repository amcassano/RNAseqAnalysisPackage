% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pairwiseDEGresults.R
\name{pairwiseDEGresults}
\alias{pairwiseDEGresults}
\title{Pairwise Differentitally Expressed Genes comparison results}
\usage{
pairwiseDEGresults(
  numerator,
  denominator,
  deseq_obj,
  contrast_str = "Group",
  tidy_result = TRUE,
  log2FCshrink = FALSE
)
}
\arguments{
\item{numerator}{string, name of the conditon of interest for comparison, cannot contain spaces}

\item{denominator}{string, name of baseline condition for comparison}

\item{deseq_obj}{DESeqdata object}

\item{contrast_str}{string, defaults to "Group", the comparison of interest for pairwise comparison}

\item{tidy_result}{boolean, defaults to TRUE, dictates if the results are returned as a dataframe or obj}

\item{log2FCshrink}{boolean, defaults to false, dictates if ashr lfc shrinking formula is used}
}
\value{
DESeq Results object, or a dataframe
}
\description{
Get the results of a deseq object
}
\examples{
pairwiseDEGresults("Tolerant", "Naive", dseqobj)
pairwiseDEGresults("Rejecting", "Tolerant", dseqobj, log2FCshrink = TRUE)

}
