% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addLabels.R
\name{add_labels}
\alias{add_labels}
\title{Add Labels}
\usage{
add_labels(dat, overlaps = 6)
}
\arguments{
\item{dat}{data frame}

\item{overlaps}{numeric, number of overlaps allowed}
}
\value{
geom label repel object
}
\description{
adds labels to each point on a ggplot
}
\examples{
pca_data <- DESeq2::plotPCA(rlognorm, intgroup = Group, returnData = TRUE)
add_labels(pca_data)
}
