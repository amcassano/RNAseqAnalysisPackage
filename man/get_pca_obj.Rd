% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getPCAobj.R
\name{get_pca_obj}
\alias{get_pca_obj}
\title{Get PCA object}
\usage{
get_pca_obj(norm_df, metadat, varianceRemove = 0.9)
}
\arguments{
\item{norm_df}{dataframe of normalized counts (rlog df)}

\item{metadat}{metadata dataframe}

\item{varianceRemove}{Remove this percentage of variables based on low variance.}
}
\value{
pca_obj
}
\description{
Creates the PCA object that is used in the PCA plot as well as other things like the Screeplot
}
\examples{
get_pca_obj(rlog_df, meta, 0.1)
}
