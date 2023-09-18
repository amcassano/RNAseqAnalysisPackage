\name{anno_rejvsnaive}
\alias{anno_rejvsnaive}
\docType{data}
\title{
annotated rejecting vs naive degs}
\description{
DEGs between naive and rejecting conditions with gene annotations}
\usage{data("anno_rejvsnaive")}
\format{
  A data frame with 55414 observations on the following 8 variables.
  \describe{
    \item{\code{MGI_Symbol}}{a character vector}
    \item{\code{MGI_Desc}}{a character vector}
    \item{\code{baseMean}}{a numeric vector}
    \item{\code{Log2FoldChange}}{a numeric vector}
    \item{\code{Log2FC_SE}}{a numeric vector}
    \item{\code{stat}}{a numeric vector}
    \item{\code{P_Value}}{a numeric vector}
    \item{\code{Adj_P_Value}}{a numeric vector}
  }
}

\source{
Alegre lab}
\examples{
data(anno_rejvsnaive)
## maybe str(anno_rejvsnaive) ; plot(anno_rejvsnaive) ...
}
\keyword{datasets}
