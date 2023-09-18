\name{anno_tolvsnaive}
\alias{anno_tolvsnaive}
\docType{data}
\title{
annotated data frame}
\description{
example results data frame with gene annotations}
\usage{data("anno_tolvsnaive")}
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

\examples{
data(anno_tolvsnaive)
## maybe str(anno_tolvsnaive) ; plot(anno_tolvsnaive) ...
}
\keyword{datasets}
