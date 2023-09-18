\name{tol_vs_naive}
\alias{tol_vs_naive}
\docType{data}
\title{
tol vs naive}
\description{
results data frame example}
\usage{data("tol_vs_naive")}
\format{
  A data frame with 55414 observations on the following 7 variables.
  \describe{
    \item{\code{GeneID}}{a character vector}
    \item{\code{baseMean}}{a numeric vector}
    \item{\code{log2FoldChange}}{a numeric vector}
    \item{\code{lfcSE}}{a numeric vector}
    \item{\code{stat}}{a numeric vector}
    \item{\code{pvalue}}{a numeric vector}
    \item{\code{padj}}{a numeric vector}
  }
}

\examples{
data(tol_vs_naive)
## maybe str(tol_vs_naive) ; plot(tol_vs_naive) ...
}
\keyword{datasets}
