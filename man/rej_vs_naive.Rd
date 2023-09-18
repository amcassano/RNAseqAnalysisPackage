\name{rej_vs_naive}
\alias{rej_vs_naive}
\docType{data}
\title{
rej vs naive sample results dataset}
\description{
example results data frame }
\usage{data("rej_vs_naive")}
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

\source{
Alegre}

\examples{
data(rej_vs_naive)
## maybe str(rej_vs_naive) ; plot(rej_vs_naive) ...
}
\keyword{datasets}
