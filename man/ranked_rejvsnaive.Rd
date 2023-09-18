\name{ranked_rejvsnaive}
\alias{ranked_rejvsnaive}
\docType{data}
\title{
Rej vs Naive ranked DEGs}
\description{
data frame containing the significantly changed DEGs between rejection and naive, for use as example data set}
\usage{data("ranked_rejvsnaive")}
\format{
  A data frame with 6204 observations on the following 6 variables.
  \describe{
    \item{\code{MGI_Symbol}}{a character vector}
    \item{\code{MGI_Desc}}{a character vector}
    \item{\code{Log2FoldChange}}{a numeric vector}
    \item{\code{Adj_P_Value}}{a numeric vector}
    \item{\code{ChangeDirection}}{a character vector}
    \item{\code{ChangeMetric}}{a numeric vector}
  }
}
\source{
Alegre lab}

\examples{
data(ranked_rejvsnaive)
## maybe str(ranked_rejvsnaive) ; plot(ranked_rejvsnaive) ...
}
\keyword{datasets}
