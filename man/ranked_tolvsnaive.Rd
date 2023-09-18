\name{ranked_tolvsnaive}
\alias{ranked_tolvsnaive}
\docType{data}
\title{
Tol vs Naive ranked DEGs
}
\description{
data frame containing the significantly changed DEGs between tolerance and naive, for use as example data set

}
\usage{data("ranked_tolvsnaive")}
\format{
  A data frame with 4274 observations on the following 6 variables.
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
Alegre Lab}

\examples{
data(ranked_tolvsnaive)
## maybe str(ranked_tolvsnaive) ; plot(ranked_tolvsnaive) ...
}
\keyword{datasets}
