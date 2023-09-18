\name{palette}
\alias{palette}
\docType{data}
\title{
Aes palette}
\description{
example of an aesthetics palette for use in plotting}
\usage{data("palette")}
\format{
  A data frame with 3 observations on the following 7 variables.
  \describe{
    \item{\code{labels}}{a factor with levels \code{Naive} \code{Tolerant} \code{Rejecting}}
    \item{\code{shapeDescription}}{a character vector}
    \item{\code{shapeID}}{a numeric vector}
    \item{\code{outlineName}}{a character vector}
    \item{\code{outlineID}}{a character vector}
    \item{\code{fillName}}{a character vector}
    \item{\code{fillID}}{a character vector}
  }
}

\examples{
data(palette)
## maybe str(palette) ; plot(palette) ...
}
\keyword{datasets}
