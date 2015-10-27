\name{ftrs}
\alias{ftrs}
\docType{data}
\title{
Well profiles extracted from images using imageHTS
}
\description{
This file contains the \code{data.frame} output of the \code{imageHTS}
pipeline.
}
\value{
Feature summaries calculated over all features extracted from the cells in each
well. The wells are defined by their \code{uname}.
}
\usage{data(ftrs)}
\format{
  A data frame with 36864 observations on 396 variables.
  The first is '\code{uname}', a character vector of the well unames.
  The remaining 395 variables are numeric vectors of feature summaries 
  for each well.
}
\examples{
data(ftrs)
names(ftrs)
}
\seealso{
  \code{\link{extractFeatures}}, \code{\link{computeFeatures}},
  \code{\link{getFeaturesAllSpots}}, \code{\link{summarizeWellsExtended}}
}
\keyword{datasets}
