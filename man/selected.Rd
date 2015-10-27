\name{selected}
\alias{selected}
\docType{data}
\title{
Result of the feature selection function
}
\description{
This file contains the result of the feature selection function for 40 
itterations.
}
\value{
A \code{list} containing the \code{vector} \code{selected}
of selected features, the \code{vector} \code{correlation} of their
correlations at the time they were selected, the \code{vector} 
\code{ratioPositive} with the fraction of positive correlations for each
iteration and a \code{list} \code{correlationAll} which contains the 
correlations of all features at each iteration step.
}
\usage{data(selected)}
\format{
\preformatted{
List of 4
 $ selected      : chr [1:40] 
 $ correlation   : num [1:40] 
 $ ratioPositive : num [1:40] 
 $ correlationAll:List of 40 Named num vectors
}
}
\examples{
data(selected)
str(selected)
}
\keyword{datasets}
