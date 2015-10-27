\name{datamatrixTransformed}
\alias{datamatrixTransformed}
\docType{data}
\title{
Intermediate data of the analysis pipeline
}
\description{
Glog transformed feature data represented as a \code{list} containing
\code{array} \code{D} and \code{list} \code{anno}.
}
\value{
The four dimensional \code{array} \code{D} contains the glog transformed feature
data which were rearanged from the original \code{ftrs} \code{data.frame}. 
The dimensions represent:
  \enumerate{
    \item{drug}
    \item{cell line}
    \item{replicate}
    \item{feature}
  }

The annotation \code{anno} is represented as a \code{list} containing a 
\code{data.frame} \code{drug} with the drug annotation, a \code{data.frame} 
\code{line} with the cell line annotation, a \code{vector} \code{repl} with the
information about the replicates and a \code{vector} \code{ftr} with the
feature names.
}
\usage{data(datamatrixTransformed)}
\format{
\preformatted{
List of 2
 $ D   : num [1:1372, 1:12, 1:2, 1:385]
 $ anno:List of 4
  ..$ drug:'data.frame':	1372 obs. of  15 variables
  ..$ line:'data.frame':	12 obs. of  4 variables
  ..$ repl: int [1:2]
  ..$ ftr : chr [1:385]
}
}
\examples{
data(datamatrixTransformed)
str(datamatrixTransformed)
}
\keyword{datasets}
