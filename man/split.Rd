\name{splitObjects}
\alias{splitAffyBatch}
\alias{splitFileVector}
\alias{splitMatrix}
\title{Functions to split objects into parts}

\description{
Functions to split an \link[affy:AffyBatch-class]{AffyBatch}, a \code{list} of files and a \code{matrix} into several objects for distributed computing. If possible objects will be of the same size.
}

\usage{
splitAffyBatch(abatch, number.part)
splitFileVector(fileVec, number.part)
splitMatrix(matrix, number.part)
}

\arguments{
  \item{abatch}{ An object of class \link[affy:AffyBatch-class]{AffyBatch}. }
  \item{fileVec}{ A \code{character} vector containing the names of the files. }
  \item{matrix}{ An object of class \link[base]{matrix}. }
  \item{number.part}{ Number of parts to split the object }
}
\details{
\describe{
  \item{\code{splitAffyBatch}}{Splits an \link[affy:AffyBatch-class]{AffyBatch} into a list of AffyBatches. }
  \item{\code{splitFileVector}}{Splits a \code{character} vector of file names into a list of \code{character} vectors with file names. }
  \item{\code{splitMatrix}}{Splits a \code{matrix} by columns into a \code{list} of matrices.}
}

These functions use the functions \link[snow:snow-internal]{splitIndices} and \link[snow:snow-internal]{splitCols} from the \code{SNOW} package.
}
\value{
A \code{list} of the split objects.
}

\author{ Markus Schmidberger \email{schmidb@ibe.med.uni-muenchen.de}, Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de} }

\examples{
library(affyPara)
if (require(affydata)) {
  data(Dilution)

  spAffyB <- splitAffyBatch(Dilution, 2)
}
}

\keyword{programming}
