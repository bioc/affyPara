\name{mergeSplitObjects}
\alias{mergeAffyBatches}
\alias{combineMatrices}
\title{Merge a list of split objects}

\description{
Functions to merge or combine a \code{list} of split objects (AffyBatch, Matrix).
}

\usage{
mergeAffyBatches(abatch.list, description = NULL, notes = character(0))
combineMatrices(matrix.list, verbose = getOption("verbose"))
}

\arguments{
  \item{abatch.list}{ A \code{list} of objects of class \link[affy:AffyBatch-class]{AffyBatch}. }
  \item{description}{ A \link[Biobase]{MIAME} object. }
  \item{notes}{ A \code{character} vector of explanatory text. }
  \item{matrix.list}{A \code{list} of objects of class \link[base]{matrix}. }
  \item{verbose}{ A logical value. If \code{TRUE} it writes out some messages. }
}
\details{
Functions to merge or combine a \code{list} of split objects.
\describe{
  \item{\code{mergeAffyBatches}}{Merges a \code{list} of AffyBatches to one AffyBatch.}
  \item{\code{combineMatrices}}{Combines a \code{list} of matrices by columns to one matrix.}
}
}

\value{
\item{\code{mergeAffyBatches}}{Returns ONE object of class \link[affy:AffyBatch-class]{AffyBatch}.}
\item{\code{combineMatrices}}{Returns ONE object of class \link[base]{matrix}.}
}

\author{ Markus Schmidberger \email{schmidb@ibe.med.uni-muenchen.de}, Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de} }

\examples{
library(affyPara)
if (require(affydata)) {
  data(Dilution)

  #split AffyBatch
  abatch.list<- splitAffyBatch(Dilution, 2)

  #Merge AffyBatch
  AffyBatch <- mergeAffyBatches(abatch.list)

  # Create matrices
  a <- matrix(1:25, nrow=5)
  b <- matrix(101:125, nrow=5)
  matrix.list <- list(a,b)

  # Combine matrices
  combineMatrices(matrix.list)
}
}

\keyword{programming}
