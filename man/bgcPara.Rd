\name{bgCorrectPara}
\alias{bgCorrectPara}
\alias{bgCorrectParaSF}
\title{Parallelized Background Correction}

\description{
Parallelized functions for background correction of probe intensities.
}
\usage{
bgCorrectPara(cluster,
	object, phenoData = new("AnnotatedDataFrame"), method, 
	verbose = FALSE)

bgCorrectParaSF(method) 
}

\arguments{
  \item{cluster}{ A cluster object obtained from the function \link[snow:snow-startstop]{makeCluster} in the \code{SNOW} package. }
  \item{object}{ An object of class \link[affy:AffyBatch-class]{AffyBatch} OR a \code{character} vector with the names of CEL files OR a (partitioned) list of \code{character} vectors with CEL file names.}
  \item{phenoData}{ An \link[Biobase]{AnnotatedDataFrame} object }
  \item{method}{ A \code{character} that defines what background correction method will be used. Available methods are given by \code{bg.correct.methods}. The name of the method to apply must be double-quoted. }
  \item{verbose}{ A logical value. If \code{TRUE} it writes out some messages. }
}
\details{
\code{bgCorrectPara} is the parallelized function for background correction of probe intensities. For serial function an more details see \link[affy]{bg.correct}.

For using this function a computer cluster using the \code{snow} package has to be started.

\code{bgCorrectParaSF} is a internal function which will be executed at slaves.
\item{\code{bgCorrectParaSF}}{Calls \code{bg.correct} at slaves.}
}
\value{
  An \code{\link[affy:AffyBatch-class]{AffyBatch}} for which the intensities have been background adjusted. For some methods (RMA), only PMs are corrected and the MMs remain the same. 
}

\author{ Markus Schmidberger \email{schmidb@ibe.med.uni-muenchen.de}, Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de} }

\examples{
\dontrun{
library(affyPara)
if (require(affydata)) {
  data(Dilution)

  c1 <- makeCluster(3)

  ##bgc will be the bg corrected version of Dilution 
  bgc <- bgCorrectPara(c1, Dilution, method="rma", verbose=TRUE)

  stopCluster(c1)
}
}
}

\keyword{programming}
\keyword{manip}
