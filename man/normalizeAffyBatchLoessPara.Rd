\name{normalizeAffyBatchLoessPara}
\alias{normalizeAffyBatchLoessPara}
\title{Parallelized loess normalization}
\description{
Parallelized loess normalization of arrays.
}
\usage{
normalizeAffyBatchLoessPara(cluster,
	object,
	phenoData = new("AnnotatedDataFrame"), cdfname = NULL,
	type=c("separate","pmonly","mmonly","together"), 
	subset = NULL,
	epsilon = 10^-2, maxit = 1, log.it = TRUE, 
	span = 2/3, family.loess ="symmetric",
	verbose=FALSE) 

}

\arguments{
  \item{cluster}{ A cluster object obtained from the function \link[snow:snow-startstop]{makeCluster} in the \code{SNOW} package. }
  \item{object}{ An object of class \link[affy:AffyBatch-class]{AffyBatch} OR a \code{character} vector with the names of CEL files OR a (partitioned) list of \code{character} vectors with CEL file names.}
 
  \item{phenoData}{ An \link[Biobase]{AnnotatedDataFrame} object. }
  \item{cdfname}{ Used to specify the name of an alternative cdf package. If set to \code{NULL}, the usual cdf package based on Affymetrix' mappings will be used. }
  \item{type}{A string specifying how the normalization should be applied.}
  \item{subset}{a subset of the data to fit a loess to.}
  \item{epsilon}{a tolerance value (supposed to be a small value - used as a stopping criterium).}
  \item{maxit}{maximum number of iterations.}
  \item{log.it}{logical. If \code{TRUE} it takes the log2 of mat}
  \item{span}{parameter to be passed the function loess}
  \item{family.loess}{parameter to be passed the function loess. "gaussian" or "symmetric" are acceptable values for this parameter.}
  \item{verbose}{ A logical value. If \code{TRUE} it writes out some messages. } 
}
\details{
Parallelized loess normalization of arrays.

For the serial function and more details see the function \code{normalize.AffyBatch.loess}.

For using this function a computer cluster using the \code{snow} package has to be started.
In the loess normalization the arrays will compared by pairs. Therefore at every node minimum two arrays have to be!
}

\value{
An \link[affy:AffyBatch-class]{AffyBatch} of normalized objects.
}

\author{ Markus Schmidberger \email{schmidb@ibe.med.uni-muenchen.de}, Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de} }

\examples{
\dontrun{
library(affyPara)

data(Dilution)

c1 <- makeCluster(3)

AffyBatch <- normalizeAffyBatchLoessPara(c1, Dilution, verbose=TRUE)

stopCluster(c1)
}
}

\keyword{programming}
\keyword{manip}
