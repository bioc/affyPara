\name{normalizeAffyBatchLoessPara}
\alias{normalizeAffyBatchLoessPara}
\title{Parallelized loess normalization}
\description{
Parallelized loess normalization of arrays.
}
\usage{
normalizeAffyBatchLoessPara(object,
	phenoData = new("AnnotatedDataFrame"), cdfname = NULL,
	type=c("separate","pmonly","mmonly","together"), 
	subset = NULL,
	epsilon = 10^-2, maxit = 1, log.it = TRUE, 
	span = 2/3, family.loess ="symmetric",
	cluster, verbose = getOption("verbose")) 

}

\arguments{
  \item{object}{ An object of class \link[affy:AffyBatch-class]{AffyBatch} 
  		OR a \code{character} vector with the names of CEL files 
  		OR a (partitioned) list of \code{character} vectors with CEL file names.}
 
  \item{phenoData}{ An \link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame} object. }
  \item{cdfname}{ Used to specify the name of an alternative cdf package. If set to \code{NULL}, 
  		the usual cdf package based on Affymetrix' mappings will be used. }
  \item{type}{A string specifying how the normalization should be applied.}
  \item{subset}{a subset of the data to fit a loess to.}
  \item{epsilon}{a tolerance value (supposed to be a small value - used as a stopping criterium).}
  \item{maxit}{maximum number of iterations.}
  \item{log.it}{logical. If \code{TRUE} it takes the log2 of mat}
  \item{span}{parameter to be passed the function loess}
  \item{family.loess}{parameter to be passed the function loess. "gaussian" or "symmetric" are acceptable values for this parameter.}
  \item{cluster}{ A cluster object obtained from the function \link[snow:snow-startstop]{makeCluster} in the SNOW package. 
  		For default \code{.affyParaInternalEnv$cl}  will be used. }
  \item{verbose}{ A logical value. If \code{TRUE} it writes out some messages. default: getOption("verbose") } 
}
\details{
Parallelized loess normalization of arrays.

For the serial function and more details see the function \code{normalize.AffyBatch.loess}.

For using this function a computer cluster using the SNOW package has to be started. 
Starting the cluster with the command \code{makeCluster} generates an cluster object in the affyPara environment (.affyParaInternalEnv) and 
no cluster object in the global environment. The cluster object in the affyPara environment will be used as default cluster object, 
therefore no more cluster object handling is required.  
The \code{makeXXXcluster} functions from the package SNOW can be used to create an cluster object in the global environment and 
to use it for the preprocessing functions.

In the loess normalization the arrays will compared by pairs. Therefore at every node minimum two arrays have to be!
}

\value{
An \link[affy:AffyBatch-class]{AffyBatch} of normalized objects.
}

\author{ Markus Schmidberger \email{schmidb@ibe.med.uni-muenchen.de}, Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de} }

\examples{
\dontrun{
library(affyPara)
if (require(affydata)) {
  data(Dilution)

  makeCluster(3)

  AffyBatch <- normalizeAffyBatchLoessPara(Dilution, verbose=TRUE)

  stopCluster()
}
}
}

\keyword{programming}
\keyword{manip}
