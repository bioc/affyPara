\name{normalizeAffyBatchQuantilesPara}
\alias{normalizeAffyBatchQuantilesPara}
\alias{normalizeQuantilesPara}
\title{Parallelized quantile normalization}
\description{
Parallelized normalization of arrays based upon quantiles.
}
\usage{
normalizeAffyBatchQuantilesPara(object,
    phenoData = new("AnnotatedDataFrame"), cdfname = NULL,
    type = c("separate", "pmonly", "mmonly", "together"), 
    cluster, verbose = FALSE)

normalizeQuantilesPara(cluster, type, object.length, verbose=FALSE)
}

\arguments{
   \item{object}{ An object of class \link[affy:AffyBatch-class]{AffyBatch} 
   		OR a \code{character} vector with the names of CEL files 
   		OR a (partitioned) list of \code{character} vectors with CEL file names.}
 
  \item{phenoData}{ An \link[Biobase]{AnnotatedDataFrame} object. }
  \item{cdfname}{ Used to specify the name of an alternative cdf package. If set to \code{NULL}, 
  		the usual cdf package based on Affymetrix' mappings will be used. }
  \item{type}{A string specifying how the normalization should be applied.}
  \item{cluster}{ A cluster object obtained from the function \link[snow:snow-startstop]{makeCluster} in the SNOW package. 
  		For default \code{.affyParaInternalEnv$cl}  will be used. }
  \item{verbose}{ A logical value. If \code{TRUE} it writes out some messages.} 
  \item{object.length}{ Number of samples, which should be normalized. }
}
\details{
Parallelized normalization of arrays based upon quantiles. 
This method is based upon the concept of a quantile-quantile plot extended to \code{n} dimensions. 
No special allowances are made for outliers.

For the serial function and more details see the function \code{normalize.AffyBatch.quantiles}.

For using this function a computer cluster using the SNOW package has to be started. 
Starting the cluster with the command \code{makeCluster} generates an cluster object in the affyPara environment (.affyParaInternalEnv) and 
no cluster object in the global environment. The cluster object in the affyPara environment will be used as default cluster object, 
therefore no more cluster object handling is required.  
The \code{makeXXXcluster} functions from the package SNOW can be used to create an cluster object in the global environment and 
to use it for the preprocessing functions.

\code{normalizeQuantilesPara} is a internal function which will be executed at all slaves.
\describe{
  \item{\code{normalizeQuantilesPara}}{Function for quantile normalization.}
}
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

  AffyBatch <- normalizeAffyBatchQuantilesPara(Dilution, verbose=TRUE)

  stopCluster()
}
}
}

\keyword{programming}
\keyword{manip}
