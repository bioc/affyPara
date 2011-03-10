\name{normalizeAffyBatchConstantPara}
\alias{normalizeAffyBatchConstantPara}
\title{Parallelized scaling normalization }

\description{
Parallelized scaling normalization of arrays.
}
\usage{
normalizeAffyBatchConstantPara(object, cluster,
      refindex = 1, FUN = mean, na.rm = TRUE,
      phenoData = new("AnnotatedDataFrame"), cdfname = NULL,
      verbose = getOption("verbose"))  
}

\arguments{
   \item{object}{ An object of class \link[affy:AffyBatch-class]{AffyBatch} 
   		OR a \code{character} vector with the names of CEL files 
   		OR a (partitioned) list of \code{character} vectors with CEL file names.}
  
  \item{refindex}{ The index of the array used as reference. }
  \item{FUN}{ A function generating a value from the intensities on an array. 
  		Typically \code{mean} or \code{median}. }
  \item{na.rm}{ Paramater passed to the function FUN. A logical value indicating whether \code{NA} values 
  		should be stripped before the computation proceeds. }
  \item{phenoData}{ An \link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame} object. }
  \item{cdfname}{ Used to specify the name of an alternative cdf package. If set to \code{NULL}, 
  		the usual cdf package based on Affymetrix' mappings will be used. }
  \item{cluster}{ A cluster object obtained from the function \link[snow:snow-startstop]{makeCluster} in the SNOW package. 
  		For default \code{.affyParaInternalEnv$cl}  will be used. }
  \item{verbose}{ A logical value. If \code{TRUE} it writes out some messages. default: getOption("verbose") } 
  
}

\details{
Parallelized scaling normalization of arrays. This means that all the array are scaled so that they have the same mean value.

For the serial function and more details see the function \code{normalize.constant}.

For using this function a computer cluster using the SNOW package has to be started. 
Starting the cluster with the command \code{makeCluster} generates an cluster object in the affyPara environment (.affyParaInternalEnv) and 
no cluster object in the global environment. The cluster object in the affyPara environment will be used as default cluster object, 
therefore no more cluster object handling is required.  
The \code{makeXXXcluster} functions from the package SNOW can be used to create an cluster object in the global environment and 
to use it for the preprocessing functions.
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

  AffyBatch <- normalizeAffyBatchConstantPara(Dilution, verbose=TRUE)

  stopCluster()
}
}
}

\keyword{programming}
\keyword{manip}
