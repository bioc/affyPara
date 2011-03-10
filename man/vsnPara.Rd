\name{vsnPara}
\alias{vsn2Para}
\alias{vsnrmaPara}
\alias{justvsnPara}
\title{Parallel fir of the vsn model}
\description{
These parallel functions fit the vsn model to intensity data in an AffyBatch. 
They hav the same functionality than the vsn methods in the \code{vsn} package but are implemented in parallel (and only supports an AffyBatch as input data).
}
\usage{
vsn2Para(object, cluster,
	phenoData = new("AnnotatedDataFrame"), cdfname = NULL,
	reference, subsample,
	..., 
	verbose = getOption("verbose")) 

justvsnPara(object, cluster,
		...,
		verbose = getOption("verbose")) 
	
vsnrmaPara(object, cluster,
	pmcorrect.method="pmonly", pmcorrect.param=list(),
	summary.method="medianpolish", summary.param=list(),
	ids=NULL,
	phenoData = new("AnnotatedDataFrame"), cdfname = NULL,
	..., 
	verbose = getOption("verbose"))  
}

\arguments{
  \item{object}{ An object of class \link[affy:AffyBatch-class]{AffyBatch} 
  		OR a \code{character} vector with the names of CEL files 
  		OR a (partitioned) list of \code{character} vectors with CEL file names.}
 
  \item{phenoData}{ An \link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame} object. }
  \item{cdfname}{ Used to specify the name of an alternative cdf package. If set to \code{NULL}, 
  		the usual cdf package based on Affymetrix' mappings will be used. }
 
  \item{subsample}{Integer of length 1. If specified, the model parameters are
    estimated from a subsample of the data of size \code{subsample}
    only, yet the fitted transformation is
    then applied to all data. For large datasets, this can substantially
    reduce the CPU time and memory consumption at a negligible loss of precision.}
  \item{reference}{Optional, a \code{\linkS4class{vsn}} object from
    a previous fit. If this argument is specified, the data
    are normalized "towards" an existing set of reference arrays whose
    parameters are stored in the object \code{reference}. If this
    argument is not specified, then the data are normalized
    "among themselves".}
  \item{\dots}{Further arguments that get passed and are similar to \code{vsn2}.}
  \item{cluster}{ A cluster object obtained from the function \link[snow:snow-startstop]{makeCluster} in the SNOW package. 
  		For default \code{.affyParaInternalEnv$cl} will be used. }
  \item{verbose}{ A logical value. If \code{TRUE} it writes out some messages. default: getOption("verbose") } 
 
  \item{pmcorrect.method}{ The name of the PM adjustement method. }
  \item{pmcorrect.param}{ A list of parameters for \code{pmcorrect.method} (if needed/wanted). }
  \item{summary.method}{ The method used for the computation of expression values }
  \item{summary.param}{ A list of parameters to be passed to the \code{summary.method} (if wanted). }
  \item{ids}{ List of \code{ids} for summarization }
  
}
\details{
For the serial function and more details see the function \code{vsn2}.

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

  AB1 <- justvsnPara(Dilution, verbose=verbose )

  stopCluster()
}
}
}

\keyword{programming}
\keyword{manip}