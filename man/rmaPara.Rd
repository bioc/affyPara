\name{rmaPara}
\alias{rmaPara}
\title{ Parallelized PMA preprocessing }
\description{Parallelized preprocessing function, which converts an \link[affy:AffyBatch-class]{AffyBatch} into an  \link[Biobase]{ExpressionSet} using the robust multi-array average (RMA) expression measure.}

\usage{
rmaPara(object,
	ids = NULL,
	phenoData = new("AnnotatedDataFrame"), cdfname = NULL,
	cluster, verbose = getOption("verbose"))
}

\arguments{
  \item{object}{ An object of class \link[affy:AffyBatch-class]{AffyBatch} 
  		OR a \code{character} vector with the names of CEL files 
  				OR a (partitioned) list of \code{character} vectors with CEL file names.}

  \item{ids}{ List of \code{ids} for summarization }

  \item{phenoData}{ An \link[Biobase]{AnnotatedDataFrame} object. }
  \item{cdfname}{ Used to specify the name of an alternative cdf package. If set to \code{NULL}, 
  		the usual cdf package based on Affymetrix' mappings will be used. }
  \item{cluster}{ A cluster object obtained from the function \link[snow:snow-startstop]{makeCluster} in the SNOW package. 
  		For default \code{.affyParaInternalEnv$cl} will be used. }
  \item{verbose}{ A logical value. If \code{TRUE} it writes out some messages. default: getOption("verbose") } 
}

\details{
Parallelized preprocessing function, which goes from raw probe intensities to expression values using the robust multi-array average (RMA) expression measure:
 Background correction: rma; Normalization: quantile; Summarization: medianpolish

For the serial function and more details see the function \code{rma}.

For using this function a computer cluster using the SNOW package has to be started. 
Starting the cluster with the command \code{makeCluster} generates an cluster object in the affyPara environment (.affyParaInternalEnv) and 
no cluster object in the global environment. The cluster object in the affyPara environment will be used as default cluster object, 
therefore no more cluster object handling is required.  
The \code{makeXXXcluster} functions from the package SNOW can be used to create an cluster object in the global environment and 
to use it for the preprocessing functions.

This is a wrapper function for \code{preproPara}.
}

\value{
An object of class \link[Biobase]{ExpressionSet}.
}

\author{ Markus Schmidberger \email{schmidb@ibe.med.uni-muenchen.de}, Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de} }

\examples{
\dontrun{
library(affyPara)
if (require(affydata)) {
  data(Dilution)

  makeCluster(3)

  esset <- rmaPara(Dilution)

  stopCluster()
}
}
}

\keyword{programming}
\keyword{manip}
