\name{justvsnPara}
\alias{justvsnPara}
\title{Parallel Fit the vsn model}
\description{
'justvsnPara' fits the vsn model to intensity data in an AffyBatch. It has the same functionallity than justvsn but is implemented in parallel (and only supports an AffyBatch as input data).
}
\usage{
justvsnPara(cluster,
	object,
	phenoData = new("AnnotatedDataFrame"), cdfname = NULL,
	subsample, ...,
	verbose=TRUE) 
}

\arguments{
  \item{cluster}{ A cluster object obtained from the function \link[snow:snow-startstop]{makeCluster} in the \code{SNOW} package. }
  \item{object}{ An object of class \link[affy:AffyBatch-class]{AffyBatch} OR a \code{character} vector with the names of CEL files OR a (partitioned) list of \code{character} vectors with CEL file names.}
 
  \item{phenoData}{ An \link[Biobase]{AnnotatedDataFrame} object. }
  \item{cdfname}{ Used to specify the name of an alternative cdf package. If set to \code{NULL}, the usual cdf package based on Affymetrix' mappings will be used. }
 
   \item{subsample}{Integer of length 1. If specified, the model parameters are
    estimated from a subsample of the data of size \code{subsample}
    only, yet the fitted transformation is
    then applied to all data. For large datasets, this can substantially
    reduce the CPU time and memory consumption at a negligible loss of precision.}
  \item{\dots}{Further arguments that get passed and are similar to \code{vsn2}.}
  \item{verbose}{ A logical value. If \code{TRUE} it writes out some messages. } 
}
\details{
For the serial function and more details see the function \code{vsn2}.

For using this function a computer cluster using the \code{snow} package has to be started.
}

\value{
An \link[affy:AffyBatch-class]{AffyBatch} of normalized objects.
}

\author{ Markus Schmidberger \email{schmidb@ibe.med.uni-muenchen.de}, Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de} }

\examples{
\dontrun{
library(affyPara)

data(affybatch.example)

c1 <- makeCluster(3)

AB1 <- justvsnPara(cluster, affybatch.example, verbose=verbose )

stopCluster(c1)
}
}

\keyword{programming}
\keyword{manip}