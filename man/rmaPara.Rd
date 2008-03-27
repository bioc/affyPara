\name{rmaPara}
\alias{rmaPara}
\title{ Parallelized PMA preprocessing }
\description{Parallelized preprocessing function, which converts an \link[affy:AffyBatch-class]{AffyBatch} into an  \link[Biobase]{ExpressionSet} using the robust multi-array average (RMA) expression measure.}

\usage{
rmaPara(cluster,
	object,
	ids = NULL,
	phenoData = new("AnnotatedDataFrame"), cdfname = NULL,
	verbose=FALSE)
}

\arguments{
  \item{cluster}{ A cluster object obtained from the function \link[snow:snow-startstop]{makeCluster} in the \code{SNOW} package. }
  \item{object}{ An object of class \link[affy:AffyBatch-class]{AffyBatch} OR a \code{character} vector with the names of CEL files OR a (partitioned) list of \code{character} vectors with CEL file names.}

  \item{ids}{ List of \code{ids} for summarization }

  \item{phenoData}{ An \link[Biobase]{AnnotatedDataFrame} object. }
  \item{cdfname}{ Used to specify the name of an alternative cdf package. If set to \code{NULL}, the usual cdf package based on Affymetrix' mappings will be used. }
  \item{verbose}{ A logical value. If \code{TRUE} it writes out some messages. }
}

\details{
Parallelized preprocessing function, which goes from raw probe intensities to expression values using the robust multi-array average (RMA) expression measure:
 Background correction: rma; Normalization: quantile; Summarization: medianpolish

For the serial function and more details see the function \code{rma}.

For using this function a computer cluster using the \code{snow} package has to be started.

This is a wrapper function for \code{preproPara}.
}

\value{
An object of class \link[Biobase]{ExpressionSet}.
}

\author{ Markus Schmidberger \email{schmidb@ibe.med.uni-muenchen.de}, Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de} }

\examples{
\dontrun{
library(affyPara)

data(affybatch.example)

c1 <- makeCluster(3)

esset <- rmaPara(cluster, affybatch.example)
    
stopCluster(c1)
}
}

\keyword{programming}
\keyword{manip}
