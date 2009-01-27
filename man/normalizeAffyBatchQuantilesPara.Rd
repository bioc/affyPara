\name{normalizeAffyBatchQuantilesPara}
\alias{normalizeAffyBatchQuantilesPara}
\alias{normalizeQuantilesPara}
\alias{normalizeQuantilesParaSF1}
\alias{normalizeQuantilesParaSF2}
\title{Parallelized quantile normalization}
\description{
Parallelized normalization of arrays based upon quantiles.
}
\usage{
normalizeAffyBatchQuantilesPara(cluster,
    object,
    phenoData = new("AnnotatedDataFrame"), cdfname = NULL,
    type = c("separate", "pmonly", "mmonly", "together"), verbose = FALSE)

normalizeQuantilesPara(cluster, type, object.length, verbose=FALSE)
normalizeQuantilesParaSF1(type)
normalizeQuantilesParaSF2(row_mean)
}

\arguments{
  \item{cluster}{ A cluster object obtained from the function \link[snow:snow-startstop]{makeCluster} in the \code{SNOW} package. }
  \item{object}{ An object of class \link[affy:AffyBatch-class]{AffyBatch} OR a \code{character} vector with the names of CEL files OR a (partitioned) list of \code{character} vectors with CEL file names.}
 
  \item{phenoData}{ An \link[Biobase]{AnnotatedDataFrame} object. }
  \item{cdfname}{ Used to specify the name of an alternative cdf package. If set to \code{NULL}, the usual cdf package based on Affymetrix' mappings will be used. }
  \item{type}{A string specifying how the normalization should be applied.}
  \item{verbose}{ A logical value. If \code{TRUE} it writes out some messages. } 
  \item{object.length}{ Number of samples, which should be normalized. }
  \item{row_mean}{ Row mean used for normalization. }
}
\details{
Parallelized normalization of arrays based upon quantiles. This method is based upon the concept of a quantile-quantile plot extended to \code{n} dimensions. No special allowances are made for outliers.

For the serial function and more details see the function \code{normalize.AffyBatch.quantiles}.

For using this function a computer cluster using the \code{snow} package has to be started.

\code{normalizeQuantilesPara}, \code{normalizeQuantilesParaSF1} and \code{normalizeQuantilesParaSF2} are internal function which will be executed at all slaves.
\item{\code{normalizeQuantilesPara}}{Function for quantil normalization.}
\item{\code{normalizeQuantilesParaSF1}}{Slavefuntion to calculate row means.}
\item{\code{normalizeQuantilesParaSF2}}{Slavefunction to do normalization.}
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

  c1 <- makeCluster(3)

  AffyBatch <- normalizeAffyBatchQuantilesPara(c1, Dilution, verbose=TRUE)

  stopCluster(c1)
}
}
}

\keyword{programming}
\keyword{manip}
