\name{normalizeAffyBatchConstantPara}
\alias{normalizeAffyBatchConstantPara}
\alias{normalizeConstantParaSF1}
\alias{normalizeConstantParaSF2}
\title{Parallelized scaling normalization }

\description{
Parallelized scaling normalization of arrays.
}
\usage{
normalizeAffyBatchConstantPara(cluster,
      object,
      refindex = 1, FUN = mean, na.rm = TRUE,
      phenoData = new("AnnotatedDataFrame"), cdfname = NULL,
      verbose = FALSE)  
       
normalizeConstantParaSF1(refindexname)
normalizeConstantParaSF2(refconstant, refindexname, FUN = mean, na.rm = TRUE)
}

\arguments{
  \item{cluster}{ A cluster object obtained from the function \link[snow:snow-startstop]{makeCluster} in the \code{SNOW} package. }
  \item{object}{ An object of class \link[affy:AffyBatch-class]{AffyBatch} OR a \code{character} vector with the names of CEL files OR a (partitioned) list of \code{character} vectors with CEL file names.}
  
  \item{refindex}{ The index of the array used as reference. }
  \item{FUN}{ A function generating a value from the intensities on an array. Typically \code{mean} or \code{median}. }
  \item{na.rm}{ Paramater passed to the function FUN. A logical value indicating whether \code{NA} values should be stripped before the computation proceeds. }
  \item{phenoData}{ An \link[Biobase]{AnnotatedDataFrame} object. }
  \item{cdfname}{ Used to specify the name of an alternative cdf package. If set to \code{NULL}, the usual cdf package based on Affymetrix' mappings will be used. }
  \item{verbose}{ A logical value. If \code{TRUE} it writes out some messages. }
  \item{refindexname}{ The name of the array used as reference. }
  \item{refconstant}{ The constant calculated by \code{FUN} from the reference array. This value will be used for scaling. }
}

\details{
Parallelized scaling normalization of arrays. This means that all the array are scaled so that they have the same mean value.

For the serial function and more details see the function \code{normalize.constant}.

For using this function a computer cluster using the \code{snow} package has to be started.

\code{normalizeConstantParaSF1} and \code{normalizeConstantParaSF2} are internal function which will be executed at all slaves.
\item{\code{normalizeConstantParaSF1}}{Calculates and returns intensities from the refindex array. }
\item{\code{normalizeConstantParaSF2}}{Normalizes the arrays at all slaves. }
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

AffyBatch <- normalizeAffyBatchConstantPara(c1, Dilution, verbose=TRUE)

stopCluster(c1)
}
}

\keyword{programming}
\keyword{manip}
