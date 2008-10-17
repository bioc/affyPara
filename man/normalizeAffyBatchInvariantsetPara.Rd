\name{normalizeAffyBatchInvariantsetPara}
\alias{normalizeAffyBatchInvariantsetPara}
\alias{normalizeInvariantsetParaSF1}
\alias{normalizeInvariantsetParaSF2}
\title{Parallelized Invariante Set normalization}
\description{
Parallelized normalization of arrays using an invariant set.
}
\usage{
normalizeAffyBatchInvariantsetPara(cluster,
    object,
    prd.td = c(0.003, 0.007), baseline.type = c("mean", "median", "pseudo-mean", "pseudo-median"),
    type = c("separate", "pmonly", "mmonly", "together"),
    phenoData = new("AnnotatedDataFrame"), cdfname = NULL,
    verbose = FALSE)
    
    normalizeInvariantsetParaSF1(FUN, type)
    normalizeInvariantsetParaSF2(refindexname, rows, prd.td, baseline.chip)
}

\arguments{
  \item{cluster}{ A cluster object obtained from the function \link[snow:snow-startstop]{makeCluster} in the \code{SNOW} package. }
  \item{object}{ An object of class \link[affy:AffyBatch-class]{AffyBatch} OR a \code{character} vector with the names of CEL files OR a (partitioned) list of \code{character} vectors with CEL file names.}
 
  \item{prd.td}{ A cutoff parameter for normalization. }
  \item{baseline.type}{ Specify how to determine the baseline array (mean, median). }
  \item{type}{ A string specifying how the normalization should be applied. }
  \item{phenoData}{ A \link[Biobase]{AnnotatedDataFrame} object. }
  \item{cdfname}{ Used to specify the name of an alternative cdf package. If set to \code{NULL}, the usual cdf package based on Affymetrix' mappings will be used. }
  \item{verbose}{ A logical value. If \code{TRUE} it writes out some messages. }
  
  \item{FUN}{A function generating a value from the intensities on an array. Typically mean or median.}
  \item{refindexname}{ The name of the array used as reference.}
  \item{rows}{ Number of rows which have to be changed.}
  \item{baseline.chip}{ The values of the array used as reference.}
}

\details{
Parallelized normalization of arrays using an invariant set. The set of invariant intensities between data and ref is found through an iterative process (based on the respective ranks the intensities). This set of intensities is used to generate a normalization curve by smoothing.

For the serial function and more details see the function \code{normalize.invariantset}.

For using this function a computer cluster using the \code{snow} package has to be started.

\code{normalizeInvariantsetParaSF1} and \code{normalizeInvariantsetParaSF2} are internal function which will be executed at all slaves.
\item{\code{normalizeInvariantsetParaSF1}}{Calculates and returns a value from the intensities of an array.}
\item{\code{normalizeInvariantsetParaSF2}}{Normalizes the arrays at slaves.}
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

AffyBatch <- normalizeAffyBatchInvariantsetPara(c1, Dilution, verbose=TRUE)

stopCluster(c1)
}
}

\keyword{programming}
\keyword{manip}
