\name{computeExprSetPara}
\alias{computeExprSetPara}

\title{ Parallel generate a set of expression values }
\description{
Parallel generation of a set of expression values from the probe pair information. The set of expression is returned as an ExpressionSet object. 
}
\usage{
computeExprSetPara(cluster,
    object,
    ids = NULL,
    pmcorrect.method, summary.method,
    summary.param = list(), pmcorrect.param = list(),
    phenoData = new("AnnotatedDataFrame"), cdfname = NULL,
    verbose = TRUE)

}
\arguments{
  \item{cluster}{ A cluster object obtained from the function \link[snow:snow-startstop]{makeCluster} in the \code{SNOW} package. }
  \item{object}{ An object of class \link[affy:AffyBatch-class]{AffyBatch} OR a \code{character} vector with the names of CEL files OR a (partitioned) list of \code{character} vectors with CEL file names.}
 
  \item{pmcorrect.method}{ The name of the PM adjustement method. }
  \item{pmcorrect.param}{ A list of parameters for \code{pmcorrect.method} (if needed/wanted). }
  \item{summary.method}{ The method used for the computation of expression values }
  \item{summary.param}{ A list of parameters to be passed to the \code{summary.method} (if wanted). }
  \item{ids}{ List of \code{ids} for summarization }

  \item{phenoData}{ An \link[Biobase]{AnnotatedDataFrame} object. }
  \item{cdfname}{ Used to specify the name of an alternative cdf package. If set to \code{NULL}, the usual cdf package based on Affymetrix' mappings will be used. }
  \item{verbose}{ A logical value. If \code{TRUE} it writes out some messages. }
}
\details{
Parallelized preprocessing function, which goes from raw probe intensities to expression values in one steps: summarization

For the serial function and more details see the function \code{computeExprSet}.

For using this function a computer cluster using the \code{snow} package has to be started.
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

esset <- computeExprSetPara(cluster,
    affybatch.example,
    pmcorrect.method = "pmonly",
    summary.method = "avgdiff",
    verbose = TRUE)

stopCluster(c1)
}
}

\keyword{programming}
\keyword{manip}
