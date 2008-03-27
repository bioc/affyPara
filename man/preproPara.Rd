\name{preproPara}
\alias{preproPara}
\title{ Parallelized preprocessing }
\description{Parallelized preprocessing function, which goes from raw probe intensities to expression values in three steps: Background correction, normalization and summarization}

\usage{
preproPara(cluster,
    object,
    bgcorrect = TRUE, bgcorrect.method = NULL, bgcorrect.param = list(),
    normalize = TRUE, normalize.method = NULL, normalize.param = list(),
    pmcorrect.method = NULL, pmcorrect.param = list(),
    summary.method = NULL, summary.param = list(),
    ids = NULL,phenoData = new("AnnotatedDataFrame"), cdfname = NULL,
    verbose = FALSE)
}

\arguments{
  \item{cluster}{ A cluster object obtained from the function \link[snow:snow-startstop]{makeCluster} in the \code{SNOW} package. }
  \item{object}{ An object of class \link[affy:AffyBatch-class]{AffyBatch} OR a \code{character} vector with the names of CEL files OR a (partitioned) list of \code{character} vectors with CEL file names.}
 
  \item{bgcorrect}{ A boolean to express whether background correction is wanted or not. }
  \item{bgcorrect.method}{ The name of the background adjustment method to use. }
  \item{bgcorrect.param}{ A list of parameters for \code{bgcorrect.method} (if needed/wanted) }

  \item{normalize}{ A boolean to express whether normalization is wanted or not. }
  \item{normalize.method}{ The name of the normalization method to use. }
  \item{normalize.param}{ A list of parameters to be passed to the normalization method (if wanted). }

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
Parallelized preprocessing function, which goes from raw probe intensities to expression values in three steps: Background correction, normalization and summarization

For the serial function and more details see the function \code{expresso}.

For using this function a computer cluster using the \code{snow} package has to be started.

Available methods:
\item{\code{bgcorrect.method}:}{see \code{bgcorrect.methods}}
\item{\code{normalize.method}:}{'quantil', 'constant', 'invariantset'}
\item{\code{summary.method}:}{see \code{generateExprSet.methods}}
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

esset <- preproPara(cluster,
    affybatch.example,
    bgcorrect = TRUE, bgcorrect.method = "rma2",
    normalize = TRUE, normalize.method = "quantil",
    pmcorrect.method = "pmonly",
    summary.method = "avgdiff",
    verbose = TRUE)
    
stopCluster(c1)
}
}

\keyword{programming}
\keyword{manip}
