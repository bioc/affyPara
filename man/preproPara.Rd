\name{preproPara}
\alias{preproPara}
\title{ Parallelized preprocessing }
\description{Parallelized preprocessing function, which goes from raw probe intensities to expression values in three steps: Background correction, normalization and summarization}

\usage{
preproPara(object,
    bgcorrect = TRUE, bgcorrect.method = NULL, bgcorrect.param = list(),
    normalize = TRUE, normalize.method = NULL, normalize.param = list(),
    pmcorrect.method = NULL, pmcorrect.param = list(),
    summary.method = NULL, summary.param = list(),
    ids = NULL,phenoData = new("AnnotatedDataFrame"), cdfname = NULL,
    cluster, verbose = FALSE)
}

\arguments{
  \item{object}{ An object of class \link[affy:AffyBatch-class]{AffyBatch} 
  		OR a \code{character} vector with the names of CEL files 
  				OR a (partitioned) list of \code{character} vectors with CEL file names.}
 
  \item{bgcorrect}{ A boolean to express whether background correction is wanted or not. }
  \item{bgcorrect.method}{ The name of the background adjustment method to use. }
  \item{bgcorrect.param}{ A list of parameters for \code{bgcorrect.method} (if needed/wanted) }

  \item{normalize}{ A boolean to express whether normalization is wanted or not. }
  \item{normalize.method}{ The name of the normalization method to use. }
  \item{normalize.param}{ A list of parameters to be passed to the normalization method (if wanted). }

  \item{pmcorrect.method}{ The name of the PM adjustment method. }
  \item{pmcorrect.param}{ A list of parameters for \code{pmcorrect.method} (if needed/wanted). }
  \item{summary.method}{ The method used for the computation of expression values }
  \item{summary.param}{ A list of parameters to be passed to the \code{summary.method} (if wanted). }
  \item{ids}{ List of \code{ids} for summarization }

  \item{phenoData}{ An \link[Biobase]{AnnotatedDataFrame} object. }
  \item{cdfname}{ Used to specify the name of an alternative cdf package. If set to \code{NULL}, 
  		the usual cdf package based on Affymetrix' mappings will be used. }
  \item{cluster}{ A cluster object obtained from the function \link[snow:snow-startstop]{makeCluster} in the SNOW package. 
  		For default \code{.affyParaInternalEnv$cl} will be used. }
  \item{verbose}{ A logical value. If \code{TRUE} it writes out some messages.} 
}

\details{
Parallelized preprocessing function, which goes from raw probe intensities to expression values in three steps: 
Background correction, normalization and summarization

For the serial function and more details see the function \code{expresso}.

For using this function a computer cluster using the SNOW package has to be started. 
Starting the cluster with the command \code{makeCluster} generates an cluster object in the affyPara environment (.affyParaInternalEnv) and 
no cluster object in the global environment. The cluster object in the affyPara environment will be used as default cluster object, 
therefore no more cluster object handling is required.  
The \code{makeXXXcluster} functions from the package SNOW can be used to create an cluster object in the global environment and 
to use it for the preprocessing functions.

Available methods:
\describe{
  \item{\code{bgcorrect.method}:}{see \code{bgcorrect.methods()}}
  \item{\code{normalize.method}:}{'quantil', 'constant', 'invariantset','loess'}
  \item{\code{summary.method}:}{see \code{generateExprSet.methods()}}
}
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

  esset <- preproPara(Dilution,
    bgcorrect = TRUE, bgcorrect.method = "rma2",
    normalize = TRUE, normalize.method = "quantil",
    pmcorrect.method = "pmonly",
    summary.method = "avgdiff",
    verbose = TRUE)
    
  stopCluster()
}
}
}

\keyword{programming}
\keyword{manip}
