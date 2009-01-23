\name{bgCorrectPara}
\alias{bgCorrectPara}
\title{Parallelized Background Correction}

\description{
Parallelized functions for background correction of probe intensities.
}
\usage{
bgCorrectPara(object,
	phenoData = new("AnnotatedDataFrame"), method,
	cluster, verbose = FALSE)
}

\arguments{
  \item{object}{ An object of class \link[affy:AffyBatch-class]{AffyBatch} OR 
  		a \code{character} vector with the names of CEL files 
  		OR a (partitioned) list of \code{character} vectors with CEL file names.}
  \item{phenoData}{ An \link[Biobase]{AnnotatedDataFrame} object }
  \item{method}{ A \code{character} that defines what background correction method will be used. 
  		Available methods are given by \code{bg.correct.methods}. 
  		The name of the method to apply must be double-quoted. }
  \item{cluster}{ A cluster object obtained from the function \link[snow:snow-startstop]{makeCluster} in the SNOW package. 
  		For default \code{.affyParaInternalEnv$cl}  will be used. }
  \item{verbose}{ A logical value. If \code{TRUE} it writes out some messages. }
}
\details{
\code{bgCorrectPara} is the parallelized function for background correction of probe intensities. 
For serial function an more details see \link[affy]{bg.correct}.

For using this function a computer cluster using the SNOW package has to be started. 
Starting the cluster with the command \code{makeCluster} generates an cluster object in the affyPara environment (.affyParaInternalEnv) and 
no cluster object in the global environment. The cluster object in the affyPara environment will be used as default cluster object, 
therefore no more cluster object handling is required.  
The \code{makeXXXcluster} functions from the package SNOW can be used to create an cluster object in the global environment and 
to use it for the preprocessing functions.
}

\value{
  An \code{\link[affy:AffyBatch-class]{AffyBatch}} for which the intensities have been background adjusted. For some methods (RMA), only PMs are corrected and the MMs remain the same. 
}

\author{ Markus Schmidberger \email{schmidb@ibe.med.uni-muenchen.de}, Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de} }

\examples{
\dontrun{
library(affyPara)
if (require(affydata)) {
  data(Dilution)

  makeCluster(3)

  ##bgc will be the bg corrected version of Dilution 
  bgc <- bgCorrectPara(Dilution, method="rma", verbose=TRUE)

  stopCluster()
}
}
}

\keyword{programming}
\keyword{manip}
