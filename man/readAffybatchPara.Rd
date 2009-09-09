\name{readAffybatchPara}
\alias{read.affybatchPara}
\title{Parallelized Read-AffyBatch function}
\description{
Parallelization of the read.affybatch function. This parallel implementation is especially 
useful for multicore machines.
}
\usage{
read.affybatchPara(object,
		phenoData = new("AnnotatedDataFrame"),
		description = NULL, notes = "",	
		cluster, verbose=getOption("verbose")) 
}

\arguments{
   \item{object}{ An object of a \code{character} vector with the names of CEL files 
   		OR a (partitioned) list of \code{character} vectors with CEL file names.}
   \item{phenoData}{ An \link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame} object. }
   \item{description}{a 'MIAME' object.}
   \item{notes}{notes.}

   \item{cluster}{ A cluster object obtained from the function \link[snow:snow-startstop]{makeCluster} in the SNOW package. 
  		For default \code{.affyParaInternalEnv$cl}  will be used. }
  \item{verbose}{ A logical value. If \code{TRUE} it writes out some messages. default: getOption("verbose") } 
 }
\details{
Parallelized creation of an AffyBatch object. Especially useful on multi-core machines 
to accelerate the creation of the AffyBatch object.

For the serial function and more details see the function \code{read.affybatch}.

For using this function a computer cluster using the SNOW package has to be started. 
Starting the cluster with the command \code{makeCluster} generates an cluster object in the affyPara environment (.affyParaInternalEnv) and 
no cluster object in the global environment. The cluster object in the affyPara environment will be used as default cluster object, 
therefore no more cluster object handling is required.  
The \code{makeXXXcluster} functions from the package SNOW can be used to create an cluster object in the global environment and 
to use it for the preprocessing functions.

}

\value{
An \link[affy:AffyBatch-class]{AffyBatch} object.
}

\author{ Markus Schmidberger \email{schmidb@ibe.med.uni-muenchen.de}, Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de} }

\examples{
\dontrun{
library(affyPara)
if (require(affydata)) {
	celpath <- system.file("celfiles", package="affydata")
    fns <- list.celfiles(path=celpath,full.names=TRUE)

	makeCluster(3)
	
	##read a text celfile
   	abatch <- read.affybatchPara(fns[2], verbose=TRUE)
	
	stopCluster()
}
}
}

\keyword{programming}
\keyword{manip}
