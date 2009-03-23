\name{removeDistributedFiles}
\alias{removeDistributedFiles}
\title{Remove distributed files from slaves}

\description{
This function removes distributed files from a special path at the disk at all slaves in a computer cluster.
}

\usage{
removeDistributedFiles(path=tempdir(), cluster, master=TRUE, verbose = getOption("verbose"))
}
\arguments{
  \item{path}{ A \code{character} that defines which path (inclusive files) should be removed at every slave. Default: tempdir() }
  \item{cluster}{ A cluster object obtained from the function \link[snow:snow-startstop]{makeCluster} in the SNOW package. 
  		For default \code{.affyParaInternalEnv$cl} will be used. }
  \item{master}{ A logical value. If \code{TRUE} the files will be removed from the master. default: TRUE}
  \item{verbose}{ A logical value. If \code{TRUE} it writes out some messages. default: getOption("verbose") } 
}
\details{
This function removes distributed files from a special path at the disk at all slaves in a computer cluster.

For using this function a computer cluster using the SNOW package has to be started. 
Starting the cluster with the command \code{makeCluster} generates an cluster object in the affyPara environment (.affyParaInternalEnv) and 
no cluster object in the global environment. The cluster object in the affyPara environment will be used as default cluster object, 
therefore no more cluster object handling is required.  
The \code{makeXXXcluster} functions from the package SNOW can be used to create an cluster object in the global environment and 
to use it for the preprocessing functions.
}

\value{
If \code{verbose = TRUE}, result of removing (successfully / not successfully) will be noticed with a message.
}

\author{ Markus Schmidberger \email{schmidb@ibe.med.uni-muenchen.de}, Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de} }

\examples{
\dontrun{
library(affyPara)

makeCluster(10)

removeDistributedFiles(verbose=TRUE)

stopCluster()
}
}

\keyword{programming}