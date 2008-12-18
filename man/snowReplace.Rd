\name{snow-startstop}
\title{Starting and Stopping SNOW Clusters}
\alias{makeCluster}
\alias{stopCluster}
\usage{
makeCluster( ...)
stopCluster(cl)
}
\arguments{
  \item{...}{cluster option specifications, for details see SNOW package.}
  \item{cl}{cluster object.}
}
\description{
  Functions to start and stop a SNOW cluster and to set default cluster
  options. Wrapper around original start-stop commands from SNOW to hide 
  the cluster object.
}
\details{
  \code{makeCluster} starts a cluster of the specified or default type
  and returns NO reference to the cluster. The reference is stored in an internal 
  environment (.affyParaInternalEnv) and will be automatically used from the functions 
  in affyPara.
  For further parameters and documentation see the SNOW package. 
  
  \code{stopCluster} should be called to properly shut down the cluster
  before exiting R.  If it is not called it may be necessary to use
  external means to ensure that all slave processes are shut down.
}

\examples{
  \dontrun{
	makeCluster(2)
	stopCluster()
  }
}
\keyword{programming}
