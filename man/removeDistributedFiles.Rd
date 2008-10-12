\name{removeDistributedFiles}
\alias{removeDistributedFiles}
\alias{removeFilesSF}
\alias{accessFilesSF}
\title{Remove distributed files from slaves}

\description{
This function removes distributed files from a special path at the disk at all slaves in a computer cluster.
}

\usage{
removeDistributedFiles(cluster, path=tempdir(), verbose = FALSE)
}
\arguments{
  \item{cluster}{ A cluster object obtained from the function \link[snow:snow-startstop]{makeCluster} in the \code{SNOW} package. }
  \item{path}{ A \code{character} that defines which path (inclusive files) should be removed at every slave. Default: tempdir() }
  \item{verbose}{ A logical value. If \code{TRUE} it writes out some messages.}
}
\details{
This function removes distributed files from a special path at the disk at all slaves in a computer cluster.

For using this function a computer cluster using the \code{SNOW} package has to be started.
}

\value{
If \code{verbose = TRUE}, result of removing (successfully / not successfully) will be noticed with a message.
}

\author{ Markus Schmidberger \email{schmidb@ibe.med.uni-muenchen.de}, Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de} }

\examples{
\dontrun{
library(affyPara)

c1 <- makeCluster(10)

removeDistributedFiles(c1, verbose=TRUE)

stopCluster(c1)
}
}

\keyword{programming}