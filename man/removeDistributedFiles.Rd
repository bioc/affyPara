\name{removeDistributedFiles}
\alias{removeDistributedFiles}
\alias{removeFilesSF}
\alias{accessFilesSF}
\title{Remove distributed files from slaves}

\description{
This function removes distributed files from a special path at the disk at all slaves in a computer cluster.
}

\usage{
removeDistributedFiles(cluster, path="/usr1/tmp/CELfiles", verbose = FALSE)

removeFilesSF(path, recursive = TRUE)
accessFilesSF(path, mode = 0)

}
\arguments{
  \item{cluster}{ A cluster object obtained from the function \link[snow:snow-startstop]{makeCluster} in the \code{SNOW} package. }
  \item{path}{ A \code{character} that defines which path (inclusive files) should be removed at every slave. Default: /usr1/tmp/CELfiles }
  \item{verbose}{ A logical value. If \code{TRUE} it writes out some messages.}
  \item{recursive}{ A logical value. Should directories be deleted recursively? }
  \item{mode}{ A integer specifying access mode required. See \link[base]{file.access}}
}
\details{
This function removes distributed files from a special path at the disk at all slaves in a computer cluster.

For using this function a computer cluster using the \code{SNOW} package has to be started.

\code{removeFilesSF} and \code{accessFilesSF} are internal functions which will be executed at slaves.
\item{\code{removeFilesSF}}{ Slavefunction for deleting the files(s) or directories specified by argument \code{path}. }
\item{\code{accessFilesSF}}{ Slavefunction to access information about files on slaves. Function will be used for checks. }
}

\value{
If \code{verbose = TRUE}, result of removing (successfully / not successfully) will be noticed with a message.
\item{removeFilesSF}{ Returns a logical value. \code{TRUE} for success. }
\item{accessFilesSF}{ Returns \code{0} for success and \code{1} for failure. }
}

\author{ Markus Schmidberger \email{schmidb@ibe.med.uni-muenchen.de}, Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de} }

\examples{
\dontrun{
library(affyPara)

c1 <- makeCluster(10)

removeDistributedFiles(c1, "/usr1/tmp/CELfiles", verbose=TRUE)

stopCluster(c1)
}
}

\keyword{programming}
