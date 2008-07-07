\name{distributeFiles}
\alias{distributeFiles}
\alias{hirarchicalDistSF}
\title{Distribute files to slaves}

\description{
This function distributes files from the master node to the disk of the slaves in the computer cluster.
}

\usage{
distributeFiles(cluster, 
	files, to = "/usr1/tmp/CELfiles",
	protocol = c("R", "RCP", "SCP"), hierarchicallyDist = FALSE, 
	master=TRUE, delExistTo=FALSE,
	full.names=FALSE, verbose = FALSE)
	
hirarchicalDistSF(to, nodes, protocol)
}

\arguments{
  \item{cluster}{ A cluster object obtained from the function \link[snow:snow-startstop]{makeCluster} in the \code{SNOW} package. }
  \item{files}{ A \code{character} vector containing the names of the files. }
  \item{to}{ A \code{character} that defines the path where the files should be stored at the slaves. Default: /usr1/tmp/CELfiles }
  \item{protocol}{ A \code{character} that defines the Copy-Protocol: "R", "RCP", "SCP" }
  \item{hierarchicallyDist}{ A logical value. If \code{TRUE} data will be hierarchically distributed to all slaves. If \code{FALSE} at every slave only a part of data is available. }
  \item{master}{ A logical value. If \code{TRUE} all data will be copied to the 'to' directory at the master node. Default = TRUE }
  \item{delExistTo}{ A logical value. If \code{TRUE} directory 'to' will be deleted at master and all nodes first. Default = FALSE}
  \item{full.names}{ A logical value. If \code{TRUE}, the directory path is prepended to the file names. If \code{FALSE}, only the file names are returned . }
  \item{verbose}{ A logical value. If \code{TRUE} it writes out some messages.} 
  \item{nodes}{ A list of \code{character} with the node names where the files have to be send.}
}
\details{
This function distributes files from the master node to the disk of the slaves in the computer cluster. First the vector of files get partitioned by the number of slaves. Then the parts will be copied to the \code{to} directory at the slaves. If \code{hierarchicallyDist} is \code{TRUE}, all slaves change the files among each other and in the end at every slave all files are located. (But this is not necessary for distributed computing with the \code{affyPara} package.

For using this function a computer cluster using the \code{SNOW} package has to be started.

\code{hirarchicalDistSF} is a internal function which will be executed at slaves.
\item{\code{hirarchicalDistSF}}{ Slavefunction for exchanging the files between all slaves.}
}

\section{Warning}{For protocol "R" hierarchically distribution not yet available.}

\value{
A list of two objects
  \item{to}{A \code{character} that defines the path where the files are located at the slaves.}
  \item{CELfiles}{A list of \code{characters}, how the files are distributed to the slaves.}
}

\author{ Markus Schmidberger \email{schmidb@ibe.med.uni-muenchen.de}, Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de} }

\examples{
\dontrun{
library(affyPara)

c1 <- makeCluster(10)

path <- "tmp/CELfiles"
CELfiles <- list.files(path,full.names=TRUE)

distList <- distributeFiles(c1, CELfiles, protocol="RCP", verbose=TRUE)

stopCluster(c1)
}
}

\keyword{programming}