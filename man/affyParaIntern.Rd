\name{affyParaIntern}
\alias{initAffyBatchSF}
\alias{getAffyBatchSF}
\alias{getIntensitySF}
\alias{getCompIntensitySF}
\alias{getCompIntensityMatrixSF}
\alias{getFUNAffyBatchSF}
\alias{writeLinesSF}
\alias{ReadHeaderSF}
\alias{getObjectType}
\alias{checkPartSize}

\title{ Internal affyPara objects / functions }

\description{
Internal functions for the \code{affyPara} package.
}

\usage{
initAffyBatchSF(object, object.type, rm.list)

getAffyBatchSF()
getIntensitySF(rows, refindexname)
getCompIntensitySF(rows)
getCompIntensityMatrixSF(rows, drop=FALSE)
getFUNAffyBatchSF(FUN)

writeLinesSF(data, fileName)
ReadHeaderSF(object)

getObjectType(object)
checkPartSize(object, number.parts)
}

\arguments{
  \item{object}{ An object of class \link[affy:AffyBatch-class]{AffyBatch} OR a \code{character} vector with the names of CEL files OR a (partitioned) list of \code{character} vectors with CEL file names.}
  \item{object.type}{ Declaration for the type of the argument \code{object}: "AffyBatch", "CELfileVec", "partCELfileList" }
  \item{rm.list}{ A logical value. If \code{TRUE} all variables in the global environment will be removed before initialization of the AffyBatch}
  \item{rows}{ Number of rows which have to be changed. }
  \item{refindexname}{ The name of the array used as a reference. }
  \item{drop}{ A logical value. If \code{TRUE} the dimensions of an array which have only one level will be deleted. } 
  \item{FUN}{ A function generating a value from a an \link[affy:AffyBatch-class]{AffyBatch}. e.g. \code{dim} }
  \item{data}{ A \code{character} vector, containing the data for the file.}
  \item{fileName}{ A \code{character} string with the file name.}
  \item{number.parts}{ Number of nodes in the computer cluster.}
}
\details{
Internal functions for the \code{affyPara} package. The functions have to be used in a cluster function (e.g. \link[snow:snow-cluster]{clusterApply}) from the SNOW package.
\item{\code{initAffyBatchSF}}{ Slavefunction for initializing an \link[affy:AffyBatch-class]{AffyBatch} at slaves. \link[affy:AffyBatch-class]{AffyBatch} will be stored in the \link[base:environment]{.GlobalEnv} with the name 'AffyBatch'.}

\item{\code{getAffyBatchSF}}{ Slavefunction to get an \link[affy:AffyBatch-class]{AffyBatch} from slaves. Gets the object \link[affy:AffyBatch-class]{AffyBatch} from the \link[base:environment]{.GlobalEnv}. }
\item{\code{getIntensitySF}}{ Slavefunction to get special values from the intensity matrix from slaves.}
\item{\code{getCompIntensitySF}}{ Slavefunction to get special rows from the intensity matrix from slaves.}
\item{\code{getCompIntensityMatrixSF}}{ Slavefunction to get complete intensity matrix from slaves.}
\item{\code{getFUNAffyBatchSF}}{ Slavefunction to get a value from an \link[affy:AffyBatch-class]{AffyBatch} at slaves.}

\item{\code{writeLinesSF}}{ Slavefunction to write data into a file at slaves.}
\item{\code{ReadHeaderSF}}{ Slavefunction to return Header-Informations from CEL Files at slaves. }

\item{\code{getObjectType}}{ Function to get type from object.}
\item{\code{checkPartSize}}{ Function to check object for length.}
}

\author{ Markus Schmidberger \email{schmidb@ibe.med.uni-muenchen.de}, Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de} }

\keyword{programming}