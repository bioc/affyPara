\name{affyParaIntern}
\alias{initAffyBatchSF}
\alias{getAffyBatchSF}
\alias{getIntensitySF}
\alias{getCompIntensitySF}
\alias{getCompIntensityMatrixSF}
\alias{getFUNAffyBatchSF}
\alias{createDirSF}
\alias{writeLinesSF}
\alias{listFilesSF}
\alias{file.name}
\alias{ReadHeader}
\alias{getObjectType}
\alias{checkPartSize}

\title{ Internal affyPara objects / functions }

\description{
Internal functions for the \code{affyPara} package.
}

\usage{
initAffyBatchSF(object, object.type)

getAffyBatchSF()
getIntensitySF(rows, refindexname)
getCompIntensitySF(rows)
getCompIntensityMatrixSF(rows, drop=FALSE)
getFUNAffyBatchSF(FUN)

createDirSF(to, showWarnings = TRUE, recursive = TRUE)
writeLinesSF(data, fileName)
listFilesSF(path, full.names=TRUE)
ReadHeader(object)

file.name(file, fsep = .Platform$file.sep)

getObjectType(object)
checkPartSize(object, number.parts)
}

\arguments{
  \item{object}{ An object of class \link[affy:AffyBatch-class]{AffyBatch} OR a \code{character} vector with the names of CEL files OR a (partitioned) list of \code{character} vectors with CEL file names.}
  \item{object.type}{ Declaration for the type of the argument \code{object}: "AffyBatch", "CELfileVec", "partCELfileList" }
  \item{rows}{ Number of rows which have to be changed. }
  \item{refindexname}{ The name of the array used as a reference. }
  \item{drop}{ A logical value. If \code{TRUE} the dimensions of an array which have only one level will be deleted. } 
  \item{FUN}{ A function generating a value from a an \link[affy:AffyBatch-class]{AffyBatch}. e.g. \code{dim} }
  \item{to}{ A \code{character} vector, that defines the path and name of the new directory.}
  \item{showWarnings}{ A logical value; should the warnings on failure be shown?}
  \item{recursive}{ A logical values; should elements of the path other than the last be created? If \code{TRUE}, like Unix's mkdir -p.}
  \item{data}{ A \code{character} vector, containing the data for the file.}
  \item{fileName}{ A \code{character} string with the file name.}
  \item{path}{ A \code{character} vector with the full directory.}
  \item{full.names}{ A logical value. If \code{TRUE}, the directory path is prepended to the file names. If \code{FALSE}, only the file names are returned.}
  \item{file}{ A \code{character} vector of file names with prepended directory path.}
  \item{fsep}{ A \code{character} string, giving the file separator used on your platform. }
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

\item{\code{createDirSF}}{ Slavefunction to create a directory at slaves.}
\item{\code{writeLinesSF}}{ Slavefunction to write data into a file at slaves.}
\item{\code{listFilesSF}}{ Slavefunction to list files from directory at slaves.}
\item{\code{ReadHeader}}{ Slavefunction to return Header-Informations from CEL Files at slaves. }

\item{\code{file.name}}{ Function to get file name from a \code{character} with complete file path (directory path and filename).}

\item{\code{getObjectType}}{ Function to get type from object.}
\item{\code{checkPartSize}}{ Function to check object for length.}
}

\author{ Markus Schmidberger \email{schmidb@ibe.med.uni-muenchen.de}, Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de} }

\keyword{programming}