\name{vsnInputPara}
\docType{class}
\alias{class:vsnInputPara}
\alias{vsnInputPara}
\alias{vsnInputPara-class}

\alias{dim,vsnInputPara-method} 
\alias{nrow,vsnInputPara-method} 
\alias{ncol,vsnInputPara-method} 

\title{
  Class to contain input data and parameters for parallel vsn functions
}

\description{
  Class extends the class \code{vsnInput}. The class contains input data and parameters for parallel vsn functions.
}

\section{Creating Objects}{
  \code{new("vsnInputPara")}
}

\section{Slots}{
 \describe{
  \item{\code{...}:}{ as class \code{vsnInput}.}
  \item{\code{dimAB}:}{ The dimension of the complete AffyBatch.}
  }
}


\section{Methods}{
  \describe{
    \item{\code{dim}}{Get dimensions of data matrix.}
    \item{\code{nrow}}{Get number of rows of data matrix.}
    \item{\code{ncol}}{Get number of columns of data matrix.}
}}


\author{ Markus Schmidberger \email{schmidb@ibe.med.uni-muenchen.de}, Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de} }

\seealso{
  \code{vsn2}
}

\keyword{classes}
