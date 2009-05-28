\name{qa}
\alias{summaryM1M2Para}
\title{Parallel Quality Assessment Summary}

\description{
Creates a Summary Matrix from parallel quality assessment results.
}
\usage{
summaryM1M2Para(method1, method2, 
		level, verbose=FALSE)
}

\arguments{
	\item{method1}{Result object form boxplotPara.}
	\item{method2}{Result object from MAplotPara.}
	\item{level}{level- numerical - indicates which level of "bad" quality arrays should be plot if plotDraw =TRUE:
  	1 - only first level "bad" quality  will be considered. First level "bad" array quality are the arrays considered as "bad" after the three possible parameter: S, loess, and sigma 
    2 - first level  "bad" quality and second level  will be considered. Second level "bad" quality Arrays are the arrays which has been classified as bad after two of the three possible parameter 
    3 - all levels will be plot : first, second and third. Third level "bad" quality Arrays are the arrays which are considered as "bad" after one of the three parameter.}
  \item{verbose}{ A logical value. If \code{TRUE} it writes out some messages. default: getOption("verbose") }
}
\details{
\code{summaryM1M2Para} creates a Summary Matrix from parallel quality assessment results. In the 
rows there are the arrays and in the colums the qa-methods: 0 = good quality, 1 = bad quality.

If the rowSum is bigger than 2, than the arrays should be considered as bad quality.
}

\value{
A matrix of all arrays (rows) and qa-methods (colums): 0 = good quality, 1 = bad quality
}

\author{Esmeralda Vicedo <e.vicedo@gmx.net>, Markus Schmidberger \email{schmidb@ibe.med.uni-muenchen.de}, Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de} }

\examples{
\dontrun{
library(affyPara)
if (require(affydata)) {
  data(Dilution)

  makeCluster(3, type='MPI')

  box1 <- boxplotPara(Dilution)
  ma1 <- MAplotPara(Dilution)
  
  summaryM1M2Para(box1, ma1, level=3)

  stopCluster()
}
}
}

\keyword{programming}
\keyword{manip}