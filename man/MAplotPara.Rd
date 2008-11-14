\name{MAplotPara}
\alias{MAplotPara}

\title{Parallelized relative M vs. A plots for Microarray Data}
\description{
  Parallelized creation of M vs A plots. Where M is determined relative to a specified chip or to a pseudo-median reference chip.
  And optimized plots for a large number of microarrays.
}
\usage{
MAplotPara(cluster,
 			object, 
            log=TRUE, type=c("both","pm","mm"),
            ref=NULL, which=NULL, 
            ref.title="vs mean ref.Array", subset=NULL,
            span=1/4, family.loess="gaussian",
            plot.method=c("normal","smoothScatter","add"),
            pch=".", plot =TRUE,
            cutoff =0.5,# add parameter to generic function ma.plot
            level=1,
            verbose=FALSE,
            ...)
}

\arguments{
  \item{cluster}{ A cluster object obtained from the function \link[snow:snow-startstop]{makeCluster} in the \code{SNOW} package. }
  \item{object}{ An object of class \link[affy:AffyBatch-class]{AffyBatch} OR a \code{character} vector with the names of CEL files OR a (partitioned) list of \code{character} vectors with CEL file names.}
  \item{...}{Additional parameters for the routine}
  \item{log}{ ??}
  \item{type}{ ??}
  \item{ref}{ ??}
  \item{which}{ ??}
  \item{ref.title}{ ??}
  \item{subset}{A set of indices to use when drawing the loess curve}
  \item{span}{span to be used for loess fit.}
  \item{family.loess}{\code{"guassian"} or \code{"symmetric"} as in \code{\link[modreg]{loess}}.}
  \item{plot.method}{a string specifying how the plot is to be
    drawn. \code{"normal"} plots points, \code{"smoothScatter"} uses the
  \code{\link[geneplotter]{smoothScatter}} function. Specifying
  \code{"add"} means that the MAplot should be added to the current
  plot}
  \item{pch}{ ??}
  \item{plot}{ A logical value. TRUE as default, MAplots will be drawn otherwise only the "bad" quality arrays will be calculated.}
  \item{cutoff}{??}
  \item{level}{level- numerical - indicates which level of "bad" quality arrays should be ploted if plotDraw =TRUE:
  	1 - only first level "bad" quality  will be considered. First level "bad" array quality are the arrays considered as "bad" after the three possible parameter: S, loess, and sigma 
    2 -  first level  "bad" quality and second level  will be considered. Second level "bad" quality Arrays are the arrays which has been classified as bad after two of the three possible parameter 
    3 -  all levels will be ploted : first, second and third. Third level "bad" quality Arrays are the arrays which are considered as "bad" after one of the three parameter.
  }
  \item{verbose}{ A logical value. If \code{TRUE} it writes out some messages. default: TRUE }
}
\details{
MAplotPara is a function based on the generic function ma.plot from affy package.
Only the following parameters are original for MAplotPAra: cluster, object, cuttof, plot, level
The parameter ref.fn=c("median","mean") is not allowed because it is not possible to calculate the reference median array as parallelized.

\code{MAplotPara} is the parallelized function for MA plots of probe intensities. It is a function to check and control the Data quality of the samples using the MA plot method. For serial function an more details see \link[affy]{boxplot}. This function is optimized for huge numbers of microarray data.
For using this function a computer cluster using the \code{snow} package has to be started.
}

\author{Esmeralda Vicedo <e.vicedo@gmx.net>, Markus Schmidberger \email{schmidb@ibe.med.uni-muenchen.de}, Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de} }

\examples{
\dontrun{
library(affyPara)

data(Dilution)

c1 <- makeCluster(3)

##MA of Dilution data (affybatch)
ma1 <- MAplotPara(c1, Dilution)

## MAplot to a pdf file
pdf(file="maplot.pdf", title="AffyBatch MAplot")
ma2 <- MAplotPara(c1, Dilution)
dev.off()

stopCluster(c1)
}
}

\keyword{programming}
\keyword{manip}