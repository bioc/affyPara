\name{boxplotPara}
\alias{boxplotPara}
\title{Parallelized Box Plots for Microarray Data}

\description{
Parallelized functions for Box Plots for Microarray Data.
}
\usage{
boxplotPara(cluster,
		object,
		nSample=if(length(object)> 200) nSample<-200 else nSample <- length(object),
		iqrMethod=TRUE, percent=0.05,
		typDef="mean", plot=TRUE,	 
		plotAllBoxes=TRUE, verbose=FALSE) 
}

\arguments{
  \item{cluster}{ A cluster object obtained from the function \link[snow:snow-startstop]{makeCluster} in the \code{SNOW} package. }
  \item{object}{ An object of class \link[affy:AffyBatch-class]{AffyBatch} OR a \code{character} vector with the names of CEL files OR a (partitioned) list of \code{character} vectors with CEL file names.}
  \item{nSample}{ A numeric value. Indicates the number of maximal samples that should be ploted at the same boxplot. default: 250 }
  \item{iqrMethod}{ A logical value, if TRUE the second method will be considered to calculate the "bad" quality Samples otherwise the first Method. See Details. default: TRUE } 
  \item{percent}{ A numeric value [0.0-1.1]. If iqrMethod=TRUE the second method will be considered to calculate the "bad" quality Samples otherwise the first Method. default: 0.05 }
  \item{typDef}{ A character value. Indicates how the default Sample should be calculated. As median or mean from all Samples value. Only three possibilities: "media", "mean", c("median", "mean"). default: "mean" }
  \item{plot}{ A logical value. Indicates if graphics should be drawn. default: TRUE}
  \item{plotAllBoxes}{ A logical value. If TRUE then all Samples will be ploted, if FALSE only bad arrays will be plotted. default: TRUE }
  \item{verbose}{ A logical value. If \code{TRUE} it writes out some messages. default: TRUE }
}
\details{
\code{boxplotPara} is the parallelized function for box plots of probe intensities. It is a function to check and control the Data quality of the samples using the boxplot methode. For serial function an more details see \link[affy]{boxplot}. This function is optimized for huge numbers of microarray data.

For using this function a computer cluster using the \code{snow} package has to be started.

We need to calculate a default Sample as reference , which has been built from all Samples data. 
Therefore the first is the calculation of the boxplot.stats, it will be made parallel at the cluster.
The calculated values are merged at the master as well as the following calculations, plots and histograms
There are two possibility to calculate the limits between the "good-bad"  quality Samples:
1. From the differences between defaultSample values(only media, HL nad HU will be considered) and all the samples. At limit will be considered as critical and it help to calculate the "bad" quality Samples.(it is fixed as parameter and thus do it not so sure. 
2. From the median and IQR obtained from a boxplot, which is calculated from all Samples values. The outliers of these boxplot are the "bad" quality Samples. It should be as default parameter.

}

\author{Esmeralda Vicedo <e.vicedo@gmx.net>, Markus Schmidberger \email{schmidb@ibe.med.uni-muenchen.de}, Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de} }

\examples{
\dontrun{
library(affyPara)

data(Dilution)

c1 <- makeCluster(3)

##boxplot of affybatch.example 
box1 <- boxplotPara(c1, Dilution)

## boxplots to a pdf file
pdf(file="boxplot.pdf", title="AffyBatch Boxplot")
box3 <- boxplotPara(c1, affybatch.example)
dev.off()

stopCluster(c1)
}
}

\keyword{programming}
\keyword{manip}