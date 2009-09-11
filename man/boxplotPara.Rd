\name{boxplotPara}
\alias{boxplotPara}
\title{Parallelized Box Plots for Microarray Data}

\description{
Parallelized functions for Box Plots for Microarray Data. And optimized plots for a large number of microarrays.
}
\usage{
boxplotPara(object,
		nSample=if(length(object)> 200) nSample<-200 else nSample <- length(object),
		iqrMethod=TRUE, percent=0.05,
		typDef="mean", plot=TRUE,	 
		plotAllBoxes=TRUE,
		cluster, verbose = getOption("verbose")) 
}

\arguments{
  \item{nSample}{ A numeric value. Indicates the number of maximal samples that should be plot at the same boxplot. 
  		default: 250 }
  \item{iqrMethod}{ A logical value, if TRUE the second method will be considered to calculate the 
  		"bad" quality Samples otherwise the first Method. See Details. default: TRUE } 
  \item{percent}{ A numeric value [0.0-1.1]. If iqrMethod=TRUE the second method will be considered to 
  		calculate the "bad" quality Samples otherwise the first Method. default: 0.05 }
  \item{typDef}{ A character value. Indicates how the default Sample should be calculated. 
  		As median or mean from all Samples value. Only three possibilities: "media", "mean", c("median", "mean"). default: "mean" }
  \item{plot}{ A logical value. Indicates if graphics should be drawn. default: TRUE}
  \item{plotAllBoxes}{ A logical value. If TRUE then all Samples will be ploted, 
  		if FALSE only bad arrays will be plotted. default: TRUE }
  \item{object}{ An object of class \link[affy:AffyBatch-class]{AffyBatch} 
  		OR a \code{character} vector with the names of CEL files 
  		OR a (partitioned) list of \code{character} vectors with CEL file names.}
  \item{cluster}{ A cluster object obtained from the function \link[snow:snow-startstop]{makeCluster} in the SNOW package. 
  		For default \code{.affyParaInternalEnv$cl}  will be used. }
  \item{verbose}{ A logical value. If \code{TRUE} it writes out some messages. default: getOption("verbose") }
}
\details{
\code{boxplotPara} is the parallelized function for box plots of probe intensities. 
It is a function to check and control the Data quality of the samples using the boxplot methode. 
For serial function an more details see \link[affy:AffyBatch-class]{boxplot}. This function is optimized for huge numbers of microarray data.

For using this function a computer cluster using the SNOW package has to be started. 
Starting the cluster with the command \code{makeCluster} generates an cluster object in the affyPara environment (.affyParaInternalEnv) and 
no cluster object in the global environment. The cluster object in the affyPara environment will be used as default cluster object, 
therefore no more cluster object handling is required.  
The \code{makeXXXcluster} functions from the package SNOW can be used to create an cluster object in the global environment and 
to use it for the preprocessing functions.

We need to calculate a default Sample as reference , which has been built from all Samples data. 
Therefore the first is the calculation of the boxplot.stats, it will be made parallel at the cluster.
The calculated values are merged at the master as well as the following calculations, plots and histograms
There are two possibility to calculate the limits between the "good-bad"  quality Samples:
 1. From the differences between defaultSample values (only media, HL nad HU will be considered) and all the samples. 
	At limit will be considered as critical and it help to calculate the "bad" quality Samples.
	(it is fixed as parameter and thus do it not so sure) 
 2. From the median and IQR obtained from a boxplot, which is calculated from all Samples values. 
	The outliers of these boxplot are the "bad" quality Samples. It should be as default parameter.
}

\value{
boxplotPara returns a list with elements from the boxplot.stats function ('stats', 'n', 'conf', 'out', 'group', 'names') and 
\code{QualityPS}, \code{values_boxP} and \code{results_boxP}.
\item{qualityPS}{is a list which contains all "bad" quality Arrays classified in levels.}
\item{values_boxP}{contains the calculated differences and limits, which are used to make the classification of the Arrays in levels.}
\item{results_boxP}{summary from \code{qualityPS} as matrix, which contains only the Arrays that are considered as "bad" quality and in which levels are they classified. Possible values are 0 if the Array is not at this levels and 1 if it is classified as "bad" sample at this level}
}

\author{Esmeralda Vicedo <e.vicedo@gmx.net>, Markus Schmidberger \email{schmidb@ibe.med.uni-muenchen.de}, Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de} }

\examples{
\dontrun{
library(affyPara)
if (require(affydata)) {
  data(Dilution)

  makeCluster(3)

  ##boxplot of Dulution data (affybatch) 
  box1 <- boxplotPara(Dilution)

  ## boxplots to a pdf file
  pdf(file="boxplot.pdf", title="AffyBatch Boxplot")
  box2 <- boxplotPara(Dilution)
  dev.off()

  stopCluster()
}
}
}

\keyword{programming}
\keyword{manip}