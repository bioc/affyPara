\name{MAplotPara}
\alias{MAplotPara}
\alias{MAplotSer}
\title{Parallelized relative M vs. A plots for Microarray Data}
\description{
  Parallelized creation of M vs A plots. Where M is determined relative to a specified chip or to a pseudo-median reference chip.
  And optimized plots for a large number of microarrays.
}
\usage{
MAplotPara(object, 
                       log=TRUE, 
                       type=c("both","pm","mm"),
                       ref=NULL,
                       which=NULL, 
                       ref.title="vs mean ref.Array",
                       subset=NULL,
                       span=1/4,
                       show.statistics=TRUE,
                       family.loess ="gaussian", 
                       pchs=".",
                       plot =TRUE,
                       cutoff =0.5,# add parameter to generic function ma.plot
                       level=1,
                       cluster, verbose = getOption("verbose"),
                       ... )
                       
}

\arguments{
	\item{object}{An object of class \link[affy:AffyBatch-class]{AffyBatch} 
			OR a \code{character} vector with the names of CEL files 
			OR a (partitioned) list of \code{character} vectors with CEL file names.}
  	\item{...}{Additional parameters for the routine.}
  	\item{log}{A logical value. TRUE as default. Computes logarithms of Array's intensities.}
  	\item{type}{Defines how the Affimetrix Array's intensities should be considered 
  			(only for Affymetrix: perfect match(pm), mismatch or both).}
  	\item{ref}{gives the possibility to define some array of the Affybatch object as 'medianchip' 
  			otherwise it will be calculate as reference median array.}
  	\item{which}{A list of index samples from the object class, which indicates the samples to be ploted.}
  	\item{ref.title}{character- gives a 'title' to label the plots by ma.plot function.}
  	\item{subset}{A set of indices to use when drawing the loess curve.}
  	\item{span}{span to be used for loess fit.}
  	\item{family.loess}{\code{"guassian"} or \code{"symmetric"} as in \link[stats]{loess}.}
  	\item{show.statistics}{A logical value. TRUE as default. 
  			If true some summary statistics of the M values are drawn.}
  	\item{pchs}{Graphical plotting 'character'. Default Value \code{"."} equivalently to \code{pchs = 46} 
  			from the function \link[graphics]{points} }
  	\item{plot}{ A logical value. TRUE as default, MAplots will be drawn otherwise only the "bad" quality arrays will be calculated.}
  	\item{cutoff}{numerical between [0.0-1.0]. As default 0.5. It is considered as limit for the sigma parameter, 
  			which is one of the three classificators for the 'bad' quality arrays by MAplotsPara.}
  	\item{level}{level- numerical - indicates which level of "bad" quality arrays should be plot if plotDraw =TRUE:
  	1 - only first level "bad" quality  will be considered. First level "bad" array quality are the arrays considered as "bad" after the three possible parameter: S, loess, and sigma 
    2 - first level  "bad" quality and second level  will be considered. Second level "bad" quality Arrays are the arrays which has been classified as bad after two of the three possible parameter 
    3 - all levels will be plot : first, second and third. Third level "bad" quality Arrays are the arrays which are considered as "bad" after one of the three parameter.
  }
    \item{cluster}{ A cluster object obtained from the function \link[snow:snow-startstop]{makeCluster} in the SNOW package. 
  		For default \code{.affyParaInternalEnv$cl}  will be used. }
  	\item{verbose}{ A logical value. If \code{TRUE} it writes out some messages. default: getOption("verbose") }
}

\details{
MAplotPara is a function based on the generic function ma.plot from affy package.
Only the following parameters are original for MAplotPAra: cluster, object, cuttof, plot, level
The parameter ref.fn=c("median","mean") is not allowed because it is not possible to calculate the reference median array as parallelized.

\code{MAplotPara} is the parallelized function for MA plots of probe intensities. 
It is a function to check and control the Data quality of the samples using the MA plot method. 
For serial function an more details see \link[affy:AffyBatch-class]{boxplot}. This function is optimized for huge numbers of microarray data.

For using this function a computer cluster using the SNOW package has to be started. 
Starting the cluster with the command \code{makeCluster} generates an cluster object in the affyPara environment (.affyParaInternalEnv) and 
no cluster object in the global environment. The cluster object in the affyPara environment will be used as default cluster object, 
therefore no more cluster object handling is required.  
The \code{makeXXXcluster} functions from the package SNOW can be used to create an cluster object in the global environment and 
to use it for the preprocessing functions.
}

\value{  
MAplotPara return a list with four elements: \code{values_MAP}, \code{loess_y}, \code{quality_MAP}, \code{results_MAP}.
\item{values_MAP}{is a matrix which contains the calculated values for all arrays of the object (\code{sampleNames}, \code{S}, \code{osc_Loess}, \code{sigma}, \code{var_sigma}).}
\item{loess_y}{contains the y values (50 points) of the loess curve, which are used to calculate the \code{osc_Loess} from \code{values_MAP}}
\item{quality_MAP}{list which contains all "bad" quality Arrays}
\item{results_MAP}{summary from \code{quality_MAP} as matrix, which contains only the Arrays that are considered as "bad" quality and in which levels are they classified. Possible values are 0 if the Array is not at this levels and 1 if it is classified as "bad" sample at this level.}
} 

\author{Esmeralda Vicedo <e.vicedo@gmx.net>, Markus Schmidberger \email{schmidb@ibe.med.uni-muenchen.de}, Ulrich Mansmann \email{mansmann@ibe.med.uni-muenchen.de} }

\examples{
\dontrun{
library(affyPara)
if (require(affydata)) {
  data(Dilution)

  makeCluster(3)

  ##MA of Dilution data (affybatch)
  ma1 <- MAplotPara(Dilution)

  ## MAplot to a pdf file
  pdf(file="maplot.pdf", title="AffyBatch MAplot")
  ma2 <- MAplotPara(Dilution)
  dev.off()

  stopCluster()
}
}
}

\keyword{programming}
\keyword{manip}