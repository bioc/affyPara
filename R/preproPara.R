# file: preProPara.R
#
# Parallelization of the Preprocessing Function
#
# History
# 28.10.2008 : ... old stuff removed ...
# 16.05.2008 : Version 0.17 - one node bug fix
# 28.05.2008 : Version 0.18 - loess normalization added
# 23.10.2008 : Version 0.19 - awfull bug in checks removed
# 28.10.2008 : Version 0.20 - doSummarizationPara imporved
# 18.12.2008 : Version 0.21 - cluster object gets default parameter: .affyParaInternalEnv$cl
# 23.03.2009 : Version 0.22 - Option verbose set to getOption("verbose") and added . to names of internatl functions
# 24.03.2009 : Version 0.23 - Summarization optimized
#
# Copyright (C) 2008 : Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de>
###############################################################################

preproPara <- function(object,
		bgcorrect=TRUE, bgcorrect.method=NULL, bgcorrect.param=list(), 
		normalize=TRUE, normalize.method=NULL, normalize.param=list(), 
		pmcorrect.method=NULL, pmcorrect.param=list(),
		summary.method=NULL, summary.param=list(), ids=NULL,
		phenoData = new("AnnotatedDataFrame"), cdfname = NULL,
		cluster, verbose=getOption("verbose")) 
{
	#################
	# Check Functions
	#################
	#Check for affy amd snow
	require(affy)
	require(snow)
	
	#Get cluster object form default environment
	if(missing(cluster))
		cluster <- .affyParaInternalEnv$cl
	
	#Check cluster and generate number.parts
	checkCluster(cluster)
	number.parts <- length(cluster)
	
    #Check Methods
    if(bgcorrect){
    	if (is.null(bgcorrect.method))
    		stop("You have to choose a BG-Method")
    	if (any(bgcorrect.method == bgcorrect.methods()) == 0)
			stop(paste("Unknown BG-method (cannot find function '", bgcorrect.method,"')",sep=""))
    } else {
    	bgcorrect.method="none"
    }
    if(normalize){
    	if (is.null(normalize.method))
    		stop("You have to choose a Normalization-Method")
    	if (any(c("quantiles", "constant", "invariantset", "loess", "none")==normalize.method) == 0)
    		stop("Unknown Normalize-Method")
    } else {
    	normalize.method="none"
    }
   	if (is.null(summary.method))
    	stop("You have to choose a Summarization-Method")
    if (any( express.summary.stat.methods()==summary.method) == 0)
    	stop("Unknown Summarization-Method")

	#Check object type
	object.type <- .getObjectType(object) 
	
	#Check size of partitions
	parts <- .checkPartSize(object, number.parts)
	number.parts <- parts$number.parts
	object.length <- parts$object.length
	
	####################
	#Partition of object
	####################
	if (verbose) cat("Partition of object ")
		t0 <- proc.time();
		if (object.type == "AffyBatch"){
			object.list <- splitAffyBatch(object, number.parts)
			samples.names <- sampleNames(object)
		} else if( object.type == "CELfileVec" ){
			object.list <- splitFileVector(object, number.parts)
			samples.names <- gsub("^/?([^/]*/)*", "", unlist(object), extended = TRUE)
		} else if( object.type == "partCELfileList" ){
			object.list <- object
			object <- unlist(object)
			object.length <- length(object)
			samples.names <- gsub("^/?([^/]*/)*", "", unlist(object), extended = TRUE)
		}		
		if (normalize.method == "loess"){
			# Check for minimum number of arrays at each node
			if ( any( lapply(object.list, function(x) if (length(x) < 2) TRUE else FALSE  ) == TRUE ) ){
				cat("Object Distribution:", unlist(lapply(object.list,length))) 
				stop("There have to be minimum TWO arrays at each node!")
			}
		}		
		t1 <- proc.time();
	if (verbose) cat(paste(round(t1[3]-t0[3],3),"sec DONE\n"))				
	
	#Info-Output for Distribution
	if (verbose){ cat("\tObject Distribution: "); cat(paste(lapply(object.list,length))); cat("\n") }
	
	#################################
	#Initialize AffyBatches at slaves
	#################################
	if (verbose) cat("Initialize AffyBatches at slaves ")
		t0 <- proc.time();
		check <- clusterApply(cluster, object.list, .initAffyBatchSF, object.type) 
		t1 <- proc.time();
	if (verbose) cat(round(t1[3]-t0[3],3),"sec DONE\n")
	
	#Check phenoData and create TMP AffyBatch
	if (verbose) cat("Create TMP AffyBatch ")
	t0 <- proc.time();
	if( object.type == "CELfileVec" || object.type == "partCELfileList" ){
		headdetails <- clusterApply(cluster, object.list, .ReadHeaderSF)[[1]]
		dim.intensity <- headdetails[[2]]
		ref.cdfName <- headdetails[[1]]
		if( dim(phenoData)[1] == 0 ){
			pData <- data.frame(sample = seq(1, length(samples.names)), row.names = samples.names)
			varMetadata <- data.frame(labelDescription = "arbitrary numbering", row.names = names(pData))
			phenoData <- new("AnnotatedDataFrame", data = pData, varMetadata = varMetadata)
		}
		if (is.null(cdfname))
			cdfname <- ref.cdfName
		#Trick: exprs Matrix mit nur einer Zeile wird initialisiert
		exprs <- matrix(data = NA, nrow=1, ncol=object.length)
		AffyBatch <- new("AffyBatch", cdfName = cdfname,
				exprs=exprs, phenoData = phenoData,
				annotation = cleancdfname(cdfname, addcdf = FALSE))	
	} else if( object.type == "AffyBatch" ){
		AffyBatch <- object
	}
	t1 <- proc.time();
	if (verbose) cat(round(t1[3]-t0[3],2),"sec DONE\n")
	
	################
    # BGC on Slaves
	################
	if (bgcorrect.method != "none"){
		if (verbose) cat("BGC on Slaves ")
			t0 <- proc.time();
			check <- clusterCall(cluster, bgCorrectParaSF, bgcorrect.method)
			t1 <- proc.time();
		if (verbose) cat(round(t1[3]-t0[3],3),"sec DONE\n")
	}	
 
	##########################
	# Normalization on Slaves
	##########################
	if (normalize.method == "quantiles"){
		
		if (verbose) cat("Quantil Normalization on Slaves\n")
			t00 <- proc.time();
			if(is.null(normalize.param$type)){
				type <- "separate"
			} else type <- normalize.param$type
			############################
			# Do quantile normalization
			############################
			normalizeQuantilesPara(cluster,	type, object.length)
			t01 <- proc.time();
		if (verbose) cat(round(t01[3]-t00[3],3),"sec DONE\n")
		
	} else if (normalize.method == "constant"){
		
		if (verbose) cat("Constant Normalization on Slaves\n")
			t00 <- proc.time();
			if(is.null(normalize.param$refindex)){
				refindex <- 1
			} else refindex <- normalize.param$refindex
			if(is.null(normalize.param$FUN)){
				FUN <- mean
			} else FUN <- normalize.param$FUN
			if(is.null(normalize.param$na.rm)){
				na.rm <- TRUE
			} else na.rm <- normalize.param$na.rm
			############################
			# Do constant normalization
			############################
			normalizeConstantPara(cluster, samples.names, refindex=refindex, na.rm=na.rm, FUN=FUN)
			t01 <- proc.time();
		if (verbose) cat(round(t01[3]-t00[3],3),"sec DONE\n")
		
	} else if (normalize.method == "invariantset"){
		
		if (verbose) cat("Invariantset Normalization on Slaves ")
			t00 <- proc.time();
			if(is.null(normalize.param$prd.td)){
				prd.td <- c(0.003,0.007)
			} else prd.td <- normalize.param$prd.td	
			if(is.null(normalize.param$baseline.type)){
				baseline.type <- "mean"
			} else baseline.type <- normalize.param$baseline.type
			if(is.null(normalize.param$type)){
				type <- "separate"
			} else type <- normalize.param$type
			###############################
			# Do invariantset normalization
			###############################
			normalizeInvariantsetPara(cluster, AffyBatch, samples.names, prd.td=prd.td, baseline.type=baseline.type, type=type, verbose=verbose)
			t01 <- proc.time();
		if (verbose) cat(round(t01[3]-t00[3],3),"sec DONE\n")
		
	} else if (normalize.method == "loess"){
		
			if (verbose) cat("Loess Normalization on Slaves ")
			t00 <- proc.time();
			if(is.null(normalize.param$subset)){
				subset <- NULL
			} else subset <- normalize.param$subset
			if(is.null(normalize.param$type)){
				type <- "separate"
			} else type <- normalize.param$type
			if(is.null(normalize.param$epsilon)){
				epsilon <- 10^-2
			} else epsilon <- normalize.param$epsilon	
			if(is.null(normalize.param$maxit)){
				maxit <- 1
			} else maxit <- normalize.param$maxit
			if(is.null(normalize.param$span)){
				span <- 2/3
			} else span <- normalize.param$span
			if(is.null(normalize.param$family.loess)){
				family.loess <- "symmetric"
			} else family.loess <- normalize.param$family.loess
			if(is.null(normalize.param$log.it)){
				log.it <- TRUE
			} else log.it <- normalize.param$log.it
			###############################
			# Do loess normalization
			###############################
			#Check for order index at slaves
			check <- clusterCall(cluster, Sys.setlocale, "LC_COLLATE", Sys.getlocale("LC_COLLATE"))
			if (type == "pmonly"){
				if (verbose) cat("PM loess normalization\n")
			}
			if(type == "mmonly"){
				if (verbose) cat("MM loess normalization\n")
			}
			if (type == "together"){
				if (verbose) cat("PM and MM loess normalization\n")
			}
			if(type == "separate"){
				if (verbose) cat("PM and MM loess separate normalization\n")
				type <- "pmonly"
				normalizeLoessPara(cluster, samples.names, type, subset, epsilon, maxit, span, family.loess, log.it, object.length, verbose)
				type <- "mmonly"
			}
			normalizeLoessPara(cluster, samples.names, type, subset, epsilon, maxit, span, family.loess, log.it, object.length, verbose)
			t01 <- proc.time();
			if (verbose) cat(round(t01[3]-t00[3],3),"sec DONE\n")
		}

	
	#################
	#Do Summarization
	#################	
	if (verbose) cat("Summarization ")
	t0=proc.time();
	eset <- .doSummarizationPara(cluster, object.length, AffyBatch, 
			samples.names, ids=ids, pmcorrect.method=pmcorrect.method, summary.method=summary.method,
			summary.param=summary.param, pmcorrect.param=pmcorrect.param)
	t1=proc.time();
	if (verbose) cat(round(t1[3]-t0[3],3),"sec DONE\n")
		
	#Return Expression Set
	return(eset)
}