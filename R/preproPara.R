# file: preProPara.R
#
# Parallelization of the Preprocessing Function
#
# History
# 01.03.2008 : ... old stuff removed ...
# 15.02.2008 : Version 0.11 - bugfixes in summarization
# 18.02.2008 : Version 0.12 - bugfixes in constant normalization and check
# 19.02.2008 : Version 0.13 - progressbar removed (Namespace Problems)
# 22.02.2008 : Version 0.14 - modularization (first part)
# 28.02.2008 : Version 0.15 - modularization (second part)
# 27.03.2008 : Version 0.16 - object.type as input removed
# 16.05.2008 : Version 0.17 - one node bug fix
#
# Copyright (C) 2008 : Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de>
###############################################################################

preproPara <- function(cluster,
		object,
		bgcorrect=TRUE, bgcorrect.method=NULL, bgcorrect.param=list(), 
		normalize=TRUE, normalize.method=NULL, normalize.param=list(), 
		pmcorrect.method=NULL, pmcorrect.param=list(),
		summary.method=NULL, summary.param=list(), ids=NULL,
		phenoData = new("AnnotatedDataFrame"), cdfname = NULL,
		verbose=FALSE) 
{
	#################
	# Check Functions
	#################
	#Check for affy amd snow
	require(affy)
	require(snow)
	
	#Check cluster and generate number.parts
	checkCluster(cluster)
	number.parts <- length(cluster)
	
    #Check Methods
    if(bgcorrect){
    	if (is.null(bgcorrect.method))
    		stop("You have to choose a BG-Method")
    	if (any(bgcorrect.method == bgcorrect.methods) == 0)
			stop(paste("Unknown BG-method (cannot find function '", bgcorrect.method,"')",sep=""))
    } else {
    	bgcorrect.method="none"
    }
    if(normalize){
    	if (is.null(normalize.method))
    		stop("You have to choose a Normalization-Method")
    	if (any(c("quantiles", "constant", "invariantset", "none")==normalize.method) == 0)
    		stop("Unknown Normalize-Method")
    } else {
    	normalize.method="none"
    }
   	if (is.null(summary.method))
    	stop("You have to choose a Summarization-Method")
    if (any(generateExprSet.methods==summary.method) == 0)
    	stop("Unknown Summarization-Method")

	#Check object type
	object.type <- getObjectType(object) 
	
	#Check size of partitions
	parts <- checkPartSize(object, number.parts)
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
		t1 <- proc.time();
	if (verbose) cat(paste(round(t1[3]-t0[3],3),"sec DONE\n"))				
	
	#Info-Output for Distribution
	if (verbose){ cat("\tObject Distribution: "); cat(paste(lapply(object.list,length))); cat("\n") }
	
	#################################
	#Initialize AffyBatches at slaves
	#################################
	if (verbose) cat("Initialize AffyBatches at slaves ")
		t0 <- proc.time();
		check <- clusterApply(cluster, object.list, initAffyBatchSF, object.type) 
		t1 <- proc.time();
	if (verbose) cat(round(t1[3]-t0[3],3),"sec DONE\n")
	
	#Check phenoData and create TMP AffyBatch
	if (verbose) cat("Create TMP AffyBatch ")
	t0 <- proc.time();
	if( object.type == "CELfileVec" || object.type == "partCELfileList" ){
		headdetails <- clusterApply(cluster, object.list, ReadHeader)[[1]]
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
			normalizeQuantilesPara(cluster,	type, object.length, verbose=verbose)
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
			normalizeConstantPara(cluster, samples.names, refindex=refindex, na.rm=na.rm, FUN=FUN, verbose=verbose)
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
	}
	
	#################
	#Do Summarization
	#################
	eset <- doSummarizationPara(cluster, object.list, AffyBatch, 
			samples.names, ids=ids, pmcorrect.method=pmcorrect.method, summary.method=summary.method,
			summary.param=summary.param, pmcorrect.param=pmcorrect.param, verbose=verbose)
		
	#Return Expression Set
	return(eset)
}