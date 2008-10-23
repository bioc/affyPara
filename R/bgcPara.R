# file: bgPara.R
# 
# Parallelization of the bg.correct function
#
# History
# 08.04.2008 : ... old stuff removed ...
# 12.11.2007 : Version 0.10 - cleaning of slave functions
# 13.11.2007 : Version 0.11 - use of getAffyBatchSF
# 27.11.2007 : Version 0.12 - merge.AffyBatches renamed to mergeAffyBatches, bugfixes variable cluster
# 10.12.2007 : Version 0.13 - bugfixes method none!
# 27.03.2008 : Version 0.14 - object.type as input removed
# 08.04.2008 : Version 0.15 - some code cleaning
# 16.05.2008 : Version 0.16 - one node bug fix
# 08.09.2008 : Version 0.17 - remove all variables from all slaves added (von Esmeralda)
# 23.10.2008 : Version 0.18 - awfull bug in checks remuved
#
# Sending AffyBatch form master to slave an back is very time consuming. Sending a list
# of CEL files from master to slave, creating the AffyBatch and do BG-Correction is faster.
# Using the right combination "size of AffyBatch on slaves" - "number of slaves" the parallelized
# version is more than ten times faster as the serial version. 
#
# Copyright (C) 2008 : Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de>
###############################################################################

bgCorrectPara <- function(cluster,
	object,	phenoData = new("AnnotatedDataFrame"),
	method, verbose=FALSE)
{
	########
	# Checks
	########
	#Check for affy amd snow
	require(affy)
	require(snow)
	
	#Check cluster and generate number.parts
	checkCluster(cluster)
	number.parts <- length(cluster)
	
	#Check methods
	if ( any(bgcorrect.methods() == method) == 0 )
		 stop(paste("Unknown method (cannot find function '", method,"')",sep=""))
	 
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
			samples.names <- sampleNames(object)
			object.list <- splitAffyBatch(object, number.parts)
		} else if( object.type == "CELfileVec" ){
			object.list <- splitFileVector(object, number.parts)
		} else if( object.type == "partCELfileList" ){
			object.list <- object
		}				
		t1 <- proc.time();
	if (verbose) cat(paste(round(t1[3]-t0[3],3),"sec DONE\n"))				

	#Info-Output for Distribution
	if (verbose){ cat("Object Distribution: "); cat(paste(lapply(object.list,length))); cat("\n") }
	
	##################################
	# Initialize AffyBatches at slaves
	##################################
	if (verbose) cat("Initialize AffyBatches at slaves ")
		t0 <- proc.time();
		#remove all variables from all slaves
		check <- clusterCall(cluster, function(){rm(list=ls())})
	    check <- clusterApply(cluster, object.list, initAffyBatchSF, object.type) 
		t1 <- proc.time();
	if (verbose) cat(paste(round(t1[3]-t0[3],3),"sec DONE\n"))
	
	###########################
	#Do BG-Correction on slaves
	###########################
	if (verbose) cat("BGC on Slaves ")
		t0 <- proc.time();
		check <- clusterCall(cluster, bgCorrectParaSF, method)
		t1 <- proc.time();
	if (verbose) cat(paste(round(t1[3]-t0[3],3),"sec DONE\n"))
	
	##############################
	#Combine / Rebuild affyBatches
	##############################
	if (verbose) cat("Rebuild AffyBatch ")
		t0 <- proc.time();
		AffyBatch.list.bgc <- clusterCall(cluster, getAffyBatchSF)
		AffyBatch <- mergeAffyBatches(AffyBatch.list.bgc)
		t1 <- proc.time();
	if (verbose) cat(paste(round(t1[3]-t0[3],3),"sec DONE\n"))

	#Check phenoData
	if( object.type == "CELfileVec" || object.type == "partCELfileList" ){
		pdata <- pData(phenoData)
		if( object.type == "partCELfileList" )
			object <- unlist(object)
		if (dim(pdata)[1] != length(object)) {
			warning("Incompatible phenoData object. Created a new one.\n")
			samples.names <- gsub("^/?([^/]*/)*", "", unlist(object), extended = TRUE)
			pdata <- data.frame(sample = 1:length(object), row.names = samples.names)
			phenoData <- new("AnnotatedDataFrame", data = pdata, varMetadata = data.frame(labelDescription = "arbitrary numbering", row.names = "sample"))
		}
		else samples.names <- rownames(pdata)
		AffyBatch@phenoData <- phenoData
	}
		
	#################
	#Return AffyBatch
	#################
	return( AffyBatch[,samples.names] )
}

###
# Slavefunction for BG-Correction at Slaves
###
bgCorrectParaSF <- function(method) 
{
	if (exists("AffyBatch", envir = .GlobalEnv)) {
		require(affy)
		#load AffyBatch
		AffyBatch <- get("AffyBatch", envir = .GlobalEnv)
		#BGC AffyBatch
		AffyBatch <- bg.correct(AffyBatch, method=method)
		#temporary save AffyBatch
		assign("AffyBatch", value=AffyBatch, envir= .GlobalEnv)
		return(TRUE)
	} else
		return( NA )
}