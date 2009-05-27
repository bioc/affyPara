# file: readAffybatchPara.R
#
# Parallelization of the read.affybatch function
#
# History
# 27.08.2009 : created - first idea
#
# Copyright (C) 2009 : Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de>
###############################################################################

read.affybatchPara <- function(object,
		phenoData = new("AnnotatedDataFrame"),
		description = NULL, notes = "",	
		cluster, verbose=getOption("verbose")) 
{
	########
	# Checks
	########
	#Check for affy amd snow
	require(affy)
	require(snow)
	
	#Get cluster object form default environment
	if(missing(cluster))
		cluster <- .affyParaInternalEnv$cl
	
	#Check cluster and generate number.parts
	checkCluster(cluster)
	number.parts <- length(cluster)
	
	#Check object type
	object.type <- .getObjectType(object) 
	if(object.type == "AffyBatch"){
		stop("This function does not work for an AffyBath!")
	}
	
	#Check size of partitions
	parts <- .checkPartSize(object, number.parts)
	number.parts <- parts$number.parts
	object.length <- parts$object.length
	
	####################
	#Partition of object
	###################
	if (verbose) cat("Partition of object ")
	t0 <- proc.time();
	if (object.type == "AffyBatch"){
		stop("This function does not work for an AffyBath!")
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
	if (verbose) cat(round(t1[3]-t0[3],3),"sec DONE\n")			
	
	#Info-Output for Distribution
	if (verbose) cat("Object Distribution:", unlist(lapply(object.list,length)), "\n") 
	
	#################################
	#Create AffyBatches at slaves
	##################################
	if (verbose) cat("Create AffyBatches at slaves ")
	t0 <- proc.time();
	check <- clusterApply(cluster, object.list, .initAffyBatchSF, object.type) 
	t1 <- proc.time();
	if (verbose) cat(round(t1[3]-t0[3],3),"sec DONE\n")
	
	##############################
	#Combine / Rebuild affyBatches
	##############################
	if (verbose) cat("Rebuild AffyBatch ")
	t0 <- proc.time();
	AffyBatch.list.norm <- clusterCall(cluster, .getAffyBatchSF)
	AffyBatch.list.norm <- .removeNA(AffyBatch.list.norm )
	AffyBatch <- mergeAffyBatches(AffyBatch.list.norm)
	t1 <- proc.time();
	if (verbose) cat(round(t1[3]-t0[3],3),"sec DONE\n")
	
	
	#################
	#Return AffyBatch
	#################
	phenoData(AffyBatch) <- phenoData
	#description(AffyBatch) <- description
	notes(AffyBatch) <- notes
	return(AffyBatch[,samples.names])
}
