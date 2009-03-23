# file: normalizeConstantPara.R
# 
# Parallelization of the normalize.AffyBatch.constant function
#
# History
# 27.03.2008 : ... old stuff removed ...
# 06.12.2007 : Version 0.5 - code cleaning and dokumentation
# 11.12.2007 : Version 0.6 - error in checks removed
# 22.02.2008 : Version 0.7 - modularization
# 27.03.2008 : Version 0.8 - object.type as input removed
# 16.05.2008 : Version 0.9 - one node bug fix
# 18.12.2008 : Version 0.10 - cluster object gets default parameter: .affyParaInternalEnv$cl
# 23.03.2009 : Version 0.11 - Option verbose set to getOption("verbose") and added . to names of internatl functions
#
# Sending AffyBatch form master to slave an back is very time consuming. Sending a list
# of CEL files from master to slave, creating the AffyBatch and do normalization is faster.
# Using the right combination "size of AffyBatch on slaves" - "number of slaves" the parallelized
# version is more than ten times faster as the serial version.
#
# Copyright (C) 2008 : Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de>
###############################################################################

normalizeAffyBatchConstantPara <- function(object,
		refindex=1, FUN=mean, na.rm=TRUE, 
		phenoData = new("AnnotatedDataFrame"), cdfname = NULL,
		cluster, verbose=getOption("verbose")) 
{
    #########
    # Checks
    #########
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
	
	#Check size of partitions
	parts <- .checkPartSize(object, number.parts)
	number.parts <- parts$number.parts
	object.length <- parts$object.length
	
	#Check refindex-Parameter
	if (! (refindex %in% 1:object.length))
		stop("invalid reference index for normalization")
	
	##################################
	#Partition of object and get data
	##################################
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
		object <- unlist(object.list)
		samples.names <- gsub("^/?([^/]*/)*", "", unlist(object), extended = TRUE)
	}				
	t1 <- proc.time();
	if (verbose) cat(paste(round(t1[3]-t0[3],3),"sec DONE\n"))				
	
	#Info-Output for Distribution
	if (verbose){ cat("Object Distribution: "); cat(paste(lapply(object.list,length))); cat("\n") }
	
	#################################
	#Initialize AffyBatches at slaves
	##################################
	if (verbose) cat("Initialize AffyBatches at slaves ")
		t0 <- proc.time();
		check <- clusterApply(cluster, object.list, .initAffyBatchSF, object.type) 
		t1 <- proc.time();
	if (verbose) cat(paste(round(t1[3]-t0[3],3),"sec DONE\n"))
	
	############################
	# Do constant normalization
	############################
	normalizeConstantPara(cluster, samples.names, refindex=refindex, na.rm=na.rm, FUN=FUN)
		
	##############################
	#Combine / Rebuild affyBatches
	##############################
	if (verbose) cat("Rebuild AffyBatch ")
		t0 <- proc.time();
		AffyBatch.list.norm <- clusterCall(cluster, .getAffyBatchSF)
		AffyBatch <- mergeAffyBatches(AffyBatch.list.norm)
		t1 <- proc.time();
	if (verbose) cat(paste(round(t1[3]-t0[3],3),"sec DONE\n"))
	
	#Return results
	return(AffyBatch[,samples.names])
}

###
# Function for constant normalization at slaves
###
normalizeConstantPara <- function(cluster,
		samples.names, refindex=1, 
		na.rm=TRUE, FUN=mean, 
		verbose=getOption("verbose"))
{
	####################################	
	#Create refconstant und refindexname
	####################################
	if (verbose) cat("Get refconstant ")
	t0 <- proc.time();
		refindexname <- samples.names[refindex]
		refconstantList <- clusterCall(cluster, normalizeConstantParaSF1, refindexname)	
		refconstant <- FUN(unlist(refconstantList[!unlist(lapply(lapply(refconstantList,is.na),any))]), na.rm=na.rm)
		
	t1 <- proc.time();
	if (verbose) cat(paste(round(t1[3]-t0[3],3),"sec DONE\n"))
	
	###########################
	#Do normalization on slaves
	###########################
	if (verbose) cat("Normalization ")
		t0 <- proc.time();
		check <- clusterCall(cluster, normalizeConstantParaSF2, refconstant=refconstant,  refindexname=refindexname, FUN=FUN, na.rm=na.rm)
		t1 <- proc.time();
	if (verbose) cat(paste(round(t1[3]-t0[3],3),"sec DONE\n"))
}

###
# Slavefunction 1
# to get intensities for Refindex
####
normalizeConstantParaSF1 <- function(refindexname)
{
	if (exists("AffyBatch", envir = .GlobalEnv)) {
		require(affy)
		#load AffyBatch
		AffyBatch <- get("AffyBatch", envir = .GlobalEnv)
		if( any(refindexname == sampleNames(AffyBatch)) ){
			return( intensity(AffyBatch)[,refindexname] )
		} else
			return(NA)
	} else
		return(NA)
}

###
# Slavefunction 2
# do normalization
####
normalizeConstantParaSF2 <- function(refconstant, refindexname, FUN = mean, na.rm = TRUE)
{
	if (exists("AffyBatch", envir = .GlobalEnv)) {
		require(affy)
		#load AffyBatch
		AffyBatch <- get("AffyBatch", envir = .GlobalEnv)
		
		n <- length(AffyBatch)
		normhisto <- vector("list", length=n)
		for (i in (1:n)) {
			#Do not normalize Ref-Array!
			if ( sampleNames(AffyBatch[,i]) != refindexname) {
				# normalize.constant can be used from library(affy)
				m <- normalize.constant(intensity(AffyBatch)[,i], refconstant, FUN=FUN, na.rm=na.rm)
				myhistory <- list(name="normalized by constant", constant=attr(m,"constant"))
				attr(m,"constant") <- NULL
				intensity(AffyBatch)[, i] <- m
				normhisto[[i]] <- myhistory
			}
		}
		attr(AffyBatch, "normalization") <- normhisto	
		assign("AffyBatch", value=AffyBatch, envir= .GlobalEnv)
		return(TRUE)
	} else
		return(NA)
}