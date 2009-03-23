# normalizeQuantilesPara.R
#
# Parallelization of the normalize.AffyBatch.quantiles function
#
# History
# 18.12.2008 : ... old stuff removed ...
# 22.02.2008 : Version 0.15 - modularization
# 26.02.2008 : Version 0.16 - error with ties fixed (affy 1.14 -> affy 1.16.2)
# 27.03.2008 : Version 0.17 - object.type as input removed
# 16.05.2008 : Version 0.18 - one node bug fix and code cleaning (tmp affyBatch removed)
# 18.12.2008 : Version 0.19 - cluster object gets default parameter: .affyParaInternalEnv$cl
# 23.03.2009 : Version 0.20 - Option verbose set to getOption("verbose") and added . to names of internatl functions
#
# Sending AffyBatch form master to slave an back is very time consuming. Sending a list
# of CEL files from master to slave, creating the AffyBatch and do normalization is faster.
# Using the right combination "size of AffyBatch on slaves" - "number of slaves" the parallelized
# version is more than ten times faster as the serial version.
#
# Copyright (C) 2008 : Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de>
###############################################################################

normalizeAffyBatchQuantilesPara <- function(object,
	phenoData = new("AnnotatedDataFrame"), cdfname = NULL,
	type=c("separate","pmonly","mmonly","together"), 
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
	
	#Check arguments
	type <- match.arg(type)
	
	#Check object type
	object.type <- .getObjectType(object) 
	
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
	if (verbose) cat(round(t1[3]-t0[3],3),"sec DONE\n")			
	
	#Info-Output for Distribution
	if (verbose) cat("Object Distribution:", unlist(lapply(object.list,length)), "\n") 
	
	#################################
	#Initialize AffyBatches at slaves
	##################################
	if (verbose) cat("Initialize AffyBatches at slaves ")
		t0 <- proc.time();
		check <- clusterApply(cluster, object.list, .initAffyBatchSF, object.type) 
		t1 <- proc.time();
	if (verbose) cat(round(t1[3]-t0[3],3),"sec DONE\n")
	
	#################################
	#Normalization depending on types
	#################################
	if (type == "pmonly"){
		if (verbose) cat("PM normalization\n")
	}
	if(type == "mmonly"){
		if (verbose) cat("MM normalization\n")
	}
	if (type == "together"){
		if (verbose) cat("PM and MM normalization\n")
	}
	if(type == "separate"){
		if (verbose) cat("PM and MM separate normalization\n")
	}
	normalizeQuantilesPara(cluster, type, object.length)
	
	##############################
	#Combine / Rebuild affyBatches
	##############################
	if (verbose) cat("Rebuild AffyBatch ")
	t0 <- proc.time();
		AffyBatch.list.norm <- clusterCall(cluster, .getAffyBatchSF)
		AffyBatch <- mergeAffyBatches(AffyBatch.list.norm)
	t1 <- proc.time();
	if (verbose) cat(round(t1[3]-t0[3],3),"sec DONE\n")
	
	
	#################
	#Return AffyBatch
	#################
	return(AffyBatch[,samples.names])
}

###
# Normalization Function for normalizeAffyBatchQuantilesPara
###
normalizeQuantilesPara <- function(cluster,
	type, object.length,
	verbose=getOption("verbose"))
{

	#Seperation for type separate
	separate = FALSE
	if (type == "separate"){
		if (verbose) cat("\tPM normalization\n")
		type = "pmonly"
		separate = TRUE
	}
	#Calculate rowMeans at parts	
	if (verbose) cat("\tCalculate row means ")
		t0 <- proc.time()	
	 	row_m <- clusterCall(cluster, normalizeQuantilesParaSF1, type)
		# remove list with length ==1 bzw. NAs
		row_m <- row_m[!unlist(lapply(row_m,length)==1)]
		#Calculate rowMean for all Parts
		rowMat <- do.call(cbind, lapply(row_m, get("/"), object.length) )	
	 	row_mean <- rowSums(rowMat)
	 	t1 <- proc.time()
	if (verbose) cat(paste(round(t1[3]-t0[3],3),"sec DONE\n"))
	
	#Normalize
	if (verbose) cat("\tQuantil Normalization ")
		t0 <- proc.time()
		check <- clusterCall(cluster, normalizeQuantilesParaSF2, row_mean)
	  	t1 <- proc.time()
	if (verbose) cat(paste(round(t1[3]-t0[3],3),"sec DONE\n"))
	
	#Second part for type separate
	if (separate==TRUE){
		type = "mmonly"
		if (verbose) cat("\tMM normalization\n")
		
		#Calculate rowMeans at parts	
		if (verbose) cat("\tCalculate row means ")
			t0 <- proc.time()	
			row_m <- clusterCall(cluster, normalizeQuantilesParaSF1, type)
			# remove list with length ==1 bzw. NAs
			row_m <- row_m[!unlist(lapply(row_m,length)==1)]
			#Calculate rowMean for all Parts
			rowMat <- do.call(cbind, lapply(row_m, get("/"), object.length) )
			row_mean <- rowSums(rowMat)
			t1 <- proc.time()
		if (verbose) cat(paste(round(t1[3]-t0[3],3),"sec DONE\n"))
		
		#Normalize
		if (verbose) cat("\tQuantil Normalization ")
			t0 <- proc.time()
			check <- clusterCall(cluster, normalizeQuantilesParaSF2, row_mean)
			t1 <- proc.time()
		if (verbose) cat(paste(round(t1[3]-t0[3],3),"sec DONE\n"))	
	}
}

###
# normalizeQuantilesPara -  Slavefunction 1
# calculates pratial rowMeans
###
normalizeQuantilesParaSF1 <- function(type)
{
	if (exists("AffyBatch", envir = .GlobalEnv)) {
		
		require(affy)
		
		#get AffyBatch and Intensity-Matirx
		AffyBatch <- get("AffyBatch", envir = .GlobalEnv)
		x <- intensity(AffyBatch)[,,drop=FALSE ]
			
		# Seperation by types
		if (type == "pmonly"){
			pms <- unlist(pmindex(AffyBatch))
			if (dim(x)[2]==1){
				noNA <- is.na(x[pms,]) == 0
			} else {
				noNA <- rowSums(is.na(x[pms,])) == 0
			}
			rows <- pms[noNA]
		}
		if(type == "mmonly"){
			mms <- unlist(mmindex(AffyBatch))
			if (dim(x)[2]==1){
				noNA <- is.na(x[mms,]) == 0
			} else {
				noNA <- rowSums(is.na(x[mms,])) == 0
			}	
			rows <- mms[noNA]
		}
		if (type == "together"){
			rows <- unlist(indexProbes(AffyBatch,"both"))
		}
		if(type == "separate"){
			#Nothing to calculate
		}
		
		# Temporary save rows at slaves
		assign("rows", value=rows, envir= .GlobalEnv)
		
		#Sort Matrix and return rowMeans
		cols <- dim(x)[2]
		if ( cols == 1 ){
			# no rowMeans to calculate
			return( sort(x[rows])*cols )
		} else if ( cols > 1 ){
			xsort <- apply(x[rows,], 2, sort)
			return( rowMeans(xsort)*cols )
		}
	} else
		return(NA)
}
	
###
# normalizeQuantilesPara -  Slavefunction 2
# normalize intensity Matrix
###
normalizeQuantilesParaSF2 <- function(row_mean)
{	
	if (exists("AffyBatch", envir = .GlobalEnv) && exists("rows", envir = .GlobalEnv)) {
		
		require(affy)
		
		#Get AffyBatch, IntensityMatrix and rows
		AffyBatch <- get("AffyBatch", envir = .GlobalEnv)
		x <- intensity(AffyBatch)[,,drop=FALSE ]
		rows <- get("rows", envir = .GlobalEnv)
		cols <- dim(x)[2]
		
		# Normalization
		if ( cols==1 ){
			xrank <- rank(x[rows])
			fix <- ( xrank-floor(xrank) > 0.4 )
			x[rows[fix]] <- 0.5 * (row_mean[floor(xrank[fix])] + row_mean[floor(xrank[fix])+1])
			x[rows[!fix]] <- row_mean[floor(xrank[!fix])]
			intensity(AffyBatch) <- x
		} else if ( cols>1 ){
			xrank <- apply(x[rows,],2,rank)
			for (j in 1:cols){
				fix <- ( xrank[,j]-floor(xrank[,j]) > 0.4 )
				x[rows[fix],j] <- 0.5 * (row_mean[floor(xrank[fix,j])] + row_mean[floor(xrank[fix,j])+1])
				x[rows[!fix],j] <- row_mean[floor(xrank[!fix,j])]
			}
			intensity(AffyBatch) <- x
		} else {
			return(NA)
		}
		assign("AffyBatch", value=AffyBatch, envir= .GlobalEnv)
		return(TRUE)
	} else
		return(NA)
}