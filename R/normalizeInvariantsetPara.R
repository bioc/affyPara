# file: normalizeInvariantsetPara.R
# 
# Parallelization of the normalize.AffyBatch.invariantset function
#
# History
# 27.03.2008 : ... old stuff removed ...
# 06.12.2007 : Version 0.5 - code cleaning and use of  splitFileVector
# 19.12.2007 : Version 0.6 - bugfixes problems with to more nodes than arrays
#							 error in checks removed, error in TMP affyBatch removed
# 28.02.2008 : Version 0.7 - modularization
# 27.03.2008 : Version 0.8 - object.type as input removed
# 16.05.2008 : Version 0.9 - one node bug fix
#
# Sending AffyBatch form master to slave an back is very time consuming. Sending a list
# of CEL files from master to slave, creating the AffyBatch and do normalization is faster.
# Using the right combination "size of AffyBatch on slaves" - "number of slaves" the parallelized
# version is more than ten times faster as the serial version.
#
# Copyright (C) 2008 : Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de>
###############################################################################

normalizeAffyBatchInvariantsetPara <- function(cluster,
		object,
		prd.td=c(0.003,0.007), 
		baseline.type=c("mean","median","pseudo-mean","pseudo-median"),
		type=c("separate","pmonly","mmonly","together"),
		phenoData = new("AnnotatedDataFrame"), cdfname = NULL,
		verbose=FALSE) 
{
    #########
    # Checks
    #########
	#Check for affy
	require(affy)
	require(snow)
	
	#Check cluster and generate number.parts
	checkCluster(cluster)
	number.parts <- length(cluster)
	
	#Check arguments
	type <- match.arg(type)
	baseline.type <- match.arg(baseline.type)

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
			object.length <- length(unlist(object))
			samples.names <- gsub("^/?([^/]*/)*", "", unlist(object), extended = TRUE)
		}				
		t1 <- proc.time();
	if (verbose) cat(paste(round(t1[3]-t0[3],3),"sec DONE\n"))				
	
	#Info-Output for Distribution
	if (verbose){ cat("\tObject Distribution: "); cat(paste(lapply(object.list,length))); cat("\n") }
	
	##################################
	#Initialize AffyBatches at slaves
	##################################
	if (verbose) cat("Initialize AffyBatches at slaves ")
		t0 <- proc.time();
		check <- clusterApply(cluster, object.list, initAffyBatchSF, object.type)
		t1 <- proc.time();
	if (verbose) cat(paste(round(t1[3]-t0[3],3),"sec DONE\n"))
	
	#########################################
	#Check phenoData and create TMP AffyBatch
	#########################################
	if (verbose) cat("Create TMP AffyBatch ")
	if( object.type == "CELfileVec" || object.type == "partCELfileList" ){
		t0 <- proc.time();
		headdetails <- clusterApply(cluster, object.list, ReadHeaderSF)[[1]]
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
		t1 <- proc.time();
	} else if( object.type == "AffyBatch" ){
		AffyBatch <- object
		dim <- dim(AffyBatch)[1]*dim(AffyBatch)[2]
	}
	if (verbose) cat(paste(round(t1[3]-t0[3],2),"sec DONE\n"))
	
	##############################
	# Normalization
	##############################
	normalizeInvariantsetPara(cluster, AffyBatch, samples.names, prd.td=prd.td, baseline.type=baseline.type, type=type, verbose=verbose)
	
	##############################
	#Combine / Rebuild affyBatches
	##############################
	if (verbose) cat("Rebuild AffyBatch ")
		t0 <- proc.time();
		AffyBatch.list.norm <- clusterCall(cluster,getAffyBatchSF)
		AffyBatch <- mergeAffyBatches(AffyBatch.list.norm)
		t1 <- proc.time();
	if (verbose) cat(paste(round(t1[3]-t0[3],3),"sec DONE\n"))
	
	#Return results
	return(AffyBatch[,samples.names])
}

##
# normalization Function
##
normalizeInvariantsetPara <- function(cluster,
		AffyBatch, samples.names,
		prd.td=c(0.003,0.007), 
		baseline.type=c("mean","median","pseudo-mean","pseudo-median"),
		type=c("separate","pmonly","mmonly","together"),
		verbose=TRUE)
{
	########################################
	#Parallel computation of means of arrays
	#and compute refindex
	#######################################
	if (verbose) cat("Calculate refindex ")
	t0 <- proc.time();
	if (baseline.type == "mean") {
		mvalues.list <- clusterCall(cluster, normalizeInvariantsetParaSF1, FUN=mean, type)
		mvalues.list <- mvalues.list[!unlist(lapply(mvalues.list,is.na))]
		refindex <- trunc(median(rank(unlist(mvalues.list))))
		refindexname <- samples.names[refindex]
	} else if (baseline.type == "median") {
		mvalues.list <- clusterCall(cluster, normalizeInvariantsetParaSF1, FUN=median, type)
		mvalues.list <- mvalues.list[!unlist(lapply(mvalues.list,is.na))]
		refindex <- trunc(median(rank(unlist(mvalues.list))))
		refindexname <- samples.names[refindex]
	} else if (baseline.type == "pseudo-mean" || baseline.type == "pseudo-median" ) {
		refindex <- 0
		refindexname <- 0
	}
	t1 <- proc.time();
	if (verbose) cat(paste(round(t1[3]-t0[3],2),"sec DONE\n"))
	
	#########################
	#Get Baseline.chip values
	########################
	if (verbose) cat("\tData from ")
	t0 <- proc.time();
	if (type == "pmonly") 
		rows <- unlist(pmindex(AffyBatch))
	else if (type == "mmonly")
		rows <- unlist(mmindex(AffyBatch))
	else if (type == "together")
		rows <- unlist(indexProbes(AffyBatch, "both"))
	else if (type == "separate")
		rows <- unlist(pmindex(AffyBatch))
	
	if (baseline.type == "mean" || baseline.type == "median") {
		baseline.chip.list <- clusterCall(cluster, getIntensitySF, rows, refindexname)
		baseline.chip <- unlist( baseline.chip.list[!unlist(lapply(lapply(baseline.chip.list,is.na),any))] )
	} else if (baseline.type == "pseudo-mean" ) {
		xpart <- clusterCall(cluster, getCompIntensitySF, rows)
		#Remove NAs
		listxpart <- unlist(xpart)
		if (any(is.na(listxpart))){
			omit <- seq_along(listxpart)[is.na(listxpart)]
			listxpart <- listxpart[-omit]
		}
		baseline.chip <- rowMeans(matrix(listxpart,ncol=length(AffyBatch)))
	} else if (baseline.type == "pseudo-median" ) {
		xpart <- clusterCall(cluster, getCompIntensitySF, rows)
		#Remove NAs
		listxpart <- unlist(xpart)
		if (any(is.na(listxpart))){
			omit <- seq_along(listxpart)[is.na(listxpart)]
			listxpart <- listxpart[-omit]
		}
		baseline.chip <- rowMedians(matrix(listxpart,ncol=length(AffyBatch)))
	}
	t1 <- proc.time();
	if (verbose) cat(paste(sampleNames(AffyBatch)[refindex], "used as baseline:", round(t1[3]-t0[3],2),"sec DONE\n"))
	
	###########################
	#Do normalization on slaves
	###########################
	if (verbose) cat(paste(type,"Normalization "))
	t0 <- proc.time();
	check <- clusterCall(cluster, normalizeInvariantsetParaSF2, refindexname, rows, prd.td, baseline.chip)
	t1 <- proc.time();
	if (verbose) cat(paste(round(t1[3]-t0[3],3),"sec DONE\n"))
	
	############################
	# second Part for seperate
	###########################
	if (type == "separate") {
		type <- "mmonly"
		########################################
		#Parallel computation of means of arrays
		#and Compute refindex
		#######################################
		if (verbose) cat("Calculate refindex ")
		t0 <- proc.time();
		if (baseline.type == "mean") {
			mvalues.list <- clusterCall(cluster, normalizeInvariantsetParaSF1, FUN=mean, type)
			mvalues.list <- mvalues.list[!unlist(lapply(mvalues.list,is.na))]
			refindex <- trunc(median(rank(unlist(mvalues.list))))
			refindexname <- samples.names[refindex]
		} else if (baseline.type == "median") {
			mvalues.list <- clusterCall(cluster, normalizeInvariantsetParaSF1, FUN=median, type)
			mvalues.list <- mvalues.list[!unlist(lapply(mvalues.list,is.na))]
			refindex <- trunc(median(rank(unlist(mvalues.list))))
			refindexname <- samples.names[refindex]
		} else if (baseline.type == "pseudo-mean" || baseline.type == "pseudo-median" ) {
			refindex <- 0
			refindexname <- 0
		}
		t1 <- proc.time();
		if (verbose) cat(paste(round(t1[3]-t0[3],2),"sec DONE\n"))
		if (verbose) cat("\tData from ")
		t0 <- proc.time();
		rows <- unlist(mmindex(AffyBatch))	
		if (baseline.type == "mean" || baseline.type == "median") {
			baseline.chip.list <- clusterCall(cluster, getIntensitySF, rows, refindexname)
			baseline.chip <- unlist( baseline.chip.list[!unlist(lapply(lapply(baseline.chip.list,is.na),any))] )
		} else if (baseline.type == "pseudo-mean" ) {
			xpart <- clusterCall(cluster, getCompIntensitySF, rows)
			#Remove NAs
			listxpart <- unlist(xpart)
			if (any(is.na(listxpart))){
				omit <- seq_along(listxpart)[is.na(listxpart)]
				listxpart <- listxpart[-omit]
			}
			baseline.chip <- rowMeans(matrix(listxpart,ncol=length(AffyBatch)))
		} else if (baseline.type == "pseudo-median" ) {
			xpart <- clusterCall(cluster, getCompIntensitySF, rows)
			#Remove NAs
			listxpart <- unlist(xpart)
			if (any(is.na(listxpart))){
				omit <- seq_along(listxpart)[is.na(listxpart)]
				listxpart <- listxpart[-omit]
			}
			baseline.chip <- rowMedians(matrix(listxpart,ncol=length(AffyBatch)))
		}
		t1 <- proc.time();
		if (verbose) cat(paste(sampleNames(AffyBatch)[refindex], "used as baseline:", round(t1[3]-t0[3],2),"sec DONE\n"))
		###########################
		#Do normalization on slaves
		###########################
		if (verbose) cat("separate2 Normalization ")
		t0 <- proc.time();
		check <- clusterCall(cluster, normalizeInvariantsetParaSF2, refindexname, rows, prd.td, baseline.chip)
		t1 <- proc.time();
		if (verbose) cat(paste(round(t1[3]-t0[3],3),"sec DONE\n"))
	}	
}

###
# Slavefunction 1
# calculate row means/medians
###
normalizeInvariantsetParaSF1 <- function(FUN, type)
{
	if (exists("AffyBatch", envir = .GlobalEnv)) {
		require(affy)
		#load AffyBatch
		AffyBatch <- get("AffyBatch", envir = .GlobalEnv)
		nc <- length(AffyBatch)
		if (type == "pmonly") 
			rows <- unlist(pmindex(AffyBatch))
		else if (type == "mmonly")
			rows <- unlist(mmindex(AffyBatch))
		else if (type == "together")
			rows <- unlist(indexProbes(AffyBatch, "both"))
		if (type == "separate") 
			rows <- unlist(pmindex(AffyBatch))
		m <- vector("numeric", length=nc)
		for (i in 1:nc)
			m[i] <- FUN(intensity(AffyBatch)[rows, i])
		return(m)
	} else
		return(NA)
}

###
# Slavefunction 2
# normalization
###
normalizeInvariantsetParaSF2 <- function(refindexname, rows, prd.td, baseline.chip)
{
	if (exists("AffyBatch", envir = .GlobalEnv)) {
		require(affy)
		#load AffyBatch
		AffyBatch <- get("AffyBatch", envir = .GlobalEnv)
		samples.names <- sampleNames(AffyBatch) 
		nc <- length(AffyBatch)
		normhisto <- vector("list", length=nc)
		for (i in (1:nc)) {
			if (sampleNames(AffyBatch[i]) != refindexname){
				#normalize.invariantset can be used from library(affy)
				tmp <- normalize.invariantset(c(intensity(AffyBatch)[rows,i]),c(baseline.chip),prd.td)
				tmp <- as.numeric(approx(tmp$n.curve$y, tmp$n.curve$x,xout=intensity(AffyBatch)[rows, i], rule=2)$y)
				attr(tmp,"invariant.set") <- NULL
				intensity(AffyBatch)[rows, i] <- tmp
			}
		}
		attr(AffyBatch, "normalization") <- normhisto
		assign("AffyBatch", value=AffyBatch[ , samples.names], envir= .GlobalEnv)
		return(TRUE)
	} else
		return(NA)
}