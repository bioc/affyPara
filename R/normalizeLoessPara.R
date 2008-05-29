# normalizeLoessPara.R
#
# Parallelization of the normalize.AffyBatch.loess function
#
# History
# 16.05.2008 : Version 0.1 - first ideas
# 26.05.2008 : Version 0.2 - working at one node
# 27.05.2008 : Version 0.3 - working on all nodes
# 28.05.2008 : Version 0.4 - working for all types
#
# Sending AffyBatch form master to slave an back is very time consuming. Sending a list
# of CEL files from master to slave, creating the AffyBatch and do normalization is faster.
# Using the right combination "size of AffyBatch on slaves" - "number of slaves" the parallelized
# version is more than ten times faster as the serial version.
#
# Copyright (C) 2008 : Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de>
###############################################################################

normalizeAffyBatchLoessPara <- function(cluster,
	object,
	phenoData = new("AnnotatedDataFrame"), cdfname = NULL,
	type=c("separate","pmonly","mmonly","together"), 
	subset = NULL,
	epsilon = 10^-2, maxit = 1, log.it = TRUE, 
	span = 2/3, family.loess ="symmetric",
	verbose=FALSE) 
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
	
	#Check arguments
	type <- match.arg(type)
	
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
		# Check for minimum number of arrays at each node
		if ( any( lapply(object.list, function(x) if (length(x) < 2) TRUE else FALSE  ) == TRUE ) ){
			cat("Object Distribution:", unlist(lapply(object.list,length))) 
			stop("There have to be minimum TWO arrays at each node!")
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
		check <- clusterApply(cluster, object.list, initAffyBatchSF, object.type) 
		t1 <- proc.time();
	if (verbose) cat(round(t1[3]-t0[3],3),"sec DONE\n")
	
	#################################
	#Normalization depending on types
	#################################
	t0 <- proc.time();
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
	t1 <- proc.time();
	if (verbose) cat(round(t1[3]-t0[3],3),"sec DONE\n")
	
	##############################
	#Combine / Rebuild affyBatches
	##############################
	if (verbose) cat("Rebuild AffyBatch ")
	t0 <- proc.time();
		AffyBatch.list.norm <- clusterCall(cluster, getAffyBatchSF)
		AffyBatch <- mergeAffyBatches(AffyBatch.list.norm)
	t1 <- proc.time();
	if (verbose) cat(round(t1[3]-t0[3],3),"sec DONE\n")
	
	#################
	#Return AffyBatch
	#################
	return(AffyBatch[,samples.names])
}

###
# Normalization Function for normalizeAffyBatchLoessPara
###
normalizeLoessPara <- function(cluster,
		samples.names, type, subset,
		epsilon, maxit, 
		span, family.loess,
		log.it, object.length,
		verbose=FALSE)
{
	if(verbose) cat("\tGenerate matrices at slaves","\n")
	#Generate matrices for loess normalization at slaves
	dimMat <- clusterCall(cluster, generateMatrizenSF, type, log.it)
	dimMat <- dimMat[!is.na(dimMat)][[1]]
	
	#Generate subset
	if (is.null(subset))
		subset <- sample(1:dimMat$nrow, min(c(5000,dimMat$nrow)))
	
	change <- epsilon + 1
	iter <- 0
	
	while(iter < maxit){
		iter <- iter + 1
		if(verbose) cat("\tIteration ",iter,"\n")
		
		#Do normalization at nodes
		if(verbose) cat("\t\tNormalization at slaves","\n")
		diag <- clusterCall(cluster, normalizeLoessParaSFnodes, subset, span, family.loess, iter, object.length)
		
		#Do normalization between nodes
		if (length(cluster)>1){
			#Matrix for allready normalized pairs
			sampleMatComp <- diag(length(samples.names))
			dimnames(sampleMatComp) <- list(samples.names,samples.names)
			for (i in 1:dim(sampleMatComp)[1]){
				for(j in 1:dim(sampleMatComp)[2])
					if (i>j) sampleMatComp[i,j]<-1
			}
			for (i in 1:length(diag))
				sampleMatComp[rownames(diag[[i]]),colnames(diag[[i]])] <- diag[[i]]
	
			if(verbose) cat("\t\tNormalization between slaves","\n")
			for (j in samples.names) {
				toNorm <- names(which(sampleMatComp[j,]==0))
				if (length(toNorm)!=0){
					if(verbose) cat("\t\tfor array: ",j,"\n")
					arrayInt <- clusterCall(cluster, getArraySF, j)	
					arrayInt<- unlist(arrayInt[lapply(arrayInt, length)>0])
					aux.list <- clusterCall(cluster, normalizeLoessParaSFbetNodes, arrayInt, j, toNorm, subset, span, family.loess, object.length)
					aux<-aux.list[[1]]
					for (i in 2:length(aux.list))
						aux <- aux + aux.list[[i]]
					#aux rewrite
					check<- clusterCall(cluster, writeMeansArraySF ,j, aux)
				}
			}
		}
		
		#Calculate change
		change <- clusterCall(cluster, writeArraySF, subset)
		change <- max(unlist(change))
			
		if(verbose) cat("\tChange: ",change,"\n")
	}
	
	#Reassign affyBatch at Slaves out of matrices form loess normalization
	check <- clusterCall(cluster, reassignMatrizenSF, log.it)
	
	if ((change > epsilon) & (maxit > 1))
		warning(paste("No convergence after", maxit, "iterations.\n"))
	
}

###
# Generate matrices for loess normalization
###
generateMatrizenSF  <- function(type, log.it)
{
	if (exists("AffyBatch", envir = .GlobalEnv)) {
		
		require(affy)
		
		#get AffyBatch and Intensity-Matirx depending on type
		AffyBatch <- get("AffyBatch", envir = .GlobalEnv)
		if (type == "pmonly"){
			Index <- unlist(indexProbes(AffyBatch,"pm"))
		}
		if(type == "mmonly"){
			Index <- unlist(indexProbes(AffyBatch,"mm"))
		}
		if (type == "together"){
			Index <- unlist(indexProbes(AffyBatch,"both"))
		}
		if(type == "separate"){
			#Nothing to do, see funktion above
		}
		mat <- intensity(AffyBatch)[Index,]
		
		J <- dim(mat)[2]
		II <- dim(mat)[1]
		newData <- mat
		if(log.it){
			mat <- log2(mat)
			newData <- log2(newData)
		}
		
		sampelMat <- diag(J)
		dimnames(sampelMat) <- list(sampleNames(AffyBatch),sampleNames(AffyBatch))
		for (i in 1:dim(sampelMat)[1]){
			for(j in 1:dim(sampelMat)[2]){
				if (i>j) sampelMat[i,j]<-1
			}
		}
		fs <- matrix(0, II, J)
		
		#Save matrices
		assign("Index", value=Index, envir= .GlobalEnv)
		assign("mat", value=mat, envir= .GlobalEnv)
		assign("newData", value=newData, envir= .GlobalEnv)
		assign("fs", value=fs, envir= .GlobalEnv)
		assign("sampelMat", value=sampelMat, envir= .GlobalEnv)
		
		return(list(nrow=II, ncol=J))
		
	} else
		return(NA)
	
}

###
# Loess Normalization at nodes
###
normalizeLoessParaSFnodes <- function(subset,
		span, family.loess, 
		iter, object.length)
{
	if ( exists("mat", envir = .GlobalEnv) &&
		 exists("newData", envir = .GlobalEnv) &&
		 exists("sampelMat", envir = .GlobalEnv) &&
		 exists("fs", envir = .GlobalEnv) ) {
		
		#Get parameters
		mat <- get("mat", envir = .GlobalEnv)
		J <- dim(mat)[2]
		II <- dim(mat)[1]
		
		sampelMat <- get("sampelMat", envir = .GlobalEnv)
		newData <- get("newData", envir = .GlobalEnv)
		fs <- get("fs", envir = .GlobalEnv)
		
		if (iter == 1){
			w <- c(0, rep(1,length(subset)), 0)
			assign("w", value=w, envir= .GlobalEnv)
		} else {
			w <- get("w", envir = .GlobalEnv)
		}
		
		means <- matrix(0,II,J) ##contains temp of what we substract
		sample_names <- colnames(newData)
		colnames(means) <- sample_names
		
		#Normalization
		for (j in 1:(J-1)){
			zeile <- sample_names[j]
			for (k in (j+1):J){
				spalte <-sample_names[k]
				sampelMat[zeile,spalte]<-1
				y <- newData[,j] - newData[,k]	
				x <- (newData[,j] + newData[,k]) / 2
				index <- c(order(x)[1], subset, order(-x)[1])
				xx <- x[index]
				yy <- y[index]
				aux <-loess(yy~xx, span=span, degree=1, weights=w, family=family.loess)
				aux <- predict(aux, data.frame(xx=x)) / object.length
				means[, j] <- means[, j] + aux
				means[, k] <- means[, k] - aux
			}
		}
		#Save results
		assign("means", value=means, envir= .GlobalEnv)
		
		return(sampelMat)
	} else
		return(NA)
}

###
# Get Intensity Array for loess normalization between nodes
###
getArraySF <- function(sampleName)
{
	if ( exists("newData", envir = .GlobalEnv) ) {
		
		#Get parameters
		newData <- get("newData", envir = .GlobalEnv)
		if ( any( colnames(newData)==sampleName ) )
			return( newData[,sampleName] )
	}else
		return(NA)	
}

###
# Loess Normalization between nodes
###
normalizeLoessParaSFbetNodes <- function(arrayInt, 
		sampleName, toNorm, 
		subset,	span, 
		family.loess, object.length)
{
	if ( exists("mat", envir = .GlobalEnv) &&
		 exists("newData", envir = .GlobalEnv) &&
		 exists("w", envir = .GlobalEnv) &&
		 exists("means", envir = .GlobalEnv) ){
		
		#Get parameters
		mat <- get("mat", envir = .GlobalEnv)
		J <- dim(mat)[2]
		II <- dim(mat)[1]
		
		newData <- get("newData", envir = .GlobalEnv)

		w <- get("w", envir = .GlobalEnv)
		means <- get("means", envir = .GlobalEnv)
		meansZw <- matrix(0,II,1) ##contains temp of what we substract
		
		for (k in toNorm){
			if (any(colnames(newData)==k)){
				y <- arrayInt - newData[,k]
				x <- (arrayInt + newData[,k]) / 2
				index <- c(order(x)[1], subset, order(-x)[1])
				xx <- x[index]
				yy <- y[index]
				aux <-loess(yy~xx, span=span, degree=1, weights=w, family=family.loess)
				aux <- predict(aux, data.frame(xx=x)) / object.length
				meansZw <- meansZw + aux
				means[, k] <- means[, k] - aux
			}
		}

		#Save results
		assign("means", value=means, envir= .GlobalEnv)
		
		return(meansZw)
		
	} else
		return(NA)
}

###
# Rewrite means array form loess normalization between arrays
###
writeMeansArraySF <- function(j, aux)
{
	if ( exists("means", envir = .GlobalEnv) ) {
		
		#Get parameters
		means <- get("means", envir = .GlobalEnv)
		
		if( any(colnames(means)==j) ){
			means[,j] <- means[,j] + aux
			#Save results
			assign("means", value=means, envir= .GlobalEnv)
			return(TRUE)
		}
		
	}else
		return(NA)	
}

###
# Rewrite Data / Matrices from loess normalization and calculate change
###
writeArraySF <- function(subset)
{
	if ( exists("mat", envir = .GlobalEnv) &&
		 exists("fs", envir = .GlobalEnv) &&
	     exists("means", envir = .GlobalEnv) ) {
		
		#Get parameters
		mat <- get("mat", envir = .GlobalEnv)
		fs <- get("fs", envir = .GlobalEnv)
		means <- get("means", envir = .GlobalEnv)
		
		fs <- fs + means
		newData <- mat - fs
		
		#Save results
		assign("fs", value=fs, envir= .GlobalEnv)
		assign("newData", value=newData, envir= .GlobalEnv)
		
		#Calculate change for iteration
		change <- max(colMeans((means[subset,])^2))
		return(change)
		
	}else
		return(NA)	
}

###
# Rebuild AffyBatch at nodes after loess normalization
###
reassignMatrizenSF <- function(log.it)
{
	if ( exists("AffyBatch", envir = .GlobalEnv) &&
			exists("newData", envir = .GlobalEnv) &&
			exists("Index", envir = .GlobalEnv)) {
		
		#Get matrices
		AffyBatch <- get("AffyBatch", envir = .GlobalEnv)
		newData <- get("newData", envir = .GlobalEnv)
		Index <- get("Index", envir = .GlobalEnv)
		
		# Rescale if necessary
		if(log.it) 
			newData <- 2^newData
		
		#Rewriting Data
		intensity(AffyBatch)[Index,] <- newData
		assign("AffyBatch", value=AffyBatch, envir= .GlobalEnv)
		
		return(TRUE)
		
	} else
		return(NA)
}
