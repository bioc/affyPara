# normalizeLoessIterPara.R
#
# Randomaized and Parallelized new version of the normalize.AffyBatch.loess function
#
# History
# 23.06.2008 : Version 0.1 - first ideas and tests
# 23.06.2008 : Version 0.2 - working version
# 07.07.2008 : Version 0.3 - Bugfix for great Matrices
# 13.07.2008 : Version 0.4 - permutation Matrix added
# 14.07.2008 : Version 0.5 - output improved
# 27.08.2008 : Version 0.6 - maxPerm removed and percentPerm introduced
# 18.10.2008 : Version 0.7 - cluster - assign substituted by clusterExport
#
# Sending AffyBatch form master to slave an back is very time consuming. Sending a list
# of CEL files from master to slave, creating the AffyBatch and do normalization is faster.
# Using the right combination "size of AffyBatch on slaves" - "number of slaves" the parallelized
# version is more than ten times faster as the serial version.
#
# Copyright (C) 2008 : Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de>
###############################################################################

normalizeAffyBatchLoessIterPara <- function(cluster,
		object,
		percentPerm = 0.75,
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
		normalizeLoessIterPara(cluster, percentPerm, samples.names, type, subset, epsilon, maxit, span, family.loess, log.it, object.length, verbose)
		type <- "mmonly"
	}
	normalizeLoessIterPara(cluster, percentPerm, samples.names, type, subset, epsilon, maxit, span, family.loess, log.it, object.length, verbose)
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
normalizeLoessIterPara <- function(cluster,
		percentPerm = 0.75,
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
	#Matrix for normalization pairs
	n <- length(samples.names)
	nenner <- ( (n-1)*n / 2 )
	normPairMat <- matrix(NA,nrow=n,ncol=n)
	dimnames(normPairMat) <- list(samples.names,samples.names)
	for (i in 1:n){
		for(j in 1:n)
			if (i<j) normPairMat[i,j]<-0
  	}
  	check <- clusterCall(cluster, assign, "normPairMat", value=normPairMat, envir= .GlobalEnv) 
	
	#Generate subset
	if (is.null(subset))
		subset <- sample(1:dimMat$nrow, min(c(5000,dimMat$nrow)))
	
	#Iteration for Permutation
	percent <- 0
	iterPerm <- 0
	samples.names.akt <- samples.names
	while(percent < percentPerm){
		iterPerm <- iterPerm +1
		if(verbose) cat("\t", iterPerm, "Permutation of Arrays\n")
                                                                                                          
		if(iterPerm > 1){
			#Permutation of arrays
			samples.names.perm <- sample(samples.names.akt)
			samples.names.akt <- permMatrix(cluster, matName="mat", samples.names.akt, samples.names.perm, verbose = verbose)
			if ( all(samples.names.perm != samples.names.akt) ) stop("Error in Permutation")
		}
				
		# Iteration for loess Normalization
		change <- epsilon + 1
		iter <- 0
		while(iter < maxit){
			iter <- iter + 1
			if(verbose) cat("\t\tIteration ",iter,"\n")
			
			#Do normalization at nodes
			if(verbose) cat("\t\tNormalization at slaves","\n")
			normPairMat_list <- clusterCall(cluster, normalizeLoessIterParaSFnodes, subset, span, family.loess, iter, object.length)
			for (l in 1:length(normPairMat_list))
		    	normPairMat <- normPairMat +  normPairMat_list[[l]]
			
			zaehler <- table(normPairMat==0)[2]
			if( is.na(zaehler) )#all filled
				percent <- 1
			else
				percent <- 1 - ( zaehler / nenner )
					
			#Calculate change
			change <- clusterCall(cluster, writeArraySF, subset)
			change <- max(unlist(change))
			
			if(verbose) cat("\t\t\tChange: ",change,"\n\t\t\tNormalization: ",percent*100,"%\n")
		}
		
	}
	
	if (iterPerm > 1){
		#RePermutieren
		if(verbose) cat("\tRe-Permutation of Arrays\n")
		permMatrix(cluster, matName="mat", samples.names.akt, samples.names, verbose = verbose)
	}
	
	#Reassign affyBatch at Slaves out of matrices form loess normalization
	check <- clusterCall(cluster, reassignMatrizenSF, log.it)
	
	if ((change > epsilon) & (maxit > 1))
		warning(paste("No convergence after", maxit, "iterations.\n"))
	
}

###
# Loess Normalization at nodes
###
normalizeLoessIterParaSFnodes <- function(subset,
		span, family.loess, 
		iter, object.length)
{
	if ( exists("mat", envir = .GlobalEnv) &&
       exists("normPairMat", envir = .GlobalEnv)  ) {
		
		#Get parameters
		mat <- get("mat", envir = .GlobalEnv)
		normPairMat <- get("normPairMat", envir = .GlobalEnv)
		J <- dim(mat)[2]
		II <- dim(mat)[1]
		
		if (iter == 1){
			w <- c(0, rep(1,length(subset)), 0)
			assign("w", value=w, envir= .GlobalEnv)
		} else {
			w <- get("w", envir = .GlobalEnv)
		}
		
		means <- matrix(0,II,J) ##contains temp of what we substract
		sample_names <- colnames(mat)
		colnames(means) <- sample_names
		
		#Normalization
		for (j in 1:(J-1)){
			for (k in (j+1):J){
			  #Set in normPairMat
			  zeile <- sample_names[j]
			  spalte <- sample_names[k]
			  #only do normalization if allowd pair
        if( is.na(normPairMat[zeile,spalte]) ){ 
          normPairMat[spalte,zeile] <- normPairMat[spalte,zeile] + 1
          #loess normalization
  				y <- mat[,k] - mat[,j]	
  				x <- (mat[,k] + mat[,j]) / 2
  				index <- c(order(x)[1], subset, order(-x)[1])
  				xx <- x[index]
  				yy <- y[index]
  				aux <-loess(yy~xx, span=span, degree=1, weights=w, family=family.loess)
  				aux <- predict(aux, data.frame(xx=x)) / object.length
  				means[, k] <- means[, k] + aux
  				means[, j] <- means[, j] - aux
 				}		  
			  else if( !is.na(normPairMat[zeile,spalte]) ){ 
          normPairMat[zeile,spalte] <- normPairMat[zeile,spalte] + 1
          #loess normalization
  				y <- mat[,j] - mat[,k]	
  				x <- (mat[,j] + mat[,k]) / 2
  				index <- c(order(x)[1], subset, order(-x)[1])
  				xx <- x[index]
  				yy <- y[index]
  				aux <-loess(yy~xx, span=span, degree=1, weights=w, family=family.loess)
  				aux <- predict(aux, data.frame(xx=x)) / object.length
  				means[, j] <- means[, j] + aux
  				means[, k] <- means[, k] - aux
 				}
			}
		}
		#Save results
		assign("means", value=means, envir= .GlobalEnv)
		
		return(normPairMat)
	} else
		return(NA)
}