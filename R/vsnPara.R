# vsnPara.R
#
# Parallelization of VSN
# we parallelize vsn2, because 'vsn' has been superseded by 'vsn2'.
#
# History
# 28.10.2008 : Version 0.1 - file vsnPara and vsnParaFunctions created
# 06.11.2008 : Version 0.2 - vsn reference implemented
# 18.12.2008 : Version 0.3 - cluster object gets default parameter: .affyParaInternalEnv$cl
# 23.03.2009 : Version 0.4 - Option verbose set to getOption("verbose") and added . to names of internatl functions
# 24.03.2009 : Version 0.5 - Summarization optimized
# 26.06.2009 : Version 0.6 - error in justvsnPara removed
# 08.03.2010 : Version 0.7 - gsub warning (extend=T) fixed
# 17.11.2010 : Version 0.8 - bug fix for cluster object
#
# Sending AffyBatch form master to slave an back is very time consuming. Sending a list
# of CEL files from master to slave, creating the AffyBatch and do normalization is faster.
# Using the right combination "size of AffyBatch on slaves" - "number of slaves" the parallelized
# version is more than ten times faster as the serial version.
#
# Copyright (C) 2009 : Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de>
###############################################################################

#######################################################
# vsn2 parallelized
#######################################################
vsn2Para <- function(object,
		phenoData = new("AnnotatedDataFrame"), cdfname = NULL,
		reference, subsample,
		..., 
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
		#samples.names <- gsub("^/?([^/]*/)*", "", unlist(object), extended = TRUE) #M.S. 8.3.2010 no more required
		samples.names <- gsub("^/?([^/]*/)*", "", unlist(object))
	} else if( object.type == "partCELfileList" ){
		object.list <- object
		object <- unlist(object)
		object.length <- length(object)
		#samples.names <- gsub("^/?([^/]*/)*", "", unlist(object), extended = TRUE) #M.S. 8.3.2010 no more required
		samples.names <- gsub("^/?([^/]*/)*", "", unlist(object))
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
	#initialize affybatch
	if(missing(reference))
		rm.list="ALL"
	else
		rm.list=c("wh", "whsel", "px")
	dimAB <- clusterApply(cluster, object.list, .initAffyBatchSF, object.type, rm.list=rm.list)
	dimAB <- .removeNA(dimAB)
	#initialize intensity matrix and remove AB
	check <- clusterCall(cluster, .setIntMatSF, rm.AB=FALSE)
	#calculate data distribution
	dist<-unlist(lapply(object.list,length))
	dist_list<-list()
	for(i in 1:length(dist)){
		if(i==1){
			dist_list[[1]]=1:dist[1]
			end <- dist_list[[i]][length(dist_list[[i]])]
		}else{
			dist_list[[i]]=(1+end):(end+dist[i])
			end <- dist_list[[i]][length(dist_list[[i]])]	
		}	
	}
	clusterApply(cluster, dist_list, function(x)assign("dist", value=x, envir= .GlobalEnv) )
	t1 <- proc.time();
	if (verbose) cat(round(t1[3]-t0[3],3),"sec DONE\n")
	
	#################################
	#Normalization depending on types
	#################################
	dimAB <- unlist(dimAB)
	ncol <- sum (dimAB[seq(2,length(dimAB),2)], na.rm=TRUE)
	nrow <- dimAB[1]
	dimAB <- c(nrow,ncol)
	
	#TODO wollen wir die 30000 wirklich?
	if(missing(subsample)){
		if( nrow > 30000L && missing(reference) )
			subsample = 30000L
		else
			subsample = 0L
	}
	fit <- vsnMatrixPara(cluster,
			dimAB, reference,  
			subsample=subsample,
			defaultpar = list(factr=5e7, pgtol=2e-4, maxit=60000L, trace=0L, cvg.niter=4L, cvg.eps=0), 
			...)
	
	#################
	#Return Parameters
	#################
	return(fit)
}

##############################################################################
# justvsn parallelized
##############################################################################
justvsnPara <- function(object,
		..., 
		cluster, verbose=getOption("verbose")) 
{
	#Get cluster object form default environment
	if(missing(cluster))
		cluster <- .affyParaInternalEnv$cl
	
	##############################
	# do vsn normalization
	##############################
	fit <- vsn2Para(object, cluster=cluster, ...) 
		
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
	return(AffyBatch)	
}

##############################################################################
# vsnrma parallelized
##############################################################################
vsnrmaPara <- function(object,
		pmcorrect.method="pmonly", pmcorrect.param=list(),
		summary.method="medianpolish", summary.param=list(),
		ids=NULL,
		phenoData = new("AnnotatedDataFrame"), cdfname = NULL,
		..., 
		cluster, verbose=getOption("verbose")) 
{
	########
	# Checks
	########
	
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
	
	#Check Parameter
	summary.method <- match.arg(summary.method, express.summary.stat.methods())
	pmcorrect.method <- match.arg(pmcorrect.method, pmcorrect.methods())
	
	####################
	#Partition of object
	###################
	if (object.type == "AffyBatch"){
		object.list <- splitAffyBatch(object, number.parts)
		samples.names <- sampleNames(object)
	} else if( object.type == "CELfileVec" ){
		object.list <- splitFileVector(object, number.parts)
		#samples.names <- gsub("^/?([^/]*/)*", "", unlist(object), extended = TRUE) #M.S. 8.3.2010 no more required
		samples.names <- gsub("^/?([^/]*/)*", "", unlist(object))
	} else if( object.type == "partCELfileList" ){
		object.list <- object
		object <- unlist(object)
		object.length <- length(object)
		#samples.names <- gsub("^/?([^/]*/)*", "", unlist(object), extended = TRUE) #M.S. 8.3.2010 no more required
		samples.names <- gsub("^/?([^/]*/)*", "", unlist(object))
	}				
	
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

	##############################
	# do vsn normalization
	##############################
	fit <- vsn2Para(object, ..., cluster=cluster, verbose=getOption("verbose")) 
	#remove log
	clusterCall(cluster, function(){
				if (exists("AffyBatch", envir = .GlobalEnv)) {
					require(affy)
					#load AffyBatch
					AffyBatch <- get("AffyBatch", envir = .GlobalEnv)
					intensity(AffyBatch) <- 2^intensity(AffyBatch)
					assign("AffyBatch", value=AffyBatch, envir= .GlobalEnv)
				} else
					return(NA)
			})
	
	#################
	#Do Summarization
	#################
	eset <- .doSummarizationPara(cluster, object.length, AffyBatch, 
			samples.names, ids=ids, pmcorrect.method=pmcorrect.method, summary.method=summary.method,
			summary.param=summary.param, pmcorrect.param=pmcorrect.param)
	
	#Return Expression Set
	return(eset)
	
}

#predict
#we use the predict function from vsn. useDataInFit=FALSE has to be set. There is no Data in the fit object stored.
