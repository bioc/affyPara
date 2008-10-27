# computeExprSetPara.R
# 
# Parallelization of the computeExprSet function
#
# History
# 27.03.2008 : ... old stuff removed ...
# 15.02.2008 : Version 0.5 - code cleaning and dokumentation
# 19.02.2000 : Version 0.6 - progressbar removed (Namespace Problems)
# 22.02.2008 : Version 0.7 - modularization
# 27.03.2008 : Version 0.8 - object.type as input removed
# 16.05.2008 : Version 0.9 - one node bug fix
# 08.08.2008 : Version 0.10 - calculate expressionSet at nodes -> speed improvement
# 22.08.2008 : Version 0.11 - different summary.methods
# 23.10.2008 : Version 0.12 - awfull bug in checks remuved
#
# Copyright (C) 2008 : Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de>
###############################################################################

computeExprSetPara <- function(cluster,
		object,	ids=NULL,
		pmcorrect.method, summary.method,
		summary.param=list(), pmcorrect.param=list(),
		phenoData = new("AnnotatedDataFrame"), cdfname = NULL,
		verbose=TRUE) 
{
	#########
	# Checks
	#########
	#Check for affy amd snow
	require(affy)
	require(snow)
	
	#Check cluster and generate number.parts
	checkCluster(cluster)
	number.parts <- length(cluster)
	
	#Check object type
	object.type <- getObjectType(object) 
	
	#Check size of partitions
	parts <- checkPartSize(object, number.parts)
	number.parts <- parts$number.parts
	object.length <- parts$object.length
	
	#Check Parameter
	pmcorrect.method <- match.arg(pmcorrect.method, pmcorrect.methods())
	upDate.express.summary.stat.methods(c(express.summary.stat.methods(), 'medianpolish_orig', 'liwong_orig', 'farms_orig', 'playerout_orig'))
	summary.method <- match.arg(summary.method, express.summary.stat.methods())
	if (summary.method == "medianpolish")
		warning("Medianpolish only at nodes, not over complete AffyBatch!\nFor same results use 'medianpolish_orig' (slow)")
	if (summary.method == "liwong")
		warning("Liwong only at nodes, not over complete AffyBatch!\nFor same results use 'liwong_orig' (slow)")
	if (summary.method == "farms")
		warning("Farms only at nodes, not over complete AffyBatch!\nFor same results use 'farms_orig' (slow)")
	if (summary.method == "playerout")
		warning("Playerout only at nodes, not over complete AffyBatch!\nFor same results use 'playerout_orig' (slow)")
	
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
	#################################
	if (verbose) cat("Initialize AffyBatches at slaves ")
	t0 <- proc.time();
	check <- clusterApply(cluster, object.list, initAffyBatchSF, object.type)
	t1 <- proc.time();
	if (verbose) cat(paste(round(t1[3]-t0[3],3),"sec DONE\n"))
	
	#Check phenoData and create TMP AffyBatch
	if (verbose) cat("Create TMP AffyBatch ")
	t0 <- proc.time();
	if( object.type == "CELfileVec" || object.type == "partCELfileList" ){
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
	} else if( object.type == "AffyBatch" ){
		AffyBatch <- object
	}
	t1 <- proc.time();
	if (verbose) cat(paste(round(t1[3]-t0[3],2),"sec DONE\n"))
	
	#################
	#Do Summarization
	#################
	if (verbose) cat("Summarization ")
	t0=proc.time();
	if ( regexpr("_orig", summary.method) == -1 ){
		eset_list <- clusterApply(cluster, object.list, calculateExprSetSF, 
				ids=ids, pmcorrect.method=pmcorrect.method, summary.method=summary.method,
				summary.param=summary.param, pmcorrect.param=pmcorrect.param, verbose=verbose)
		
		#Generat Exprs-Matrix
		eset_mat <- exprs(eset_list[[1]])[ sort(rownames(exprs(eset_list[[1]]))) ,]
		
		if ( length(eset_list) > 1){
			for( k in 2:length(eset_list)){
				eset_mat <- cbind(eset_mat,exprs(eset_list[[k]])[ sort(rownames(exprs(eset_list[[k]]))) ,] )
				eset_list[[k]]<-NA
			}
		}
		if (is.null(ids)){
			ids<-geneNames(AffyBatch)
		}
		dimnames(eset_mat) <- list(ids, samples.names)
	
		eset <- new("ExpressionSet",
				phenoData=phenoData(AffyBatch),
				## featureData picked up from exprs
				experimentData=experimentData(AffyBatch),
				exprs=eset_mat,
				annotation=annotation(AffyBatch)
		)
		
	#playerout_orig farms_orig medianpolish_orig liwong_orig
	}else {
		summary.method <- substr(summary.method,1,regexpr("_orig", summary.method)-1)
		eset <- doSummarizationPara(cluster, object.list, AffyBatch, 
				samples.names, ids=ids, pmcorrect.method=pmcorrect.method, summary.method=summary.method,
				summary.param=summary.param, pmcorrect.param=pmcorrect.param, verbose=verbose)		
	}
	
	t1=proc.time();
	if (verbose) cat(round(t1[3]-t0[3],3),"sec DONE\n")
	
	#######################
	#Return Expression Set
	#######################
	return(eset)
}

###
# Call computeExprSet at slaves
###
calculateExprSetSF <- function( cluster, object.list,
		samples.names, ids=NULL,
		pmcorrect.method, summary.method,
		summary.param=list(), pmcorrect.param=list(),
		verbose=verbose )
{
	if ( exists("AffyBatch", envir = .GlobalEnv) ) {
		require(affy)
		#if (summary.method=="farms") library(farms)
		#load AffyBatch
		AffyBatch <- get("AffyBatch", envir = .GlobalEnv)
		eset <- computeExprSet(AffyBatch, pmcorrect.method=pmcorrect.method, pmcorrect.param=pmcorrect.param,
				summary.method=summary.method, summary.param=summary.param, ids=ids)
		return( eset )
	}else
		return(NA)
}

###
# Function for parallel summarization
###
doSummarizationPara <- function(cluster, 
		object.list, AffyBatch, 
		samples.names, ids=NULL,
		pmcorrect.method, summary.method,
		summary.param=list(), pmcorrect.param=list(),
		verbose=verbose)
{
	#################################
	#Get Parameters for Summarization
	##################################
	if (verbose) cat("Get Parameters ")
	t0=proc.time();
	## Anzahl von samples
	n <- sum(unlist(lapply(object.list,length)))
	
	## if 'ids' is NULL compute for all ids
	if (is.null(ids)){
		ids <-geneNames(AffyBatch)
	}
	
	m <- length(ids)
	
	pps.warnings <- vector("list", length=m)
	
	## cheap trick to (try to) save time
	c.pps <- new("ProbeSet", pm=matrix(), mm=matrix())
	
	## matrix to hold expression values
	exp.mat <- matrix(NA, m, n)
	se.mat <- matrix(NA, m, n)
	
	CDFINFO <- getCdfInfo(AffyBatch)
	
	t1=proc.time();
	if (verbose) cat(round(t1[3]-t0[3],3),"sec DONE\n")
	
	############################
	# Summarization
	############################
	t0=proc.time();
	if (verbose) 
		cat("Sumarization:",m, "ids to be processed ")
	
	## loop over the ids
	mycall <- as.call(c(getMethod("express.summary.stat",
							signature=c("ProbeSet","character", "character")),
					list(c.pps, pmcorrect=pmcorrect.method, summary=summary.method,
							summary.param=summary.param, pmcorrect.param=pmcorrect.param))
	)
	
	
	
	for (i in seq(along=ids)) {
		
		if (verbose && i%%100==0) {
			cat(".")
		}
		
		id <- ids[i]
		
		if (! exists(id, envir=CDFINFO)) {
			pps.warnings[[i]] <- paste("Unknown id", id)
		} else {
			## locations for an id
			loc <- get(id, envir=CDFINFO)
			l.pm <- loc[, 1]
			if (ncol(loc) == 2)
				l.mm <- loc[ ,2]
			else
				l.mm <- integer()
			
			zw <- clusterCall(cluster, getCompIntensityMatrixSF, rows=l.pm)
			zw <- zw[!is.na(zw)]
			c.pps@pm <- do.call(cbind, zw)
			zw <- clusterCall(cluster, getCompIntensityMatrixSF, rows=l.mm)
			zw <- zw[!is.na(zw)]
			c.pps@mm <- do.call(cbind, zw)
			
			## generate expression values
			## (wrapped in a sort of try/catch)
			mycall[[2]] <- c.pps
			ev <- try(eval(mycall), silent = TRUE)
		}
		if (! inherits(ev, "try-error")) {
			exp.mat[i, ] <- ev$exprs
			se.mat[i,] <- ev$se.exprs
		} else {
			pps.warnings[[i]] <- ev[1]
		}
		
	}
	
	t1=proc.time();
	if (verbose) cat(" ",round(t1[3]-t0[3],3),"sec DONE\n")
	
	######################
	#Build ExpressionSet
	######################
	if (verbose) cat("Build ExpressionSet ")
	t0=proc.time();
	dimnames(exp.mat) <- list(ids, samples.names)
	dimnames(se.mat) <- list(ids, samples.names)
	eset <- new("ExpressionSet",
			phenoData=phenoData(AffyBatch),
			## featureData picked up from exprs
			experimentData=experimentData(AffyBatch),
			exprs=exp.mat,
			se.exprs=se.mat,
			annotation=annotation(AffyBatch)
	)
	attr(eset, "pps.warnings") <- pps.warnings
	t1=proc.time();
	if (verbose) cat(paste(round(t1[3]-t0[3],3),"sec DONE\n"))
	
	return(eset)
	
}