# file: affyParaIntern.R
# 
# generall support SLAVE - functions for affyPara
#
# History
# 17.07.2008 : ... old stuff removed ...
# 28.05.2008 : Version 0.12 - never used functions removed
# 23.06.2008 : Version 0.13 - permMatrix added
# 23.06.2008 : Version 0.14 - error for Repermutation fixed
# 07.07.2008 : Version 0.15 - initAffyBatchSF, error for small numbers of nodes fixed
# 17.07.2008 : Version 0.16 - bugFix in setArraySF
# 11.09.2008 : Version 0.17 - file.name removed -> basename
# 15.09.2008 : Version 0.18 - output von initAffyBatchSF set to dimAB
# 
# Copyright (C) 2008 : Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de>
###############################################################################

###
# Permutation of Arrays (over Nodes)
###
permArrays <- function(cluster, sample.names, verbose=TRUE)
{
	t0 <- proc.time();
	# Permutation
	perm <- sample(sample.names)
	
	for (i in perm){
		if (verbose) cat("Move array: ", i, "\n")
		#get array from slaves
		array <- clusterCall(cluster, getArraySF, i)
		array <- array[[which(!is.na(array))]]
		# move array to new positon
		alter_name_neue_pos <- sample.names[which(perm==i)]
		check <- clusterCall(cluster, setArraySF, array, i, alter_name_neue_pos)
	}
	#move new matrix to affyBatch at slave
	check <- clusterCall(cluster, resetABSF)
	cat("New sample list: ", unlist(check)[!is.na(unlist(check))], "\n")

	t1 <- proc.time();
	if (verbose) cat(round(t1[3]-t0[3],3),"sec DONE\n")		
	
	return(perm)
}
#####
# save intensities of array for special sample-Name at slaves
#####
setArraySF <- function(colNEU, col_name_Neu, alter_name_neue_pos)
{
	
	if (exists("AffyBatch", envir = .GlobalEnv)) {
		
		require(affy)
		#load AffyBatch
		AffyBatch <- get("AffyBatch", envir = .GlobalEnv)
		
		#load or initialize tmp mat
		if ( exists("mat", envir = .GlobalEnv) && 
				exists("sample_NAMES_neu", envir = .GlobalEnv)	){
			mat <- get("mat", envir = .GlobalEnv)
			sample_NAMES_neu <- get("sample_NAMES_neu", envir = .GlobalEnv)
		} else {
			dimAB <- dim( exprs(AffyBatch) )
			mat <- matrix( ncol=dimAB[2], nrow=dimAB[1], dimnames=list(c(),sampleNames(AffyBatch)) )
			sample_NAMES_neu <- vector(length=dimAB[2])
			
		}
		
		#write new array intensities to right position
		if ( any(sampleNames(AffyBatch)==alter_name_neue_pos) ){
			mat[, alter_name_neue_pos] <- colNEU
			sample_NAMES_neu[which(sampleNames(AffyBatch)==alter_name_neue_pos)] <- col_name_Neu
		}
		
		#save data
		assign("mat", value=mat, envir= .GlobalEnv)
		assign("sample_NAMES_neu", value=sample_NAMES_neu, envir= .GlobalEnv)
		
		return(sample_NAMES_neu)
	} else
		return(NA)
}
#####
# get intensities of array for special sample name from slaves
#####
getArraySF <- function(sample_name)
{
	if (exists("AffyBatch", envir = .GlobalEnv)) {
		require(affy)
		#load AffyBatch
		AffyBatch <- get("AffyBatch", envir = .GlobalEnv)
		if ( any(sampleNames(AffyBatch)==sample_name) ){
			mat <- as.matrix( exprs(AffyBatch)[,sample_name] )
			colnames(mat) <- sample_name
			return(mat)	
		}
		else
			return(NA)
	} else
		return(NA)
}
####
# Rewrite tmp matrix to affyBatch
####
resetABSF <- function()
{
	if (exists("AffyBatch", envir = .GlobalEnv) &&
			exists("mat", envir = .GlobalEnv) && 
			exists("sample_NAMES_neu", envir = .GlobalEnv)	) {
		
		require(affy)
		#load data
		AffyBatch <- get("AffyBatch", envir = .GlobalEnv)
		mat <- get("mat", envir = .GlobalEnv)
		sample_NAMES_neu <- get("sample_NAMES_neu", envir = .GlobalEnv)
		
		#Rewrite data
		exprs(AffyBatch) <- mat
		sampleNames(AffyBatch) <- sample_NAMES_neu
		
		#remove data and save AffyBatch
		rm(mat,sample_NAMES_neu, envir= .GlobalEnv )
		assign("AffyBatch", value=AffyBatch, envir= .GlobalEnv)
		
		return(sampleNames(AffyBatch))
		
	} else
		return(NA)
}

#####################################################################################
###
# Permutation of Matrix by column (over Nodes)
###
permMatrix <- function(cluster, matName="mat", sample.names, sample.names.perm, verbose=TRUE)
{
	t0 <- proc.time();

	for (i in sample.names.perm){
		alter_name_neue_pos <- sample.names[which(sample.names.perm==i)]
		if (verbose) cat("\t\tMove column: ", alter_name_neue_pos, " -> ",i,"\n")
		#get column from slaves
		col <- clusterCall(cluster, getColSF, matName, i)
		col <- col[[which(!is.na(col))]]
		# move column to new positon
		check <- clusterCall(cluster, setColSF, matName, col, i, alter_name_neue_pos)
	}

	#write new matrix over old matrix at slave
	col.names <- clusterCall(cluster, resetMatSF, matName)

	t1 <- proc.time();
	if (verbose) cat("\t\t",round(t1[3]-t0[3],3),"sec DONE\n")		
	
	return(unlist(col.names))
}
#####
# get col-intensities of matix for special name from slaves
#####
getColSF <- function(matName="mat", col_name)
{
	if (exists(matName, envir = .GlobalEnv)) {
		#load mat
		mat <- get(matName, envir = .GlobalEnv)
		if ( any(colnames(mat)==col_name) ){
			col <- as.matrix( mat[,col_name] )
			colnames(col) <- col_name
			return(col)	
		}
		else
			return(NA)
	} else
		return(NA)
}
#####
# save col-intensities of matrix for special col-Name at slaves
#####
setColSF <- function(matName = "mat", colNEU, col_name_Neu, alter_name_neue_pos)
{
	if (exists(matName, envir = .GlobalEnv)) {
		#load mat
		mat <- get(matName, envir = .GlobalEnv)
		
		if ( any(colnames(mat)==alter_name_neue_pos) ){
		
			#load or initialize tmp mat
			if ( exists("matZw", envir = .GlobalEnv) && 
				 exists("col_NAMES_neu", envir = .GlobalEnv)	){
					matZw <- get("matZw", envir = .GlobalEnv)
					col_NAMES_neu <- get("col_NAMES_neu", envir = .GlobalEnv)
			} else {
				dimMat <- dim( mat )
				matZw <- matrix( ncol=dimMat[2], nrow=dimMat[1], dimnames=list(c(),colnames(mat)) )
				col_NAMES_neu <- vector(length=dimMat[2])		
			}		
			#write new mat intensities to correct position
			matZw[, alter_name_neue_pos] <- colNEU
			col_NAMES_neu[which(colnames(mat)==alter_name_neue_pos)] <- col_name_Neu
		
			#save data
			assign("matZw", value=matZw, envir= .GlobalEnv)
			assign("col_NAMES_neu", value=col_NAMES_neu, envir= .GlobalEnv)
			
			return(TRUE)
		} else
			return(NA)
	} else
		return(NA)
}
####
# Rewrite tmp matrix
####
resetMatSF <- function(matName="mat")
{
	if ( exists(matName, envir = .GlobalEnv) &&
		 exists("matZw", envir = .GlobalEnv) && 
		 exists("col_NAMES_neu", envir = .GlobalEnv)	) {
		
		#load data
		mat <- get(matName, envir = .GlobalEnv)
		matZw <- get("matZw", envir = .GlobalEnv)
		col_NAMES_neu <- get("col_NAMES_neu", envir = .GlobalEnv)
		
		#Rewrite data
		mat <- matZw
		colnames(mat) <- col_NAMES_neu
		
		#remove data and save matrix
		rm(matZw,col_NAMES_neu, envir= .GlobalEnv)
		assign(matName, value=mat, envir= .GlobalEnv)
		
		return(colnames(mat))
		
	} else
		return(NA)
}

##########################################################################
###
# Initializing AffyBatch at Slaves
###
initAffyBatchSF <- function(object, object.type)
{
	require(affy)
	#remove old AffyBatches
	if (exists("AffyBatch", envir = .GlobalEnv))
		rm(AffyBatch, envir = .GlobalEnv)
	
	#create AffyBatch
	if (object.type == "AffyBatch")
		AffyBatch <- object
	else if( length(object) != 0 )
		AffyBatch <- ReadAffy(filenames=object)
	else
		return(FALSE)

	#temporary save AffyBatch
	assign("AffyBatch", value=AffyBatch, envir= .GlobalEnv)
	return(dim(exprs(AffyBatch)))
}

########################################
###
# Get AffyBatches from Slaves
###
getAffyBatchSF <- function()
{
	if (exists("AffyBatch", envir = .GlobalEnv)) {
		require(affy)
		#load and return AffyBatch
		return( get("AffyBatch", envir = .GlobalEnv) )
	}else
		return(NA)
}

###
# Excecution of a function for AffyBatch from slaves
###
getFUNAffyBatchSF <- function(FUN)
{
	if (exists("AffyBatch", envir = .GlobalEnv)) {
		require(affy)
		#load AffyBatch
		AffyBatch <- get("AffyBatch", envir = .GlobalEnv)
		return( FUN(AffyBatch) )
	} else
		return(NA)
}

###
# Get special values from intensity matrix from slaves
###
getIntensitySF <- function(rows, refindexname)
{
	if (exists("AffyBatch", envir = .GlobalEnv)) {
		require(affy)
		#load AffyBatch
		AffyBatch <- get("AffyBatch", envir = .GlobalEnv)
		if( any(sampleNames(AffyBatch) == refindexname) )
			return( c(intensity(AffyBatch)[rows, refindexname]) )
		else
			return(NA)
	} else
		return(NA)
}

###
#  Get special rows from intensity matrix from slaves
###
getCompIntensitySF <- function(rows)
{
	if (exists("AffyBatch", envir = .GlobalEnv)) {
		require(affy)
		#load AffyBatch
		AffyBatch <- get("AffyBatch", envir = .GlobalEnv)
		return( c(intensity(AffyBatch)[rows, ]) )
	} else
		return(NA)
}

###
# Get complete intensity matrix from all slaves
###
getCompIntensityMatrixSF <- function(rows, drop=FALSE)
{
	if (exists("AffyBatch", envir = .GlobalEnv)) {
		require(affy)
		#load AffyBatch
		AffyBatch <- get("AffyBatch", envir = .GlobalEnv)
		return( intensity(AffyBatch)[rows, , drop=drop] )
	} else
		return(NA)
}

#########################################
###
# Write data into a file at slaves
###
writeLinesSF <- function(data, fileName)
{
	newFile <- file(fileName, "w")
	writeLines(data,newFile)
	close(newFile)
}

###
# Get HeaderDetails from slaves
###
ReadHeaderSF <- function(object)
{
	if( length(object) != 0 )
		return( .Call("ReadHeader", as.character(object[[1]]), PACKAGE = "affyio") )
	else
		return(NA)
}

####################################
###
# getObjectType gets type from object
###
getObjectType <- function(object)
{
	if ( class(object) == "AffyBatch" )
		object.type <- "AffyBatch"
	else if ( class(object) == "character" && is.vector(object) )
		object.type <- "CELfileVec"
	else if ( class(object) == "list" && is.vector(object) )
		object.type <- "partCELfileList"
	else
		stop("Object has an unknown type")
	
	return(object.type)
}

###
# checkPartSize checks object for length
###
checkPartSize <- function(object, number.parts)
{
	if ( class(object) == "list" && is.vector(object) ){
		object.length <- length(object)
		if (object.length < number.parts)
			number.parts = object.length
		else if ( object.length  > number.parts)
			stop("Partitioned CEL file list is longer as dimension of Cluster")
		else {
			number.parts <- length(object)
			object.length <- length(unlist(object))
		}	
		return( list(number.parts=number.parts, object.length=object.length) )
	} else
		return( list(number.parts=number.parts, object.length=length(object)) )
}
	
	