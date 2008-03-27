# file: affyParaIntern.R
# 
# generall support SLAVE - functions for affyPara
#
# History
# 27.03.2008 : ... old stuff removed ...
# 11.12.2007 : Version 0.4 - file.name improved
# 28.12.2007 : Version 0.5 - readHeader added
# 03.01.2008 : Version 0.6 - getFUNAffyBatchSF added
# 27.03.2008 : Version 0.7 - checkObjectType and checkPartSize added
# 
# Copyright (C) 2008 : Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de>
###############################################################################

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
	else 
		AffyBatch <- ReadAffy(filenames=object)

	#temporary save AffyBatch
	assign("AffyBatch", value=AffyBatch, envir= .GlobalEnv)
	return(TRUE)
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
# Create Directory at slaves
###
createDirSF <- function(to,
                        showWarnings = TRUE, recursive = TRUE)
{
	dir.create(to, showWarnings = showWarnings, recursive = recursive)
}

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
# List files from directory at slaves
###
listFilesSF <- function(path, full.names=TRUE)
{
	list.files(path=path, full.names=full.names)
}

###
# Get HeaderDetails from slaves
###
ReadHeader <- function(object)
{
	if( length(object) != 0 )
		return( .Call("ReadHeader", as.character(object[[1]]), PACKAGE = "affyio") )
	else
		return(NA)
}

#################################
###
# Get file name from filecharacter with path
###
file.name <- function(file, fsep = .Platform$file.sep)
{
	i <- nchar(file)
	while(i>0){
		string <- substr(file,i,i)
		if (string==fsep){
			output <- substr(file,i+1,nchar(file))
			return(output)	
		}
		i=i-1
	}
	return(file)
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
	
	