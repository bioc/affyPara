# split.R
# 
# Functions to split objects (affyBatch, file vector and matrix) into parts
# for distributed computing
#
# History
# 22.08.2007 : Version 0.1
# 06.12.2007 : Version 0.2 - code cleaning and dokumentation
# 13.05.2008 : Version 0.3 - splitAffyBatch improved for small clusters
# 16.05.2008 : Version 0.4 - splitFileVector improved for small clusters
#
# We do not use clusterSplit, because there is no way to set the amount of parts.
#
# Copyright (C) 2008 - 2010 : Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de>
###############################################################################

# for AffyBatch
splitAffyBatch <- function(abatch, number.part)
{
	require(snow)

	#Number of Samples
	samples.count <- length(abatch)
	
	if (samples.count < number.part)
		number.part = samples.count
	
	#Partitioning 
	samples.parts <- splitIndices(samples.count, number.part)
	abatch.list <- list()
	if (number.part == 1)
		abatch.list[[1]]<- abatch
	else for (i in 1:number.part)
			abatch.list[[i]] <- abatch[,samples.parts[[i]]]
	
	return( abatch.list )
}

# for file vector
splitFileVector <- function(fileVec, number.part)
{
	require(snow)
	
	#Number of Files
	file.count <- length(fileVec)
	
	if (file.count < number.part)
		number.part = file.count

 	#Partitioning 		
	file.parts <- splitIndices(file.count, number.part)
	file.list <- list()
	if (number.part == 1)
		file.list[[1]]<- fileVec
	else for (i in 1:number.part)
		file.list[[i]] <- fileVec[file.parts[[i]]]
		
	return( file.list )
}

# for matrix
splitMatrix <- function(matrix, number.part)
{
	require(snow)
	
	#Partitioning by Columns = Samples 	
	return( splitCols(matrix, number.part) )
}
