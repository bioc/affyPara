# split.R
# 
# Functions to split objects (affyBatch, file vector and matrix) into parts for distributed computing
#
# History
# 22.08.2007 : Version 0.1
# 06.12.2007 : Version 0.2 - code cleaning and dokumentation
#
# We do not use clusterSplit, because there is no way to set the amount of parts.
#
# Copyright (C) 2008 : Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de>
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
	for (i in 1:number.part)
		abatch.list[[i]] <- abatch[,samples.parts[[i]]]
	
	return( abatch.list )
}

# for file vector
splitFileVector <- function(FileVec, number.part)
{
	require(snow)
	
	#Number of Files
	File.count <- length(FileVec)
	
	if (File.count < number.part)
		number.part = File.count

 	#Partitioning 		
	File.parts <- splitIndices(File.count, number.part)
	File.list <- list()
	for (i in 1:number.part)
		File.list[[i]] <- FileVec[File.parts[[i]]]
		
	return( File.list )
}

# for matrix
splitMatrix <- function(matrix, number.part)
{
	require(snow)
	
	#Partitioning by Columns = Samples 	
	return( splitCols(matrix, number.part) )
}
