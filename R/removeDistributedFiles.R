# removeDistributedFiles.R
# 
# Function for removing data from local slave directory
#
# History
# 03.08.2007 : Version 0.1
# 03.09.2007 : Version 0.2 - use of "snow", copying by different protocols
# 20.11.2007 : Version 0.3 - code cleaning and global functions
# 27.11.2007 : Version 0.4 - bugFixes for variable cluster, seperation for os not necessary
#                            codeCleaning
# 05.12.2007 : Version 0.5 - bugFixes files.remove -> unlink
# 06.12.2007 : Version 0.6 - code dokumentation
# 12.12.2007 : Version 0.7 - output fixed
#
# Copyright (C) 2008 : Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de>
###############################################################################

removeDistributedFiles <- function(cluster, 
		path, verbose=FALSE)
{

	#Check for snow
	require(snow)	
	
	#Check cluster
	checkCluster(cluster)
    
    #Remove Files at Master
	if (verbose) cat("Remove files on Master ")
		unlink(path, recursive=TRUE)
	if (verbose) cat("... DONE\n")

	#Remove Files at Slaves
	if (verbose) cat("Remove files on Slaves ")
		clusterCall(cluster, removeFilesSF, path)
	if (verbose) cat("... DONE\n")
		
	#Check if cleanup was successfull
	if (verbose) cat("Check: ")
		masterCheck <- file.access(path, mode = 0)
		slaveCheck <- unlist(clusterCall(cluster, accessFilesSF, path))
		anzSlaves <- length(cluster)
		if (sum(slaveCheck) == (-1 * anzSlaves) && masterCheck == -1) {
      output <- TRUE
	    if (verbose) cat("Files successfully removed!\n");
		} else {
		  output <- FALSE
			warning("Files not successfully removed!\n")
    }

  #Return result
  return( output )
}

###
# Slavefunction
# remove Files at slaves
###
removeFilesSF <- function(path, recursive=TRUE) {
	unlink(path, recursive=recursive)
}

###
# Slavefunction
# access Files at slaves
###
accessFilesSF <- function(path, mode = 0){
	file.access(path, mode = mode)
}