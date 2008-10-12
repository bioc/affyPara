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
# 02.07.2008 : Version 0.8 - default value for path added
# 10.10.2008 : Version 0.9 - path set to tempdir()
# 12.10.2008 : Version 0.10 - code cleaning and bugfix in file removing (only files, not directory)
#
# Copyright (C) 2008 : Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de>
###############################################################################

removeDistributedFiles <- function(cluster, 
		path=tempdir(), verbose=FALSE)
{

	#Check for snow
	require(snow)	
	
	#Check cluster
	checkCluster(cluster)
    
    #Remove Files at Master
	if (verbose) cat("Remove files on Master ")
		check<-unlink(paste(path,"*",sep="/"), recursive=TRUE)
		if (check == 1)
			check<-unlink(paste(path,"*",sep="/"), recursive=TRUE)
	if (verbose) cat("... DONE\n")

	#Remove Files at Slaves
	if (verbose) cat("Remove files on Slaves ")
		clusterCall(cluster, removeFilesSF, path)
	if (verbose) cat("... DONE\n")
		
	#Check if cleanup was successfull
	if (verbose) cat("Check: ")
		masterCheck <- file.access(path, mode = 0)
		slaveCheck <- unlist(clusterCall(cluster, file.access, path, mode=0))
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
	unlink(paste(path,"*",sep="/"), recursive=recursive)
}