# removeDistributedFiles.R
# 
# Function for removing data from local slave directory
#
# History
# 27.10.2008 : ... old stuff removed ...
# 12.10.2008 : Version 0.10 - code cleaning and bugfix in file removing (only files, not directory)
# 23.10.2008 : Version 0.11 - cleanup improved
# 24.10.2008 : Version 0.12 - Improvements for multiprocessor machines
# 27.10.2008 : Version 0.13 - Output improved to work at multiprocessor machines and clusters
# 18.12.2008 : Version 0.14 - cluster object gets default parameter: .affyParaInternalEnv$cl
# 23.03.2009 : Version 0.15 - Option verbose set to getOption("verbose") and added . to names of internatl functions
#
# Copyright (C) 2009 : Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de>
###############################################################################

removeDistributedFiles <- function(path=tempdir(),
		cluster, verbose=getOption("verbose"))
{
	#Check for snow
	require(snow)
	
	#Get cluster object form default environment
	if(missing(cluster))
		cluster <- .affyParaInternalEnv$cl
	
	#Check cluster
	checkCluster(cluster)
	
	#check for multiprocessor machine -> no removement!
	nodenames <-unlist(clusterEvalQ(cluster, Sys.info()["nodename"]))
	if( all(nodenames == nodenames[1]) ){
		if(verbose) cat("No data removement: you use a multiprocessor machine!\n")
	} else{	
	    
	    #Remove Files at Master, but nor directroy
		if (verbose) cat("Remove files on Master ")
			check<-unlink(paste(path,"*",sep="/"), recursive=TRUE)
			if (check == 1)
				check<-unlink(paste(path,"*",sep="/"), recursive=TRUE)
		if (verbose) cat("... DONE\n")
	
		#Remove files and directories at slaves
		if (verbose) cat("Remove files on Slaves ")
			clusterCall(cluster, unlink, path, recursive=TRUE)
		if (verbose) cat("... DONE\n")
			
		#Check if cleanup was successfull
		if (verbose) cat("Check: ")
			masterCheck <- length(list.files(path))
		 	slaveCheck <- sum( unlist(clusterCall(cluster, function(x) length(list.files(x)), path) ) )
			if ( slaveCheck == 0 && masterCheck == 0) {
				output <- TRUE
				if (verbose) cat("Files successfully removed!\n");
			} else {
				output <- FALSE
				if (verbose) cat("Files not successfully removed!\n");
				warning("Files not successfully removed!\n")
			}
	
		  #Return result
		  return( output )
	  }
}