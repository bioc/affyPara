# file: distributeFiles.R
# 
# Function for (hierarchical) distribution of data from the master node to every slave.
#
# History
# 27.10.2008 : ... old stuff removed ...
# 02.07.2008 : Version 0.15 - delExistTo added
# 11.09.2008 : Version 0.16 - file.name removed -> basename
# 10.10.2008 : Version 0.17 - to set to tempdir()
# 23.10.2008 : Version 0.18 - Improvements in directory creation
# 24.10.2008 : Version 0.19 - Improvements for multiprocessor machines
# 27.10.2008 : Version 0.20 - Output improved to work at multiprocessor machines and clusters
# 18.12.2008 : Version 0.21 - cluster object gets default parameter: .affyParaInternalEnv$cl
#
# Copyright (C) 2008 : Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de>
###############################################################################

distributeFiles <- function( files, to=tempdir(),
		protocol=c("R","RCP","SCP"), hierarchicallyDist=FALSE,
		master=TRUE, delExistTo=FALSE,
		full.names=TRUE, 
		cluster, verbose=FALSE)
{
	#Check for snow
	require(snow)
	
	#Get cluster object form default environment
	if(missing(cluster))
		cluster <- .affyParaInternalEnv$cl
	
	#Check cluster
	checkCluster(cluster)
	
	#check for multiprocessor machine -> no distribution!
	nodenames <-unlist(clusterEvalQ(cluster, Sys.info()["nodename"]))
	if( all(nodenames == nodenames[1]) ){
		if(verbose) cat("No data distribution: you use a multiprocessor machine!\n")
		if(full.names==TRUE)
			return( list(to=dirname(files)[1], CELfiles=paste(dirname(files)[1],basename(files),sep="/") ))
		else
			return( list(to=dirname(files)[1], CELfiles=basename(files)) )
	} else{	
		
		protocol <- match.arg(protocol)
		
		###################
		#Delete directories
		###################
		if (delExistTo == TRUE){
			if (verbose>0) cat("Remove Directories ")
			t1 <- proc.time()
			check <- removeDistributedFiles(cluster, path=to, verbose=FALSE)
			t2 <- proc.time()
			if (verbose>0) cat(round(t2[3]-t1[3],3),"sec DONE\n")
		}
		
	    ###################
		#Create directories
		###################
		if (verbose>0) cat("Create Directories ")
			t1 <- proc.time()
			#at master: only if not tempdir
			if(to != tempdir()){
				error <- system(paste("mkdir ",to,sep=""))
				if ( error != 0) warning("Directory ",to," already exists at master.")
			}
			#at slaves
			error <- clusterCall(cluster, dir.create, to, showWarnings = TRUE, recursive = TRUE)
			check <- lapply(error, function(x){ if(x!=TRUE) warning("Directory ",to," already exists at slaves.")} )
			
			t2 <- proc.time()
		if (verbose>0) cat(round(t2[3]-t1[3],3),"sec DONE\n")
		
		###################
		#Partition of files
		###################
		if (verbose) cat("Partition of files ")
			t0 <- proc.time();
			if ( length(cluster) == 1 ){
				filesPart<-list()
				filesPart[[1]] <- files
			} else
				filesPart <- clusterSplit(cluster,files)
			t1 <- proc.time();
		if (verbose) cat(round(t1[3]-t0[3],3),"sec DONE\n")
		#Info-Output for Distribution
		if (verbose){ cat("\tFile Distribution: "); cat(paste(lapply(filesPart,length))); cat("\n") }
		
		#####################
		#Move Files to master
		#####################
		#Copy all files to master node
		if (master == TRUE)
		{
			if (verbose>0) cat("Move Files to Master ")
			t1 <- proc.time()
			MasterNode <- Sys.info()["nodename"]
			if (protocol == "R") {
				for (j in 1:length(files)){
					if (verbose>0) cat(".")
					data <- readLines(files[j], n=-1)
					filename <- basename(files[j])
					newFile <- file(paste(to,filename,sep="/"), "w")
					writeLines(data,newFile)
					close(newFile)
				}
			} else if (protocol == "RCP" || protocol == "SCP") {
				for (j in files){
					if (verbose>0) cat(".")
					if (protocol == "RCP") {
						if(verbose>1) cat(paste("rcp ", j," ",to,"/",sep=""),"\n")
						error <- system(paste("rcp ", j," ",to,"/",sep=""))
					} else if (protocol == "SCP") {
						if(verbose>1) cat(paste("scp ", j," ",to,"/",sep=""),"\n")
						error <- system(paste("scp ", j," ",to,"/",sep=""))
					}
				}
			}
			t2 <- proc.time()
			if (verbose>0) cat(round(t2[3]-t1[3],3),"sec DONE\n")
		}
		
		#####################
		#Move Parts to slaves
		#####################
		#Copy every Part to one slave
		if (verbose>0) cat("Move Parts to Slaves ")
			t1 <- proc.time()
			#list of nodes
			nodes <- clusterEvalQ(cluster,Sys.info()["nodename"])
			anzNodes <- length(nodes)
			if (verbose>1) print(filesPart)
			
			#Copy parts of data to slaves
			if (protocol == "R") {	
				for( i in 1:length(filesPart)) {
					for (j in 1:length(filesPart[[i]])) {
						if( length(filesPart[[i]][j])!= 0 && file.exists(filesPart[[i]][j]) ){
							if (verbose>0) cat(".")
							data <- readLines(filesPart[[i]][j], n=-1)
							filename <- basename(filesPart[[i]][j])
							error <- clusterCall(cluster[i], writeLinesSF, data, paste(to,filename,sep="/"))
						}
					}
				}
			} else if (protocol == "RCP" || protocol == "SCP") {
				for (i in 1:anzNodes){
					host <- paste(nodes[[i]],":",to,"/",sep="")
					if (verbose>0) cat(".")
					for (j in filesPart[[i]]){
						if (protocol == "RCP") {
							if(verbose>1) cat(paste("rcp ", j," ",host,sep=""),"\n")
							error <- system(paste("rcp ", j," ",host,sep=""))
						} else if (protocol == "SCP") {
							if(verbose>1) cat(paste("scp ", j," ",host,sep=""),"\n")
							error <- system(paste("scp ", j," ",host,sep=""))
						}
					}
				}
			}	
			t2 <- proc.time()
		if (verbose>0) cat(round(t2[3]-t1[3],3),"sec DONE\n")
	
	    ##########################################
		#Hierarchical distribution between slaves
		##########################################
		# copy parts from slaves to all other slaves
		if (hierarchicallyDist==TRUE){
			if (verbose>0) cat("Hierarchical distribution ")
				t1 <- proc.time()
				if (protocol == "RCP" || protocol == "SCP") {	
					out <- clusterCall(cluster, hirarchicalDistSF, to, nodes, protocol=protocol)
					if (verbose>1) print(out)		
				} else if (protocol == "R") {
					warning("TODO hierarchically dist with R not yet implemented")
				}
				t2 <- proc.time()
			if (verbose>0) cat(round(t2[3]-t1[3],3),"sec DONE\n")
		}
		
		#######
		#Output
		#######
		if (hierarchicallyDist==TRUE){
			filesPart <- clusterCall(cluster, list.files,to,full.names=TRUE)
		}		
		if(full.names==TRUE)
			filesSlaves <-lapply(filesPart, function(x,to){
						if(length(x)!=0) paste(to,basename(x),sep="/")
					},to)
		else
			filesSlaves <-lapply(filesPart, basename)
		return( list(to=to, CELfiles=filesSlaves) )
	}
}

###
# Slavefunction
# Hirarchical distribution
###
hirarchicalDistSF <- function(to,nodes, protocol)
{
	nodename <- Sys.info()["nodename"]
	#listing of existing files on slave
	files <- list.files(path=to, full.names=TRUE)
	# nodes randomication to reduce load at nodes
	nodes <- nodes[sample(length(nodes))]
	out<-list()
	#Distribution
	for (i in files){
		out[[i]]<-list()						
		for (j in nodes){
			if (nodename != j){
				if (protocol == "RCP") {
					out[[i]][j] <- paste("rcp ",i," ",j,":",to,"/",sep="")
					system(paste("rcp ",i," ",j,":",to,"/",sep=""))
				} else if (protocol == "SCP") {
					out[[i]][j] <- paste("scp ",i," ",j,":",to,"/",sep="")
					system(paste("scp ",i," ",j,":",to,"/",sep=""))
				}
			}
		}
	}
	return(out)
}