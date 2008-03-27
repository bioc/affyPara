# file: distributeFiles.R
# 
# Function for (hierarchical) distribution of data from the master node to every slave.
#
# History
# 01.03.2008 : ... old stuff removed ...
# 27.11.2007 : Version 0.6 - bugFixes for variable cluster
#                            type R partly fixed, foldername removed
# 05.12.2007 : Version 0.7 - output fixed, full.names added
# 06.12.2007 : Version 0.8 - input changed to files and more Infooutput added
# 12.12.2007 : Version 0.9 - hierarchicallyDist set to FALSE
# 19.12.2007 : Version 0.10 - bugfixes for protocol = R
# 07.03.2008 : Version 0.11 - bugfixes in verbose and for more nodes at one workstation
# 17.03.2008 : Version 0.12 - extending protocol for scp, adding master distribution
#
# Copyright (C) 2008 : Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de>
###############################################################################

distributeFiles <- function(cluster, 
		files, to="/usr1/tmp/CELfiles",
		protocol=c("R","RCP","SCP"), hierarchicallyDist=FALSE,
		master=TRUE,
		full.names=FALSE, verbose=FALSE)
{

	#Check for snow
	require(snow)
	
	protocol <- match.arg(protocol)

    ###################
	#Create directories
	##################
	if (verbose>0) cat("Create Directories ")
		t1 <- proc.time()
		error <- system(paste("mkdir ",to,sep=""))
		error <- clusterCall(cluster, createDirSF, to)
		t2 <- proc.time()
	if (verbose>0) cat(round(t2[3]-t1[3],3),"sec DONE\n")
	
	####################
	#Partition of files
	###################
	if (verbose) cat("Partition of files ")
		t0 <- proc.time();
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
				filename <- file.name(files[j])
				newFile <- file(paste(to,filename,sep="/"), "w")
				writeLines(data,newFile)
				close(newFile)
			}
		} else if (protocol == "RCP" || protocol == "SCP") {
			for (j in files){
				if (verbose>0) cat(".")
				if (protocol == "RCP") {
					error <- system(paste("rcp ", j," ",to,"/",sep=""))
				} else if (protocol == "SCP") {
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
						filename <- file.name(filesPart[[i]][j])
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
						error <- system(paste("rcp ", j," ",host,sep=""))
					} else if (protocol == "SCP") {
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
				#TODO machbar über c[1]
				warning("hierarchically dist wirh R not yet implemented")
			}
			t2 <- proc.time()
		if (verbose>0) cat(round(t2[3]-t1[3],3),"sec DONE\n")
	}
	
	#######
	#Output
	#######
    filesSlaves <- clusterCall(cluster, listFilesSF, to, full.names=full.names)
	return( list(to=to, CELfiles=filesSlaves) )
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