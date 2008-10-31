# file: MAplotPara.R
# 
# Parallelization of the maplot function
#
# History
# 31.10.2008 : Version 0.1 - added to package
#
# Sending AffyBatch form master to slave an back is very time consuming. Sending a list
# of CEL files from master to slave, creating the AffyBatch and do BG-Correction is faster.
# Using the right combination "size of AffyBatch on slaves" - "number of slaves" the parallelized
# version is more than ten times faster as the serial version. 
#
#
# Copyright (C) 2008 : Esmeralda Vicedo <e.vicedo@gmx.net>, Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de> 
###############################################################################
MAplotPara <- function(cluster,
                       object, 
                       log=TRUE, 
                       type=c("both","pm","mm"),
                       ref=NULL,
                       which=NULL, 
                       ref.title="vs pseudo-median reference chip",
                       subset=sample(1:length(M),min(c(10000, length(M)))), 
                       #ref.fn=c("median","mean"),  parameter will be not used because we calculate the reference as mean
                       show.statistics=TRUE,
                       span=2/3,
                       family.loess="gaussian",
                       cuttof =0.5,
                       verbose =FALSE
                       )
{    
  #Check for affy
	require(affyPara)
	#get total number of samples   
  lobject<- length(object)
 
   #Check cluster and generate number.parts
	checkCluster(cluster)
	number.parts <- length(cluster)
	
	#Check object type
	object.type <- getObjectType(object) 
	
	#Check size of partitions
	parts <- checkPartSize(object, number.parts)
	number.parts <- parts$number.parts
	object.length <- parts$object.length
	
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
	#return(samples.names)  			
	t1 <- proc.time();
	if (verbose) cat(paste(round(t1[3]-t0[3],3),"sec DONE\n"))				
	
	#Info-Output for Distribution
	if (verbose){ cat("Object Distribution: "); cat(paste(lapply(object.list,length))); cat("\n") }
	
	#################################
	#Initialize AffyBatches at slaves
	##################################
	if (verbose) cat("Initialize AffyBatches at slaves \n")
	t0 <- proc.time();
	#remove all variables from all slaves
  check <- clusterCall(cluster, function(){rm(list=ls(envir = .GlobalEnv), envir = .GlobalEnv)})
  if (verbose) print(check)	  
  #send the Cel files liste to the slaves
	check <- clusterApply(cluster, object.list, initAffyBatchSF, object.type) 
	if (verbose) print(check)	  
	t1 <- proc.time();
	if (verbose) cat(paste(round(t1[3]-t0[3],3),"sec DONE\n"))
		
	############################
	# Do constant MAplot
	############################
	
  # check the local settigs to avoid problems with sort
	check <- clusterCall(cluster, Sys.setlocale, "LC_COLLATE", Sys.getlocale("LC_COLLATE"))
	
  if (verbose) cat("Box-Parameter parallel calculation \n")
	t2 <- proc.time()


	#check the parameter which
  if (is.character(which)){
                which.indices <- match(which,sampleNames(object))
                if (all(is.na(which.indices))){
                  stop("No known sampleNames in which")
                }
                
                if (any(is.na(which.indices))){
                  warning(paste("Omitting the following from which:",which[is.na(which.indices)], "because they can not be found."))
                }
                which <- which.indices[!is.na(which.indices)]
              }
              
              
              if (is.null(which)){
                which <- 1:dim(exprs(object))[2]
              }
  
  type <- match.arg(type)
  # if a reference chip is given take it as median chip otherwise calculated the median-chip 
  # at the slaves. 
    if(!is.null(ref)){
      if (type == "both"){
              pms <- unlist(indexProbes(object, "both"))
            } else if (type == "pm"){
              pms <- unlist(pmindex(object))
            } else if (type == "mm"){
              pms <- unlist(mmindex(object))
            }
            
      if(log){
              x <- log2(intensity(object)[pms, ])
            } else {
              x <- intensity(object)[pms, ]
            }     
      #meanchip from reference chip gives as parameter  
      meanchip <- x[,ref]                   
    }else{
      #function "MAplotParaSFgMean" calculate the mean default Array at the slaves       
       meanClusterChip  <- clusterCall(cluster, MAplotParaSFgMean, type,log) 
   
	
       if (verbose) print( meanClusterChip)	 
	
       if (verbose) cat("Mean refence Chip merge\n")
  # from every slaves is returned a reference calculated array (intensities). It will be merged as only one. 
       TmpMeanChip <- NULL
       lmChip <- length( meanClusterChip)
      
  #merge all obtained reference Chips in one column
       for(i in 1:lmChip){
        if(length(TmpMeanChip) < 1){
         TmpMeanChip <-  meanClusterChip[[i]] 
        }else {
         TmpMeanChip <- cbind(TmpMeanChip, meanClusterChip[[i]]) 
        }
      }
   } 
  #calculate the mean for all column. The result is only a value for every intensity and so we obtain the mean reference array
  meanchip <- rowMeans(TmpMeanChip, na.rm=TRUE)
  #meanchip values are too necessary at the slaves to calculate the M
  badQualChips<- clusterCall(cluster,getBadQualityChips, meanchip, type, log, span, family.loess, cuttof)
  #to give out the name of the "bad" Quality arrays 
  checkBadQC.S<- getBoxplot(badQualChips)
  checkBadQC.loess<- getBadQCLoess(badQualChips)
  maplot.BadQC <-  
 
  return (badQualChips) 
}