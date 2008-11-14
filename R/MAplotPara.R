# file: MAPara.R
# 
# Parallelization of the affybatch.MAplot function
#
# History
# 31.10.2008 : Version 0.1 - added to package
# 14.11.2008 : Version 0.2 - code cleaning and further improvements
#							 RD file added
#
# Sending AffyBatch form master to slave an back is very time consuming. Sending a list
# of CEL files from master to slave, creating the AffyBatch and do BG-Correction is faster.
# Using the right combination "size of AffyBatch on slaves" - "number of slaves" the parallelized
# version is more than ten times faster as the serial version. 
#
# Copyright (C) 2008 : Esmeralda Vicedo <e.vicedo@gmx.net>, Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de> 
###############################################################################

MAplotPara <- function(cluster,
                       object, 
                       log=TRUE, 
                       type=c("both","pm","mm"),
                       ref=NULL,
                       which=NULL, 
                       ref.title="vs mean ref.Array",
                       subset=NULL,
                       span=1/4,
                       family.loess="gaussian",
                       plot.method=c("normal","smoothScatter","add"),
                       pch=".",
                       plot =TRUE,
                       cutoff =0.5,# add parameter to generic fucntion ma.plot
                       level=1,
                       verbose =FALSE,
                       ...                       
                       )
{
	 #Check for affy
	require(affyPara)
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
  	#send the Cel files liste to the slaves and remove all variables from all slaves
	check <- clusterApply(cluster, object.list, initAffyBatchSF, object.type, rm.list="ALL")
	if (verbose) print (check)	  
	t1 <- proc.time();
	if (verbose) cat(paste(round(t1[3]-t0[3],3),"sec DONE\n"))
		
	############################
	# Do constant MAplot
	############################
	# check the local settigs to avoid problems with sort
	check <- clusterCall(cluster, Sys.setlocale, "LC_COLLATE", Sys.getlocale("LC_COLLATE"))
	
  	if (verbose) cat("Mean reference Chip parallel calculation \n")
	t2 <- proc.time()
    # to hold the optional parameter "..." of the function
  	fn.call <- list(...)  
	#check the parameter which ??????
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
    
  #to hold the parameter type 
  type <- match.arg(type)
  # if a reference chip is given take it as median chip otherwise calculated the mean-chip 
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
      #meanchip is taken from reference Array, as parameter in the function 
      meanchip <- x[,ref]                   
    }else{
      #function "MAplotParaSFgMean" calculate the mean default Array at the slaves       
       meanClusterChip  <- clusterCall(cluster, MAplotParaSFgMean, type,log) 
   
       if (verbose) cat("Mean refence Array merge\n")
      #from every slaves is returned a reference calculated array (as intensities). It will be merged as only one Array. 
       tmpMeanChip <- NULL
       lmChip <- length( meanClusterChip)
      
  #merge all obtained reference Chips in one column
       for(i in 1:lmChip){
        if(length(tmpMeanChip) < 1){
         tmpMeanChip <-  meanClusterChip[[i]] 
        }else {
         tmpMeanChip <- cbind(tmpMeanChip, meanClusterChip[[i]]) 
        }
      }
   } 
  #calculate the mean for all column. The result is only a value for every intensity and so we obtain the mean reference array
  meanchip <- rowMeans(tmpMeanChip, na.rm=TRUE)
 	t3 <- proc.time();
	if (verbose) cat(paste(round(t3[3]-t2[3],3),"sec DONE\n"))
  #meanchip values are necessary at the slaves to calculate the "bad" quality samples
 	t4 <- proc.time();
  if (verbose) cat("Bad Quality samples will be computed")
  calValuesChips<- clusterCall(cluster, getValuesChips,meanchip, type, log, cutoff, subset, span, family.loess, verbose)
  #                                                      
  #to detect the name of the "bad" Quality arrays. Two methods are used to classify the Arrays(samples):
  # 1 - as calculated with quality Metrics methods: S_i= Sum_k abs(M_ik) M_ik : log(ratio) of probe k on array i. These samples are
  #     classified as checkBadQC.S
  # 2 -  througth the loess Smoother line from the MA-plot. The samples with a oscillated loess line or 
  #      IQR(M)(written as sigma in the MA-statisc) or both are classified as checkBadQC.loess 
  # 
  
  if (verbose >1) save(calValuesChips, file="calculateValuesChips.Rda")
  if (verbose) print(calValuesChips)
  lCV <- length(calValuesChips)
  # to take only the necessary values to get the "bad" quality Arrays from the total of calValuesChips
  calQCS <- NULL
  calQCL<- NULL
  calQCsigma <- NULL
  calQCSampleName <- NULL
   for(i in 1: lCV){
      calQCS <- c(calQCS, calValuesChips[[i]]$S) 
      calQCL<- c(calQCL, calValuesChips[[i]]$Var_loess)
      calQCsigma<- c(calQCsigma, calValuesChips[[i]]$Var_sigma)
      calQCSampleName<- c(calQCSampleName, calValuesChips[[i]]$sampleName)
   }
   calQC<- list(calQCS,calQCL, calQCsigma, calQCSampleName)
  #Samples classified as "bad" after the S value (outliers) 
  checkBadQC.s<- getBoxplot(calQCS, verbose)
  #Samples classified as "bad" after the loess.smooth line 
  checkBadQC.loess<- getBadQCLoessSigma(calQCL, verbose)
  #Samples classified as "bad" after the sigma value
  checkBadQC.sigma <- getBadQCLoessSigma(calQCsigma, verbose)
  
  #to give out the index/Name of the arrays in affybatch, which are classified as "bad" quality. The samples will be classified
  # in three levels: 
  # 1 - badQC.sLoessSigma : samples which are classified as "bad" in checkBadQC.S, checkBadQC.loess and checkBadQC.sigma
  # 2 - badQC.sloess badQC.sSigma / badQC.loessSigma :  samples which are classified as "bad" only in two of three checkBAdCQ group (S- Loess, S-sigma, Loess-sigma)
  # 3 - badQC.loess/ badQC.S / badQC.sigma(samples which are classified as "bad" only in checkBadQC.loess 
  
  #Samples of the second level:
  badQC.sLoess <- sort(intersect(checkBadQC.s, checkBadQC.loess))
  badQC.sSigma<- sort(intersect(checkBadQC.s, checkBadQC.sigma))
  badQC.loessSigma<- sort(intersect(checkBadQC.loess, checkBadQC.sigma))
  badQC.secondL <- list( badQC.sLoess, badQC.sSigma, badQC.loessSigma)
  names(badQC.secondL) <- c("s-loess", "s-sigma", "loess-sigma")
  #Samples of the first level:
  badQC.sLoessSigma <- sort(intersect(badQC.sLoess ,checkBadQC.sigma))
 
  #Samples of the third level:
  badQC.unLoessSigma <- union(checkBadQC.loess, checkBadQC.sigma)
  badQC.s_LoessSigma <- sort(setdiff(checkBadQC.s, badQC.unLoessSigma))
  
  badQC.unSSigma <-  union(checkBadQC.s, checkBadQC.sigma)
  badQC.l_SSigma <- sort(setdiff(checkBadQC.loess, badQC.unSSigma))
  
  badQC.unSLoess <-  union(checkBadQC.s, checkBadQC.loess)
  badQC.sigma_SLoess <- sort(setdiff(checkBadQC.sigma, badQC.unSLoess))
  
  badQC.thirdL <- list(badQC.s_LoessSigma,  badQC.l_SSigma, badQC.sigma_SLoess)
  names(badQC.thirdL) <- c("s", "loess", "sigma")
  
  #the three levels are grouped to be ploted
  badQC.MAplots <- list(badQC.sLoessSigma, badQC.secondL, badQC.thirdL)
  names( badQC.MAplots) <- c("firstLevel", "secondLevel", "thirdLevel")
  
  
  if(verbose >1 ) save(badQC.MAplots, file="badQC_MAplots.Rda")
 	t5 <- proc.time();
  if (verbose) cat(paste(round(t4[3]-t5[3],3),"sec DONE\n"))
  if(plot){                                                                                  
     if (verbose) cat("MA-plots with the calculated bad quality samples will be done")  
     # to plot only the "bad" quality level samples specified by parameter level 
    	t6 <- proc.time();                                                                        
      badQC <- NULL
      if(level==1){  
        badQC <- unlist( badQC.MAplots$firstLevel)
      }else if (level == 2){
      
          badQC <- c(unlist(badQC.MAplots$firstLevel), unlist(badQC.MAplots$secondLevel))
          
      }else {
          badQC <- unlist(badQC.MAplots)
      }
  
      lbadQC <- length(badQC)
      if( !is.null(badQC)){
        lbadQC <- length(lbadQC)
        n <- ceiling(lbadQC/2)
    
        #to draw plots which are classified at the same level together.
        op<-par(mfrow=c(2,n))
        #the fuction drawMaplot select the data to be ploted
       check <- drawMAplot(object, badQC, meanchip, type, log, subset, span, ref.title, pch, verbose)
       if (verbose)print(check)
      }
     
     
    par(op)
    
 	  t7 <- proc.time();
   if (verbose) cat(paste(round(t7[3]-t6[3],3),"sec DONE\n"))
   
  } 
  #return the Sample Name of the Arrays which are classified and ploted as "bad" quality Arrays.  
  return (list(calQC, badQC.MAplots))
 	t8 <- proc.time();
  if (verbose) cat(paste("bad-quality Samples are calculated in ", round(t8[3]-t0[3],3)," sec")) 
}