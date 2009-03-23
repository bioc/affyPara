# file: MAPara.R
# 
# Parallelization of the affybatch.MAplot function
#
# History
# 31.10.2008 : Version 0.1 - added to package
# 18.11.2008 : Version 0.2 - drawMAplot function improved (any parameter are optimized)- fix small bug whit the list badQC.
# 19.11.2008 : Version 0.3 - results of the MAplotPara function are returned as list with 3 elements: 1-matrix with the values of MAplot, 2-matrix contained the loess values for the axis y,  3- matrix wiht the result of the MAplot analysis  
# 28.11.2008 : Version 0.4 - add $qualityMA to the result list to a better control during analysis 
# 01.12.2008 : Version 0.5 - fix small bugs: parameter ref get a numerical value,  the additional parameters for the routine(...), level 3 to plot only unique samples.
# 04.12.2008 : Version 0.6 - fix small bugs when the number of Samples are smaller than number of slaves.
# 05.12.2008 : Version 0.7 - fix small bug when parameter plot=FALSE.
# 11.12.2008 : Version 0.8 - fix small bug for large Affybatch object.
# 18.12.2008 : Version 0.9 - cluster object gets default parameter: .affyParaInternalEnv$cl
# 23.03.2009 : Version 0.10 - Option verbose set to getOption("verbose") and added . to names of internatl functions
#
# Sending AffyBatch form master to slave an back is very time consuming. Sending a list
# of CEL files from master to slave, creating the AffyBatch and do BG-Correction is faster.
# Using the right combination "size of AffyBatch on slaves" - "number of slaves" the parallelized
# version is more than ten times faster as the serial version. 
#
# Copyright (C) 2008 : Esmeralda Vicedo <e.vicedo@gmx.net>, Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de> 
###############################################################################
MAplotPara <- function(object, 
                       log=TRUE, 
                       type=c("both","pm","mm"),
                       ref=NULL,
                       which=NULL, 
                       ref.title="vs mean ref.Array",
                       subset=NULL,
                       span=1/4,
                       show.statistics=TRUE,
                       family.loess ="gaussian", 
                       pch=".",
                       plot =TRUE,
                       cutoff =0.5,# add parameter to generic function ma.plot
                       level=1,
					   cluster, verbose=getOption("verbose"),
                       ...                       
                       )
{
	#Check for affy
	require(affy)
	require(snow)
	
	#Get cluster object form default environment
	if(missing(cluster))
		cluster <- .affyParaInternalEnv$cl
	
  	#Check cluster and generate number.parts
	checkCluster(cluster)
	number.parts <- length(cluster)
	
	#Check object type
	object.type <- .getObjectType(object) 
	
	#Check size of partitions
	parts <- .checkPartSize(object, number.parts)
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
  #send the Cel files liste to the slaves
	check <- clusterApply(cluster, object.list, .initAffyBatchSF, object.type)
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
  
   
	#to hold the parameter type 
  type <- match.arg(type)
  # if a reference chip is given take it as median chip otherwise calculated the mean-chip 
  # at the slaves. 
    if(!is.null(ref) & is.numeric(ref)){
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
      #function "MAplotParaSFgMean" calculate the mean reference Array at the slaves       
       meanClusterChip  <- clusterCall(cluster, MAplotParaSFgMean, type,log) 
   
       if (verbose) cat("Mean refence Array merge\n")
      #from every slaves is returned a reference calculated array (as intensities). It will be merged as only one Array. 
       tmpMeanChip <- NULL
       lmChip <- length( meanClusterChip)
      
      #merge all obtained reference Chips in one column
       if(lmChip > 1){
       for(i in 1:lmChip){
        if(length(tmpMeanChip) < 1){
         tmpMeanChip <-  meanClusterChip[[i]] 
        }else {
         tmpMeanChip <- cbind(tmpMeanChip, meanClusterChip[[i]]) 
        }
       }
      }else{
        tmpMeanChip <- meanClusterChip
      } 
     #calculate the mean for all column. The result is only a value for every intensity and so we obtain the mean reference array
     meanchip <- rowMeans(tmpMeanChip, na.rm=TRUE)
   } 
  
  if (verbose >1 ) save(meanchip, file="meanchipMAplot.Rdata")
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
  if (verbose >1) save(calValuesChips, file="calValueSamplesMA.Rdata")
  if (verbose) print(calValuesChips)
  lCV <- length(calValuesChips)
  # to take only the necessary values to get the "bad" quality Arrays from the total of calValuesChips
  calQCS <- NULL
  calQCL<- NULL
  calQCsigma <- NULL
  calQCSampleName <- NULL
  calQCaproxAM_Y <- NULL
  calQCVarsigma <- NULL
  
  if(is.null(names(calValuesChips[[1]]))){
      calQCSampleName<-  calValuesChips$sampleName
      calQCS <-  calValuesChips$S 
      calQCsigma<- calValuesChips$sigma
      calQCaproxAM_Y <- calValuesChips$aproxAM_Y
      calQCL <- calValuesChips$Var_loess
      calQCVarsigma <-  calValuesChips$Var_sigma
   
  }else{
     for(i in 1: lCV){
     
      if(length(calValuesChips[[i]])>2){ 
      calQCSampleName<- c(calQCSampleName, calValuesChips[[i]]$sampleName)
      calQCS <- c(calQCS, calValuesChips[[i]]$S) 
      calQCsigma<- c(calQCsigma, calValuesChips[[i]]$sigma)
      calQCaproxAM_Y <- cbind(calQCaproxAM_Y, calValuesChips[[i]]$aproxAM_Y)
      calQCL<- c(calQCL, calValuesChips[[i]]$Var_loess)
      calQCVarsigma<- c(calQCVarsigma, calValuesChips[[i]]$Var_sigma)
      }
      
     }
   }
   calQC<- matrix(c(calQCSampleName,calQCS,calQCL, calQCsigma, calQCVarsigma), nrow=length(calQCSampleName), ncol=5)
   colnames(calQC) <- c("sampleNames","S","osc_Loess","sigma", "var_sigma")
  #Samples classified as "bad" after the S value (outliers) 
  checkBadQC.s<- getBoxplot(calQCS,plot)
  #Samples classified as "bad" after the loess.smooth line 
  checkBadQC.loess<- getBadQCLoessSigma(calQCL)
  #Samples classified as "bad" after the sigma value
  checkBadQC.sigma <- getBadQCLoessSigma(calQCVarsigma)
  
  #to give out the index/Name of the arrays in affybatch, which are classified as "bad" quality. The samples will be classified
  # in three levels: 
  # 1 - badQC.sLoessSigma : samples which are classified as "bad" in checkBadQC.S, checkBadQC.loess and checkBadQC.sigma
  # 2 - badQC.sloess badQC.sSigma / badQC.loessSigma :  samples which are classified as "bad" only in two of three checkBAdCQ group (S- Loess, S-sigma, Loess-sigma)
  # 3 - badQC.loess/ badQC.S / badQC.sigma(samples which are classified as "bad" only in checkBadQC.loess 
  
  badQC.MAplots <- getLevelsBQ(checkBadQC.s, checkBadQC.loess, checkBadQC.sigma )
  
  
  if(verbose >1 ) save(badQC.MAplots, file="qualityProbleMA.Rdata")
 	t5 <- proc.time();
  if (verbose) cat(paste(round(t4[3]-t5[3],3),"sec DONE\n"))
  if(plot){                                                                                  
     if (verbose) cat("MA-plots with the calculated bad quality samples will be done \n")  
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
  
      if (verbose) cat(paste("bad quality Array:", badQC, "\n"))  
     
      if(length(badQC) > 0){
        badQCu<- unique(badQC)      
        check <- drawMAplot(object, badQCu, meanchip,type, log, subset, span, pch=pch, ref.title,show.statistics, family.loess, verbose, ...)
        if(verbose) print(check)
      }else{
         
        print(paste("MAplots aren't built. Reason: 'bad ' quality samples at the level: ",level, " aren't found \n" ))
      
      }
    
   t7 <- proc.time();
   if (verbose) cat(paste(round(t7[3]-t6[3],3),"sec DONE\n"))
   
  }else{ if (verbose) cat(paste("MAplots aren't built. Reason: 'plot' parameter is given as: ", plot, "\n" ))}
  
  levelBQMatrix <- getMatrixBQLevels(calQCSampleName, badQC.MAplots) 
  #return list of the values and results obteined from the MAplotPAra function
  resultQC <- list(calQC, calQCaproxAM_Y, badQC.MAplots, levelBQMatrix)
  names(resultQC)<- c("values_MAP","loess_y", "quality_MAP", "results_MAP") 
  
  t8 <- proc.time();
  if (verbose) cat(paste("bad-quality Samples are calculated in ", round(t8[3]-t0[3],3)," sec"))
  
  return(resultQC)
}