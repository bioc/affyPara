# file: boxplotPara.R
# 
# Parallelization of the affybatch.boxplot function
#
# History
# 10.10.2008 : Version 0.1 - added to package
# 17.10.2008 : Version 0.2 - parameter plot - gives the possibility to choose output between only statistical computation 
#                            or both: statical and graphics
# 18.10.2008 : Version 0.3 - code and output cleaning
# 18.10.2008 : Version 0.4 - Documentation File move to .Rd file
# 12.11.2008 : Version 0.5 - rename von plotDraw to plot
#
# Sending AffyBatch form master to slave an back is very time consuming. Sending a list
# of CEL files from master to slave, creating the AffyBatch and do BG-Correction is faster.
# Using the right combination "size of AffyBatch on slaves" - "number of slaves" the parallelized
# version is more than ten times faster as the serial version. 
#
# Copyright (C) 2008 : Esmeralda Vicedo <e.vicedo@gmx.net>, Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de> 
###############################################################################

boxplotPara <- function(cluster,
		object,
 	  nSample=if(length(object)> 200) nSample<-200 else nSample <- length(object),
		iqrMethod=TRUE,
    percent=0.05,
		typDef="mean",
    	plot=TRUE,	 
		plotAllBoxes=TRUE,
		verbose=FALSE) 
{
   
	#Check for affy
	require(affyPara)
	#get total number of samples   
  lobject<- length(object)
  ##############################
  #Check Functions parameter 
  ##############################
  
  # check typDef :c("median" , "mean")	
  if(!boxplotParaCheckTypDef(typDef)) stop("Please check the Parameter TypDef. It should be : median, mean")
  # check iqrMethod : logical, TRUE as default,                   
  if(!iqrMethod){ 
    if(!boxplotParaCheckPercent(percent)) stop("Please check the Parameter percent. It should be a number > 0 and <1")
    else{
    if(boxplotParaCheckPercentnSample(percent,nSample)== 0)stop("Please check the value of percent.It gives 0 Samples to be considered as critical") 
    if(boxplotParaCheckPercentnSample(percent,nSample)== 2){
      message("There arent enought Samples to calculate the ", percent, " The Quality Problem Samples will be calculated with the IQR method")
      iqrMethod <- TRUE
      }
    }  
  }        
  if(!boxplotParaChecknSamples(nSample, lobject)) stop("Please check the Parameter nSample. It should be a number < 200 and < total samples: ", lobject)
 
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
	t1 <- proc.time();
	if (verbose) cat(paste(round(t1[3]-t0[3],3),"sec DONE\n"))				
	
	#Info-Output for Distribution
	if (verbose){ cat("Object Distribution: "); cat(paste(lapply(object.list,length))); cat("\n") }
	
	#################################
	#Initialize AffyBatches at slaves
	##################################
	if (verbose) cat("Initialize AffyBatches at slaves ")
	t0 <- proc.time();
 	#send the Cel files liste to the slaves# and remove data from slaves
	check <- clusterApply(cluster, object.list, initAffyBatchSF, object.type, rm.list="ALL")  
	if (verbose > 1) print(check)	  
	t1 <- proc.time();
	if (verbose) cat(paste(round(t1[3]-t0[3],3),"sec DONE\n"))
		
	############################
	# Do constant boxplot
	############################
	
	# check the local settigs to avoid problems with sort
	check <- clusterCall(cluster, Sys.setlocale, "LC_COLLATE", Sys.getlocale("LC_COLLATE"))
	
  	if (verbose) cat("Box-Parameter parallel calculation ")
	t0 <- proc.time()
	#send the function boxplotParaSFgParBox to clusters to calculate the boxplot(object, plot=false)  
	box.list <- clusterCall(cluster, boxplotParaSFgParBox)
	t1 <- proc.time()
  	if (verbose>1) print(box.list)	 
	if (verbose) cat(paste(round(t1[3]-t0[3],3),"sec DONE\n"))
  	
	if (verbose) cat("Box-Parameter merge ")
    t0 <- proc.time()
  	# function to summarize the boxplot-stats values from slaves 
  	boxpl.st <- boxplotParaComStatsM(box.list, verbose)
    t1 <- proc.time()
	if (verbose) cat(paste(round(t1[3]-t0[3],3),"sec DONE\n"))
	
  	#to save the calculated boxplot.stats values in a file when verbose = 2)
  	if(verbose > 1) save(boxpl.st, file="boxplSt.Rdata") 
 
  #Calculate the default Sample from boxplot.stats$stats object : boxpl.st$stats
  #first we check which type of Sample should be calculated : "median", "mean" or "both", as default we take the both   
  if(verbose>1) cat(paste("Length of typ Default: ",  length(typDef), " Name: ", typDef, "\n"))                                                                                                                             
 
  mstats <- 0                                                                                                            
  defaultS <- 0                                                                                                          
                  
  defaultS <- boxplotParagMdMnDef(boxpl.st, typDef)                                                                      
                                                                                                                            
  if(length(defaultS) < 3){                                                                                               
    mstats <- matrix(defaultS[["median"]]$stats, ncol=1, nrow=5)                                                          
                                                                                                                           
    #add mean values to the matrix stats and vector names                                                                 
    mstats <- cbind(mstats, defaultS[["mean"]]$stats)                                                                     
    colnames(mstats) <- c("median", "mean")                                                                               
                                                                                                                            
  }else{                                                                                                                  
    if(defaultS$names =="median")                                                                                         
      mstats <- defaultS$stats                                                                                           
                                                                                                                            
    if(defaultS$names =="mean")                                                                                           
      mstats <- defaultS$stats                                                                                           
                                                                                                                            
  }                                                                                                                       
  if(verbose > 1) save(defaultS, file="defaultSample.Rdata") 
    
    ###########################################################################################
    #                                                                                         #
    #  CALCULATE QUALITY PROBLEM SAMPLES from - IQRs BOXPLOT (median and IQR) for all Samples #
    #  as parameter iqrMethod = TRUE (as default)at the function boxplotPara()                #
    ###########################################################################################

 #parameter iqrMethod = TRUE (default)  
 if(iqrMethod){
    if(verbose) cat("Calculate the Bad Quality Samples after the IQRs Samples ")
    t9 <- proc.time()
    #matrix to save the media, HL and HU from the calculated default sample
    medIQRAllM <- matrix(c(boxpl.st$stats[2,], boxpl.st$stats[4,], boxpl.st$stats[3,]), ncol=3, nrow=length(boxpl.st$n)) 
    colnames(medIQRAllM) <- c("HL", "HU", "median")
    #calculate the median and HL , HU values to be considered as limit between "normal" samples and bad quality samples
    limitSamplesBxp <- boxplotParagMedIQR(medIQRAllM, plot, verbose)
    if(verbose > 1) save( limitSamplesBxp, file="limitSamplesBxp.Rdata")
    #take only the median and IQR values of all samples
    diffMedIQRALLM <- matrix(c(abs(medIQRAllM[,1] - medIQRAllM[,2]), medIQRAllM[,3]), ncol=2, nrow=length(boxpl.st$n))
    colnames(diffMedIQRALLM) <- c("Diff_IQR", "median")
    #calculate the critical samples from the limits and all samples 
    criticSamplesBxp <- boxplotParagCrSamp(diffMedIQRALLM, limitSamplesBxp, iqrMethod, verbose)
    #check the different obtained values and classified them in 3 classes(mdIQR, md, IQR) 
    qualityProblemBxp <- boxplotParacheckCritSamp(criticSamplesBxp, iqrMethod, verbose)
    t10 <- proc.time()
    if (verbose) cat(paste(round(t10[3] - t9[3],3)," sec DONE\n"))
  	t11 <- proc.time()
   	#if plot parametere== TRUE , drawn the boxplot 
	  if (plot) { 	
    	if(verbose) cat("Drawn the boxplots with Bad Quality Samples ")     
    	boxplotParaDrawn(boxpl.st, defaultS, qualityProblemBxp, limitSamplesBxp, nSample, plotAllBoxes, verbose)
		}
	  t12 <- proc.time()
	  if (verbose) cat(paste(round(t12[3] - t11[3],3),"sec DONE\n"))
	  if(verbose>1) cat("Total Time necessary to calculate and draw boxplotPara :", round((t12[3] +t11[3] +t10[3] +t9[3])- t1[3],3), "sec\n")  
    # to merge the statistical results for all Samples whit the problematic samples together
    boxpl.st$QualityPS.IQR <- qualityProblemBxp     
    return(boxpl.st)
  } 
   
   
     
    ########################################################################################
    #                                                                                      #
    #  CALCULATE QUALITY PROBLEM SAMPLES as DIFFERENCEN of medians and IQRs                #
    #  (THE ALTERNATIVE METHODE, only if the Data quality aren�t good or bad)              #
    #                                                                                      #
    #  If the iqrMethod parameter equal FALSE then the Quality problem samples will be     # 
    #  calculate from the difference between samples and DefaultSample values              #
    ########################################################################################

  # parameter iqrMethod = FALSE 
 if(!iqrMethod){ 
   if(verbose) cat("Calculate the Bad Quality Samples between the standart sample and the rest ")
   t5 <- proc.time()
   #calculate the differences between the Default calculated Sample(defaultS) and the rest of the samples   
   differencen<- boxplotParaCalAllDiff(boxpl.st$stats, mstats, verbose)
   if(verbose) cat(differencen)
   # to save the difference values
   if(verbose > 1) save(differencen, file="differencen.Rdata")
   #calculate  the levels for a concretes percent which is gives from boxplotPara function as parameter
   limits<- boxplotParagLimits(differencen, percent, plot, verbose)
   if(verbose) cat(limits) 
   #calculate the critical samples from the limits       
   criticSamples <- boxplotParagCrSamp(differencen, limits, iqrMethod, verbose)
   #check the different obtained values and classified them in 3 classes(mdIQR, md, IQR) 
   qualityProblem <- boxplotParacheckCritSamp(criticSamples, iqrMethod, verbose)
    
   if(verbose > 1) save(qualityProblem, file="qualProblem.Rdata")
   t6 <- proc.time()
   if(verbose) cat(paste(round(t6[3] - t5[3],3),"sec DONE\n"))    
   if(verbose) cat("Drawn the boxplots with Bad Quality Samples")
   t7 <- proc.time()
   #drawn the boxplot
   if (plot) 
   	boxplotParaDrawn(boxpl.st,defaultS, qualityProblem, limits, nSample, plotAllBoxes, verbose)
   t8 <- proc.time()
   if(verbose) cat(paste(round(t8[3] - t7[3],3), "sec DONE\n"))  
	 if(verbose > 1)cat("Total Time necessary to calculate and drawn boxplotPara : ", round((t8[3] +t7[3] + t6[3] +t5[3]) - t1[3],3), " sec\n") 
   #to merge the statistical results for all Samples whit the problematic samples together
   boxpl.st$QualityPS.Diff <- qualityProblem
   
   return(boxpl.st)
 }


   
  }                            