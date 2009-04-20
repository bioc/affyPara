# file: MAplotParaFunctions.R
# 
# Additional functions for the parallelization of the maplot function
#
# History
# 10.10.2008 : Version 0.1 - added to package
# 14.11.2008 : Version 0.2 - fix small bug and plot parameters improvements
# 18.11.2008 : Version 0.3 - add function to convert the results of the 'bad' quality samples as matrix 
# 28.11.2008 : Version 0.4 - getMatrixBQLevels function improved and fix small bugs to plot
# 01.12.2008 : Version 0.5 - function getBoxplot improvd for the parameter plot=FALSE
# 04.12.2008 : Version 0.6 - fix small bugs for the function getMatrixBQLevels  when number of Samples are smaller than number of slaves
# 04.12.2008 : Version 0.7 - set.seed for verbose > 3 added in getValuesChips 
# 11.12.2008 : Version 0.8 - improve MAplot graphics - number of plots per page is limited to 8.
# 23.03.2009 : Version 0.9 - Option verbose set to getOption("verbose") and added . to names of internatl functions
#
# Copyright (C) 2008 : Esmeralda Vicedo <e.vicedo@gmx.net>, Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de> 
###############################################################################
 

########################################################
# function to calculate the mean Array at the slaves
#
#########################################################
MAplotParaSFgMean <- function(type, 
                               log)
{
  if (exists("AffyBatch", envir = .GlobalEnv)) {
	         require(affy)
           affyBatch <- get("AffyBatch", envir = .GlobalEnv)
           #there are 3 possibilties to take the intensities, as : pm, mm and both
           if (type == "both"){
              pms <- unlist(indexProbes(affyBatch, "both"))
            } else if (type == "pm"){
              pms <- unlist(pmindex(affyBatch))
            } else if (type == "mm"){
              pms <- unlist(mmindex(affyBatch))
            }
           if(log){
              x <- log2(intensity(affyBatch)[pms, ])
            } else {
              x <- intensity(affyBatch)[pms, ]
            }
          
            
            assign("x", x, envir= .GlobalEnv)
            #it possible that only a chip at the slave has been sand
           if(!is.null(dim(x))){
               
           #calculate reference Chip as mean from every Chip at the slave
            aB.meanchip <- rowMeans(x)
          }else{  #whne onbly a chip thne we take it as reference Chip and not mean is necessary        
            aB.meanchip <- x
           } 

           return(aB.meanchip)
	 }else
       return(NA)
}



########################################################
# function to get the "bad" Quality Chips at every slave
#
#########################################################
getValuesChips <- function(meanchip, 
                               type, 
                               log, 
                               cuttof, 
                               subset, 
                               span, 
                               family.loess, 
                               verbose=getOption("verbose"))
{                            
   M<- NULL
   affyBatch <- NULL
   intensChip <- NULL
   s_absSumM <- NULL
   if (exists("AffyBatch", envir = .GlobalEnv)) {
     require(affy)
     affyBatch <- get("AffyBatch", envir = .GlobalEnv)
     if (exists("x", envir = .GlobalEnv)){
       intensChip <-  get("x", envir = .GlobalEnv)
    
     }else{
       if (type == "both"){
              pms <- unlist(indexProbes(affyBatch, "both"))
            } else if (type == "pm"){
              pms <- unlist(pmindex(affyBatch))
            } else if (type == "mm"){
              pms <- unlist(mmindex(affyBatch))
            }
       if(log){
              intensChip <- log2(intensity(affyBatch)[pms, ])
            } else {
              intensChip <- (intensity(affyBatch)[pms, ])
            }    
        
     }
     #Calculate the M for one or more Chips
     if(!is.null(dim(intensChip))){
       M <- sweep(intensChip,1,meanchip,FUN='-')
       A <- 1/2*sweep(intensChip,1,meanchip,FUN='+')
      
       # to calculate the loess.Smoother    
       J <- dim(intensChip)[2]          
         
     } else {
       M <- (intensChip - meanchip)
       A <- 1/2*(intensChip + meanchip)
       J <- 1            
     }
     
     #The absolute sum of M for every Prove on a Chip(array) will be calculated as described in the outlier
     # detection from quality Metrics : S_i= Sum_k abs(M_ik) M_ik : log(ratio) of probe k on array i
    
     qualityValues<-  vector("list",6)
     names(qualityValues) <-c("sampleName","S", "sigma", "aproxAM_Y","Var_loess", "Var_sigma")
     subsets <- NULL
     if(J > 1){
      for(j in 1:J){
		  #For Debugging
		  if(verbose>2) set.seed(1234)
       if(is.null(subset)) subsets=sample(1:length(M[,j]),min(c(10000, length(M[,j]))))       
       else subsets = subset
       #Median and sigma for every Sample 
       sigma <- IQR(M[,j])                                                                                       
       mean <- median(M[,j])
       # Summe absolute values ob M for every sample
       s_absSumM <- sum(abs(M[,j]))
                                                                                           
       #check for every calculate loess.smother                                                              
       aux <- loess(M[,j][subsets]~A[,j][subsets],degree=1,span=span,family=family.loess)$fitted                     
       o <- order(A[,j][subsets])                                                                                
       A.sub <- A[,j][subsets][o]
       M.sub <- aux[o]                                                                                        
                                                                                  
       o <-which(!duplicated(A.sub))                                                  
       # to drawn the loess.smooth line                                             
       aproxAM<- approx(A.sub[o],M.sub[o])                                          
                                                                                    
       # aproxAM delivered two list of values:                                      
       # $x -> exes x (A values)                                                    
       # $y <- exes y (M values)                                                    
       # to decide if the loess.smooth line oscilated, we take the $y value         
                                                                                    
       osc.loess <- is.unsorted(aproxAM$y)                                          
       if(verbose) print(paste("osc.loess : ", osc.loess, "\n"))                    
       # the sigam value of the statics MAplot will be considered as a              
       Var.sigma <- sigma > cuttof 
                                                        
       if(verbose) print(paste(" Var.sigma : ",  Var.sigma, "\n" ))                 
       #to write all calculates values in a list and return it as result
       if(length(qualityValues$sampleName) > 0)                         
         qualityValues$sampleName <- c(qualityValues$sampleName, sampleNames(affyBatch[,j]))
       else  
          qualityValues$sampleName <- sampleNames(affyBatch[,j])
       if(length(qualityValues$S) > 0)                         
         qualityValues$S<- c(qualityValues$S, s_absSumM)
       else  
          qualityValues$S <- s_absSumM 
       if(length(qualityValues$sigma) > 0)                         
         qualityValues$sigma <- c(qualityValues$sigma, sigma)
       else  
          qualityValues$sigma <- sigma
       if(length(qualityValues$aproxAM_Y )> 1)                         
         qualityValues$aproxAM_Y <- cbind(qualityValues$aproxAM_Y, aproxAM$y)
       else  
          qualityValues$aproxAM_Y <- aproxAM$y                                                                                                                                                                                                                                         
       if(length(qualityValues$Var_loess) > 0)                         
         qualityValues$Var_loess <- c(qualityValues$Var_loess, osc.loess)
       else  
         qualityValues$Var_loess <- osc.loess                                                                                                                                                                                                                                         
                                                                                                                                                                                         
       if(length(qualityValues$Var_sigma ) > 0)                         
         qualityValues$Var_sigma <- c(qualityValues$Var_sigma, Var.sigma)
       else  
         qualityValues$Var_sigma <- Var.sigma                                                                         
           
     }
     return(qualityValues) 
    }else{    
       
         if(is.null(subset)) subsets=sample(1:length(M),min(c(10000, length(M))))        
         else subsets = subset 
         sigma <- IQR(M)                                                                                            
         mean <- median(M)
         # Summe absolute values ob M 
         s_absSumM <- sum(abs(M))                                                                                          
         #check for every calculate loess.smother                                                                   
         aux <- loess(M[subsets]~A[subsets],degree=1,span=span,family=family.loess)$fitted                          
         o <- order(A[subsets])                                                                                     
         A.sub <- A[subsets][o]                                                                                     
         M.sub <- aux[o]                                                                                            
        
        
         o <-which(!duplicated(A.sub))                                             
         # to drawn the loess.smooth line                                          
         aproxAM<- approx(A.sub[o],M.sub[o])                                       
                                                                                   
         # aproxAM delivered two list of values:                                   
         # $x -> exes x (A values)                                                 
         # $y <- exes y (M values)                                                 
         # to decide if the loess.smooth line oscilated, we take the $y value      
                                                                                   
         osc.loess <- is.unsorted(aproxAM$y)                                       
         if(verbose) print(paste("osc.loess : ", osc.loess, "\n"))                 
         # the sigam value of the statics MAplot will be considered as a           
         Var.sigma <- sigma > cuttof                                               
         if(verbose) print(paste(" Var.sigma : ",  Var.sigma, "\n" ))              
         #to write TRUE if the array is a "bad" quality array 
                           
          qualityValues$sampleName <- sampleNames(affyBatch)                 
          qualityValues$S <- s_absSumM
          qualityValues$sigma <- sigma                      
          qualityValues$aproxAM_Y <- aproxAM$y                      
          qualityValues$Var_loess <- osc.loess                      
          qualityValues$Var_sigma <- Var.sigma    
      
      
          return(qualityValues)
        
     }  
        
   }else
     return(NA)
}

########################################################
# function to get the "bad" Quality Chips from the S values 
#       
#########################################################
getBoxplot <- function(sValues, 
                      verbose=getOption("verbose"), plot=FALSE)
{
 S.badQC<-NULL
 #calculate the statistical values
 if(plot) op<-par(mfrow=c(2,2)) 
 Sboxpl <- boxplot(sValues, plot=plot, main="S boxplot to detect outliers")
 
 S.stats <- Sboxpl$stats
 #take out only the Lower 'whisker' and Upper 'whisker'. These values will 
 # be considered as limit to classified the Samples as "outliers"("bad" quality Arrays)
 S.statsWL <- Sboxpl$stats[1,] 
 S.statsWU <- Sboxpl$stats[5,]
 if(plot){
 abline(h= S.statsWL, lty=1, lwd=1, col="red")
 abline(h= S.statsWU, lty=1, lwd=1, col="red")
 par(op)
 }
 # index of the "outliers" Samples are saved together
 S.badQC <- which(sValues < S.statsWL)  
 S.badQC<- c(S.badQC,  which(sValues > S.statsWU))
 # return only the name of the samples which are considered "outliers" from its S values
 return(S.badQC)
 
 
}

########################################################
# function to get the "bad" Quality Chips from:
# a) Loess Smooth (TRUE/FALSE) or 
# var.sigma (TRUE/FALSE)
#########################################################
getBadQCLoessSigma <- function(lSValues, 
                          verbose=getOption("verbose"))
{
    lS.badQC<- which(lSValues == 1) 
    return(lS.badQC)
}

 ########################################################
# function to give out the index/Name of the arrays in affybatch, which are classified as "bad" quality. The samples will be classified
# in three levels: 
#
# 1 - badQC.sLoessSigma : samples which are classified as "bad" in checkBadQC.S, checkBadQC.loess and checkBadQC.sigma
# 2 - badQC.sloess badQC.sSigma / badQC.loessSigma :  samples which are classified as "bad" only in two of three checkBAdCQ group (S- Loess, S-sigma, Loess-sigma)
# 3 - badQC.loess/ badQC.S / badQC.sigma(samples which are classified as "bad" only in checkBadQC.loess 
#return a list with only the index of the samples
#########################################################


getLevelsBQ<- function(checkBadQC.s, checkBadQC.loess, checkBadQC.sigma, verbose=getOption("verbose") ){

 #Samples of the second level but the first level isn�t considered yet:
  badQC.sLoess <- sort(intersect(checkBadQC.s, checkBadQC.loess))
  badQC.sSigma<- sort(intersect(checkBadQC.s, checkBadQC.sigma))
  badQC.loessSigma<- sort(intersect(checkBadQC.loess, checkBadQC.sigma))
  
  #Samples of the first level:  intersect (sLoes(s , loess) , sigma) the three methods
  badQC.sLoessSigma <- sort(intersect(badQC.sLoess ,checkBadQC.sigma))
 
  #Samples of the second level: the samples of the second level whitout first level: 
  #badQC.Sec_sLoess <- setdiff(badQC.sLoess, badQC.sLoessSigma) 
 # badQC.Sec_sSigma <- setdiff(badQC.sSigma, badQC.sLoessSigma) 
  #badQC.Sec_loessSigma <- setdiff(badQC.loessSigma, badQC.sLoessSigma) 
  
  #badQC.secondL <- list(badQC.Sec_sLoess, badQC.Sec_sSigma, badQC.Sec_loessSigma)
  badQC.secondL <-list(badQC.sLoess, badQC.sSigma, badQC.loessSigma)
  names(badQC.secondL) <- c("s-loess", "s-sigma", "loess-sigma")
  
  #Samples of the third level:
  #badQC.unLoessSigma <- union(checkBadQC.loess, checkBadQC.sigma)
  #badQC.s_LoessSigma <- sort(setdiff(checkBadQC.s, badQC.unLoessSigma))
  
 # badQC.unSSigma <-  union(checkBadQC.s, checkBadQC.sigma)
  #badQC.l_SSigma <- sort(setdiff(checkBadQC.loess, badQC.unSSigma))
  
  #badQC.unSLoess <-  union(checkBadQC.s, checkBadQC.loess)
  #badQC.sigma_SLoess <- sort(setdiff(checkBadQC.sigma, badQC.unSLoess))
  
 # badQC.thirdL <- list(badQC.s_LoessSigma,  badQC.l_SSigma, badQC.sigma_SLoess)
  
  badQC.thirdL<- list(checkBadQC.s, checkBadQC.loess, checkBadQC.sigma)
  names(badQC.thirdL) <- c("s", "loess", "sigma")
  
  #the three levels are grouped to be ploted
  badQC.MAplots <- list(badQC.sLoessSigma, badQC.secondL, badQC.thirdL)
  names( badQC.MAplots) <- c("firstLevel", "secondLevel", "thirdLevel")

  return(badQC.MAplots)

}

########################################################
# drawMAplot function to draw MAplot from only the "bad quality Arrays.
# object: Affybatch object 
# index: parameter which contains the index of the samples to be ploted
# meanchip : calculated reference chip as mean 
# type: which kind of intensities should be considered :  "PM", "MM", "both"
# log : logical - if true the intensities should be calculated as log2
# subset: used to calculate the loess.Smoother line 
#########################################################
drawMAplot <- function(object, 
                       indexSamp, 
                       meanchip, 
                       type, 
                       log, 
                       subset,
                       spans, 
                       ref.title, 
                       pchs,
                       show.statistics,
                       family.loess,
                       verbose=getOption("verbose"),
                       ...)
{
      
    lindex<- length(indexSamp)
    namesInd<- names(indexSamp)
    lobject<- length(object) 
    if (verbose) cat(paste("Samples to be ploted: ", lindex, "\n"))   
    # it is needed to calculate the M and A values only for the Arrays to be ploted     
    if (type == "both"){                                 
      pms <- unlist(indexProbes(object, "both"))        
    } else if (type == "pm"){                           
      pms <- unlist(pmindex(object))                    
    } else if (type == "mm"){                           
      pms <- unlist(mmindex(object))                    
    }                                                   
   
    
     
    
     if (verbose) cat(paste("Samples to be ploted : " , lindex,  "\n"))                          
        
   
     #to prepare the graphic window for the number of necessary MAplots
     frame()
     old.par <- par(no.readonly = TRUE)
     on.exit(par(old.par))
     par(mfrow=c(4,2),mgp=c(0,.2,0),mar=c(1,1,1,1),oma=c(1,1.4,2,1))
     #to send individually the 'bad' quality index to the ma.plot function
     
   for(i in 1:lindex){
              
       #calculate the intensity for every array. it is necessary to calculate the M and A values    
      if(log){
          x <- log2(intensity(object[,indexSamp[i]])[pms, ])
      } else {
          x <- intensity(object[,indexSamp[i]])[pms, ]
      }
           
      M <- x - meanchip                                                        
      A <- (x + meanchip)/2 
      
      if(verbose) cat("Sample ", i , " will be ploted \n")
      if(is.null(subset)) subsets=sample(1:length(M),min(c(10000, length(M))))
      else subsets = subset
     
      main.plot<-paste("Sample : ", indexSamp[i], "-", ref.title)
      
      
      
      #call the generic function ma.plot for only the array i 
      check<- ma.plot(A, M, subset=subsets, show.statistics =show.statistics, span=spans, family.loess=family.loess, pch=pchs, main=main.plot, ...)
      if(verbose) print(check)
          
    }
     invisible()
   
}


########################################################                     
# function convert the results of the 'bad' quality samples as matrix                              
# to facility the further possible calculations and analysis                                        
#                                        
#########################################################                    
 

getMatrixBQLevels <- function(calQCSampleName, badQC.MAplots){

     lSN <- length(calQCSampleName)
     levelsNam <- names(badQC.MAplots)
     lbQC <- length(levelsNam)
     default<-rep(0,lSN) 
     #built matrix
     sampleLevelAll <- matrix( c(calQCSampleName, rep(default,lbQC)),  nrow=lSN , ncol=2)
     #second and third level Quality will be splited for a better analysis
    
     
     #remplace the default values with the correct value: 0 is FALSE and 1 TRUE
     tempNamCol <- NULL
     temp<- 1
     for(i in 1: lbQC){
         
         NamIndexL<- sub("^[a-zA-Z]*[.]","",names(unlist(badQC.MAplots[levelsNam[i]])))
         
         NamIndex <- unique(sub("[1-9]*[0-9]*$","", NamIndexL)) 
         li <- length(NamIndex)
         #for subclasses of a class
         if(li > 1){         
              for( j in 1: li){
                 
                 temp<- temp +1
                 if(dim(sampleLevelAll)[2]< temp) sampleLevelAll<- cbind(sampleLevelAll, rep(0,lSN))                
                 sampleLevelAll[badQC.MAplots[[i]][[NamIndex[[j]]]], temp] <- 1
                 tempNamCol <- c(tempNamCol, NamIndex[j])
                 
              }   
         }else if(li ==1){
                temp<- temp +1
                if(dim(sampleLevelAll)[2]< temp) sampleLevelAll<- cbind(sampleLevelAll, rep(0,lSN))
                tempNamCol <- c(tempNamCol, NamIndex)
                if(data.class(badQC.MAplots[[i]]) == "list")  sampleLevelAll[badQC.MAplots[[i]][[NamIndex]],temp] <- 1 
                else sampleLevelAll[badQC.MAplots[[i]],temp] <- 1 
         }
         
                    
      }
     colnames(sampleLevelAll) <- c("sampleNames", tempNamCol)
     #second and third level Quality will be splited for a better analysis
     return(sampleLevelAll) 
     
 }
 
 
