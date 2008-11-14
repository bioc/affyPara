# file: MAplotParaFunctions.R
# 
# Additional functions for the parallelization of the maplot function
#
# History
# 10.10.2008 : Version 0.1 - added to package
# 14.11.2008 : Version 0.2 - code cleaning and further improvements
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
                               cutoff, 
                               subset, 
                               span, 
                               family.loess, 
                               verbose)
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
       Var.sigma <- sigma > cutoff 
                                                        
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
         Var.sigma <- sigma > cutoff                                               
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
getBoxplot<- function(sValues, 
                      verbose)
{
 S.badQC<-NULL
 #calculate the statistical values 
 Sboxpl <- boxplot(sValues, plot=FALSE)
 S.stats <- Sboxpl$stats
 #take out only the Lower 'hinge' and Upper 'hinge'. These values will 
 # be considered as limit to classified the Samples as "outliers"("bad" quality Arrays)
 S.statsHL <- Sboxpl$stats[2,] 
 S.statsHU <- Sboxpl$stats[4,]

 # index of the "outliers" Samples are saved together
 S.badQC <- which(sValues < S.statsHL)  
 S.badQC<- c(S.badQC,  which(sValues > S.statsHU))
 # return only the name of the samples which are considered "outliers" from iths S values
 return(S.badQC)
 
 
}

########################################################
# function to get the "bad" Quality Chips from:
# a) Loess Smooth (TRUE/FALSE) or 
# var.sigma (TRUE/FALSE)
#########################################################
getBadQCLoessSigma <- function(lSValues, 
                          verbose)
{
    lS.badQC<- which(lSValues == 1) 
    return(lS.badQC)
}

########################################################
# drawMAplot function to draw MAplot from only the "bad quality Arrays.
# object: Affybatch object 
# index: parameter which contains the index of the samples to be ploted
# meanchip : calculated reference chip as mean 
# type: which kind of intensities should be considered :  "PM", "MM", "both"
# log : logical - if true the intensities should be calculated as log2
# subset: ?? used to calculate the loess.Smoother line 
#########################################################
drawMAplot <- function(object, 
                       indexSamp, 
                       meanchip, 
                       type, 
                       log, 
                       subset,
                       span, 
                       ref.title, 
                       pch, 
                       verbose)
{
      
    lindex<- length(indexSamp)
    namesInd<- names(indexSamp)
    lobject<- length(object) 
    print(paste("length index: ", lindex, "\n"))   
    print(paste("length Obj :" , lobject))
   
      
    if (type == "both"){                                 
      pms <- unlist(indexProbes(object, "both"))        
    } else if (type == "pm"){                           
      pms <- unlist(pmindex(object))                    
    } else if (type == "mm"){                           
      pms <- unlist(mmindex(object))                    
    }                                                   
    lobject<- length(object) 
    
       
    for(i in 1:lindex){
          
      if(log){
          x <- log2(intensity(object[,indexSamp[i]])[pms, ])
      } else {
          x <- intensity(object[,indexSamp[i]])[pms, ]
      }
       
      
      M <- x - meanchip                                                        
      A <- (x + meanchip)/2   
      cat(length(M))
      cat(length(A))  
      if(verbose) cat("Sample ", i , "will be ploted")
      if(is.null(subset)) subsets=sample(1:length(M),min(c(10000, length(M))))
      else subsets = subset
      main.plot<-paste("Sample :", indexSamp[i], ref.title) 
      check<- ma.plot(A, M, subset=subsets, show.statistics =TRUE, span=2/3, plot.method="normal", family.loess ="gaussian",add.loess =TRUE,pch, main=main.plot)
      
      print(check)
          
    }
   
} 