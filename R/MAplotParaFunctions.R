# file: MAplotParaFunctions.R
# 
# Additional functions for the parallelization of the maplot function
#
# History
# 10.10.2008 : Version 0.1 - added to package
#
# Copyright (C) 2008 : Esmeralda Vicedo <e.vicedo@gmx.net>, Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de> 
###############################################################################

########################################################
# function to calculate the mean Array at the slaves
#
#########################################################
MAplotParaSFgMean <- function(type, log)
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

getBadQualityChips <- function(meanchip, type, log, span, family.loess, cuttof ){
  M<- NULL
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
     } else {
       M <- (intensChip - meanchip)
       A <- 1/2*(intensChip + meanchip)
     }
     
     #The absolute sum of M for every Prove on a Chip(array) will be calculated as described in the outlier
     # detection from quality Metrics : S_i= Sum_k abs(M_ik) M_ik : log(ratio) of probe k on array i
     s_absSumM <- colSums(abs(M))
     qualityProb <- NULL
     arraysNames <- NULL
     #calculate the loess.Smoother
     J <- dim(intensChip)[2]
     for(j in 1:J){
      sigma <- IQR(M[,j])
      mean <- median(M[,j])
      print(paste(j, "sigma: ", sigma, "\n"))
      subset=sample(1:length(M[,j]),min(c(10000, length(M[,j]))))
      aux <- loess(M[,j][subset]~A[,j][subset],degree=1,span=span,family=family.loess)$fitted
      #check for every calculate loess.smother 
      o <- order(A[,j][subset])
      A.sub <- A[,j][subset][o]
      M.sub <- aux[o]
      o <-which(!duplicated(A.sub))
      # to drawn the loess.smooth line
      aproxAM<- approx(A.sub[o],M.sub[o])
      
      osc.loess <- is.unsorted(aproxAM$x)
      print(paste(j, "osc.loess : ", osc.loess, "\n"))
      Var.sigma <- sigma > cuttof
      print(paste(j, " Var.sigma : ",  Var.sigma, "\n" ))
      #to write TRUE if the array is a "bad" quality array
      if( osc.loess | Var.sigma){
        qualityProb <- c(qualityProb,TRUE)
      }else{
        qualityProb <- c(qualityProb,FALSE)
      }
       osc.loess <- NULL
       Var.sigma <- NULL    
     }
     
     return(list(s_absSumM, qualityProb))
  
   }else
     return(NA)
}


getBoxplot<- function(badQualChips){
 lbad <- length(badQualChips)
 S <- NULL

 S.badQC <- NULL     
 for(i in 1:lbad){
  S <- c(S, unlist(badQualChips[[i]][1]))
 
 } 

 Sboxpl <- boxplot(S, plot=FALSE)
 S.stats <- Sboxpl$stats

 S.statsHL <- Sboxpl$stats[2,] 
 S.statsHU <- Sboxpl$stats[4,]


 S.badQC <- which(S < S.statsHL)  
 S.badQC<- c(S.badQC,  which(S > S.statsHU))

 return(names(S.badQC))
 
}


getBadQCLoess <- function(badQualChips){
  loessSm <- NULL         
  namesChips <- NULL
  badQC <- NULL
  badQCIndex<- NULL
  lbad <- length(badQualChips)
           
  for(i in 1:lbad){   
   loessSm <- c(loessSm, unlist(badQualChips[[i]][2])) 
   namesChips<- c(namesChips, names(unlist(badQualChips[[i]][1])))    
  }
    badQC<- which(loessSm == 1)      
    badQCIndex<- namesChips[badQC]
    return(badQCIndex)
}