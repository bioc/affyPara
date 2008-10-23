# file: boxplotParaFunctions.R
# 
# Additional functions for the parallelization of the affybatch.boxplot function
#
# History
# 10.10.2008 : Version 0.1 - added to package
# 17.10.2008 : Version 0.2 - parameter plot, !plotAllBoxes (improved with colors), function boxplotParaChecknSamples
# 18.10.2008 : Version 0.3 - code and output cleaning
#
#
# Copyright (C) 2008 : Esmeralda Vicedo <e.vicedo@gmx.net>, Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de> 
# 
###############################################################################

#########################################################################
# This is the slave function to create the boxplot-stats at the slaves
# return the boxplot.stats object calculate from AffyBatch object and
# the function boxplot( plot=false)
#########################################################################
boxplotParaSFgParBox <- function()
{       
  if (exists("AffyBatch", envir = .GlobalEnv)) {
	         require(affy)
           affyBatch <- get("AffyBatch", envir = .GlobalEnv)

           #calculate boxplot statistics for every sample
           affyBatch.stats <- boxplot(affyBatch, plot=FALSE)
           #but despite a error with the $names slot, it will be extra saved from the Affybatch object
           affyBatch.stats$names <- sampleNames(affyBatch)
           return(affyBatch.stats)
	 }else
       return(NA)
}


#########################################################################
# This is the function to  summarize the boxplot-stats values from slaves
# return the boxplot.stats object calculate from AffyBatch object and
# the function boxplot( plot=false)
#########################################################################
boxplotParaComStatsM <- function(matrix.list, 
		verbose=TRUE)
{
  if (verbose) cat("\n\tRebuild Statistic matrix\n")
  #create the summarized affyBatch.stat as null list with the name from boxplot.stats
  sumAffyBatch.stats <- vector("list",6)
  names(sumAffyBatch.stats) <-c("stats", "n","conf", "out", "group", "names")
  index <- length(matrix.list)
  #Dimension of matrix
  if(index > 1){
    #Rebuilding
    for(i in 1:index){
      if(!is.atomic(matrix.list[[i]])){
      if(length( sumAffyBatch.stats$stats)> 1)
        sumAffyBatch.stats$stats <- cbind(sumAffyBatch.stats$stats, matrix.list[[i]]$stats)
  	  else
        sumAffyBatch.stats$stats <- matrix.list[[i]]$stats

      if(length( sumAffyBatch.stats$n)> 0)
        sumAffyBatch.stats$n <- c(sumAffyBatch.stats$n, matrix.list[[i]]$n)
      else
        sumAffyBatch.stats$n <- matrix.list[[i]]$n

      if(length( sumAffyBatch.stats$conf) > 1)
        sumAffyBatch.stats$conf <- cbind(sumAffyBatch.stats$conf, matrix.list[[i]]$conf)
      else
        sumAffyBatch.stats$conf <- matrix.list[[i]]$conf

      if(length( sumAffyBatch.stats$out)> 0)
        sumAffyBatch.stats$out <- c(sumAffyBatch.stats$out, matrix.list[[i]]$out)
      else
        sumAffyBatch.stats$out <- matrix.list[[i]]$out

      if(length( sumAffyBatch.stats$group)> 1)
        sumAffyBatch.stats$group <- c( sumAffyBatch.stats$group, matrix.list[[i]]$group)
      else
        sumAffyBatch.stats$group <- matrix.list[[i]]$group

      if(length( sumAffyBatch.stats$names)> 0)
        sumAffyBatch.stats$names <- c( sumAffyBatch.stats$names, matrix.list[[i]]$names)
      else
        sumAffyBatch.stats$names <- matrix.list[[i]]$names
       }
     }
  }
  return(sumAffyBatch.stats)
 }

#########################################################################
# function to calculate the default sample from boxplot.stats$stats object
#
#########################################################################
boxplotParagMdMnDef <- function(statsObject, 
		typSample, verbose=TRUE)
{

 if( length(typSample) == 1){

  default <- vector("list",6)
  names(default) <-c("stats", "n","conf", "out", "group", "names")

  if( typSample == "median"){
    default$stats <- matrix(apply(statsObject$stats, 1,median), ncol=1, nrow=5)
    colnames(default$stats) <- typSample
    default$n <- sum(statsObject$n)
    default$conf <- matrix(apply(statsObject$conf, 1,median),ncol=1, nrow=2)
    colnames(default$conf) <- typSample

    if( length(statsObject$out) > 0 ){
      default$out <- median(statsObject$out)
    }else{
      default$out <- NULL
    }

    if( length(statsObject$group) > 0 ){
      default$group <- median(statsObject$group)
    }else{
      default$group <- NULL
    }

    default$names <- typSample
   }

   else if(typSample == "mean"){
    default$stats<- matrix(rowMeans(statsObject$stats), ncol=1, nrow=5)
    colnames(default$stats) <- typSample
    default$n <- sum(statsObject$n)
    default$conf <- matrix(rowMeans(statsObject$conf),ncol=1, nrow=2)
    colnames(default$conf) <- typSample
    if( length(statsObject$out) > 0 ){
     default$out <-mean(statsObject$out)
    }else{
     default$out <- NULL
    }

   if( length(statsObject$group) > 0 ){
     default$group <- mean(statsObject$group)
   }else{
     default$group <- NULL
   }
   default$names <- typSample
  }
 }else {
    default <- vector("list",2)
    names(default)<- c("median", "mean")
    default$median <- vector("list",6)
    default$mean <- vector("list",6)
     names(default$median) <- c("stats", "n","conf", "out", "group", "names")
     names(default$mean) <- c("stats", "n","conf", "out", "group", "names")

    default$median$stats<- matrix(apply( statsObject$stats, 1,median), ncol=1, nrow=5)
    default$mean$stats<-matrix(rowMeans(statsObject$stats),ncol=1, nrow=5)

    default$median$n <- sum(statsObject$n)
    default$mean$n <- sum(statsObject$n)

    default$median$conf <- matrix(apply(statsObject$conf, 1,median), ncol=1, nrow=2)
    default$mean$conf <- matrix(rowMeans(statsObject$conf), ncol=1, nrow=2)

    if( length(statsObject$out) > 0 ){
      default$median$out <- median(statsObject$out)   
      default$mean$out <- mean(statsObject$out) 
    }else{
      default$median$out <- 0
      default$mean$out <- 0
    }

    if( length(statsObject$group) > 0 ){
      default$median$group <- median(statsObject$group)  
      default$mean$group <- mean(statsObject$group)
    }else{
      default$median$group <- 0
      default$mean$group <- 0
    }

    default$mean$names <- "mean"
    default$median$names <- "median"

   }

   return (default)
 }



###############################################################################
# Function to calculate the medians and IQRs Difference between default Sample
# and the other Samples. These medians and IQRs are obtained from the boxplot$stats
# parameter: objecte boxplot.stats.$stats
################################################################################
boxplotParaCalAllDiff <- function(boxplotStats, 
		defaultSstats, verbose=TRUE)
{

  defaultSname <- colnames(defaultSstats)

  if (verbose) cat("Start Difference calculation trought ", defaultSname,  "\n")
  nS<- dim(boxplotStats)[2]
  if (verbose) cat(paste("Calculate the Difference on:  ", nS, "samples\n"))

  if(length(defaultSname) < 2 ){
     samp1<- boxplotStats[c(2,3,4),1]
     defaultMd <- defaultSstats[c(2,3,4)]
     getDif<- getDifferenceBox(samp1[1], samp1[2], samp1[3], defaultMd[1], defaultMd[2], defaultMd[3], defaultSname)

      if(nS > 1){
        for(i in 2:nS){
          samp1<- boxplotStats[c(2,3,4),i]

          getDif<- rbind( getDif, getDifferenceBox(samp1[1], samp1[2], samp1[3], defaultMd[1], defaultMd[2], defaultMd[3], defaultSname ))

        }
      }
     return(getDif)

  }else{
      samp1<- boxplotStats[c(2,3,4),1]
      if (verbose) cat(paste("Calculate the Difference of Median-Mean:", defaultSstats[c(2,3,4)], "\n"))
      defaultMd <- defaultSstats[c(2,3,4),"median"]
      defaultMn <- defaultSstats[c(2,3,4),"mean"]
      getDif.Md<- getDifferenceBox(samp1[1], samp1[2], samp1[3], defaultMd[1], defaultMd[2], defaultMd[3],"median")
      getDif.Mn<- getDifferenceBox(samp1[1], samp1[2], samp1[3], defaultMn[1], defaultMn[2], defaultMn[3], "mean")

      if(nS > 1){
      for(i in 2:nS){
            samp1<- boxplotStats[c(2,3,4),i]

          getDif.Md <- rbind( getDif.Md, getDifferenceBox(samp1[1], samp1[2], samp1[3], defaultMd[1], defaultMd[2], defaultMd[3], "median"))
          getDif.Mn <- rbind( getDif.Mn, getDifferenceBox(samp1[1], samp1[2], samp1[3], defaultMn[1], defaultMn[2], defaultMn[3], "mean"))
      }
    }

  #to test with the both :median and mean
  getDifT<-NULL
  getDifT<-cbind(getDif.Md, getDif.Mn)
  return(getDifT)

  }

}


###############################################################
# aid function used in boxplotParaCalAllDiff function to calculate separatly differences
# hl1:Lower hinge ,md1: median  hu1:upper hinge of any boxplot
# hld:Lower hinge ,mdd: median  hud:upper hinge of default boxplot calculated with the function getMedianDefault
# hl: boxplot.stats$stats[2] , md: boxplot.stats$stats[3]  hu: boxplot.stats$stats[4]
###############################################################
getDifferenceBox <- function(hl1, 
		md1, hu1, 
		hld, mdd, 
		hud, nm)
{

  mdDif = abs(md1 - mdd)
  iqr1 =abs(hu1 - hl1)
  iqrd = hud - hld
  iqrDif = abs(iqr1 - iqrd)

  tDif <- cbind(mdDif, iqrDif)
  iqrnm <- paste("IQRDif_",nm, sep="")
  nmd <- paste("MedDif_", nm ,sep="")
  colnames(tDif)<-c(nmd, iqrnm)

  return (tDif)
  }


boxplotParagCrSamp <- function(differencen, 
		limits, iqrL, 
		verbose=TRUE)
{
     n<- dim(differencen)[2]
     l<-1
     nam<- colnames(limits)
 #to calcualte the Quality Problem Index under the first methode - whith percent

     if(iqrL){
        index <- vector("list",n+2)

       for(i in 1:n){
          index[l]<- list(which(differencen[,i] < limits[1,l]))
           l<-l+1
          index[l]<- list(which(differencen[,i] > limits[1,l]))
           l<- l+1
       }
        names(index) <- nam
     }else{
         index <- vector("list",n)
         for(i in 1:n){
            index[i]<- list(which(differencen[,i] >= limits[1,i]))
         }
         names(index) <- nam

     }
    return(index)
}

#########################################################################################################
# function to class the critical samples in 3 types :
# 1- Samples with med and IQR classfied as critical
# 2- Samples with only med classified as critical
# 3- Samples with only IQR classified as critical
# when the critical samples are classified based on the median or mean default Sample (for the differences Method)
##########################################################################################################
boxplotParacheckCritSamp <- function(index, 
		iqrL, verbose=TRUE)
{
 n <- length(index)
 # when the critical samples are classified based on the median or mean (for the differences Method)
 # media and mean should be considered
 criticTmd <- NULL
 criticTmn <- NULL
 # when used method to calculate the problem samples is the IQR-Method
 if(iqrL){
 #save the critical Samples as list from the colnames
 criticIQR<- unlist(index[grep("IQR_[a-zA-Z]+", names(index))])
 criticMd<- unlist(index[grep("med_[a-zA-Z]+", names(index))])
 criticTmd <- vector("list",3)
 #class the samples in the 3 types
 criticTmd[1] <- list(intersect(criticIQR, criticMd))
 criticTmd[2]<- list(setdiff(criticMd, criticIQR))
 criticTmd[3]<- list(setdiff(criticIQR, criticMd))
 #gives the name of the different classes
 names(criticTmd) <- c("crit.mdIQR", "crit.md",  "crit.IQR")

 }else{ # when used method to calculate the problem samples is the difference-Method

 if(length(index[[1]])>0 & length(index[[2]])> 0){
   criticTmd <- vector("list",3)
   criticTmd[1] <- list(intersect(index[[1]], index[[2]]))
   criticTmd[2] <- list(setdiff(index[[2]], index[[1]]))
   criticTmd[3] <- list(setdiff(index[[1]], index[[2]]))

   if(names(index[1]) == "MedDif_median" ){
     names(criticTmd) <- c("critMd.mdIQR", "critMd.md",  "critMd.IQR")
   }else{
     names(criticTmd) <- c("critMn.mdIQR", "critMn.md",  "crititMn.IQR")
   }
 }else if(length(index[[1]])>0  & length(index[[2]])== 0){
   criticTmd <- list(index[[1]])
    if(names(index[1]) == "MedDif_median" ){
     names(criticTmd) <-  "critMd.md"
    }else{
     names(criticTmd) <-  "critMn.md"
    }
 }else if(length(index[[1]])==0  & length(index[[2]])>0){
    criticTmd<- list(index[[2]])

    if(names(index[1]) == "MedDif_median" ){
      names(criticTmd) <- "critMd.md"
    }else{
     names(criticTmd) <- "critMn.mdIQR"
    }
 }else{return(NA)}
 # When default samples has been calculated from the median and mean
 if(n>3){
   if(length(index[[3]])>0 & length(index[[4]]) > 0) {
     criticTmn <-vector("list", 3)
     criticTmn[1] <- list(intersect(index[[3]], index[[4]]))
     criticTmn[2] <- list(setdiff(index[[3]], index[[4]]))
     criticTmn[3] <- list(setdiff(index[[4]], index[[3]]))
     names(criticTmn) <- c("critMn.mdIQR", "critMn.md",  "critMn.IQR")
   }else if (length(index[[3]])>0  & length(index[[4]])== 0){
     criticTmn<- list(index[[3]])
     names(criticTmn) <-  "critMd.md"
   }else if (length(index[[3]])==0  & length(index[[4]])>0){
     criticTmn<- list(index[[4]])
     names(criticTmn) <-  "critMd.mdIQR"
   }else{return(NA)}
  }

 }
 # to return the values of the critical samples if the median and mean are been used
 if(length(criticTmd)>0 & length(criticTmn)> 0 ){
     criticTotal <- list(criticTmd, criticTmn)
   }else if(length(criticTmd)>0 & length(criticTmn) == 0 ){ #when only media has been used
     criticTotal <- criticTmd
   }else{ criticTotal <- criticTmn} #when only mean has been used

 return(criticTotal)

}


#########################################################################################################
# function to calculate the histogram with the bevor calculated differences between samples and
# calculate default and the limit line (calculated from the percent) which indicate
# where are the problematic Samples
# parameter : differencen =matrix with the calculate differences at the  histMdWhisker function,
# percent= to be considered as limit between critical Samples and ok samples
########################################################################################################
boxplotParagLimits <- function(differencen, 
		percent, plotDrawn, 
		verbose=TRUE)
{

  n<- length(differencen[1,])
  nam<-colnames(differencen)[c(1:n)]
  #built the matrix to save the limit values from percent and total samples
  r<- matrix(rep(0, (n*3)),c(3,n))
  rownames(r) <- c(paste("limit values for", percent, "%"), paste(((1-percent)*100),"% Samples"), paste((percent*100),"% Samples"))
  if(plotDrawn){
  file_txt <- paste("Hist_",percent,"percent", sep="")
  op <- par(mfrow=c(2,n/2))
  }
  #loop to drawn the histogram from the calculated differences(2 or 4 columns)
  for(i in 1:n){
	  #sort differences as decreasing (greader - smaller) values
	  differencenS <-sort(differencen[,i], decreasing=T)
	
	  #calculate from porcent and total samples, the index value of the sample
	  # to be considered as limit between problematic samples and "normal" samples
	  ltH<- length(differencenS)
	  ltHp<-  ltH*percent
	  per.l <- ceiling(length(differencenS)*percent)
	  r[1,i] <- differencenS[per.l]
	  if (verbose) print(paste("limit value: ", r[1,i]))
	  #amount of values to be considered as "normal"
	  r[2,i] <- ltH - per.l
	  #amount of values to be considered as "normal"
	  r[3,i] <- per.l
	  #function to drawn the histogram
	  if(plotDrawn)
	   drawHistDiff(differencenS,r[1,i], nam[i])
 
  }
  #End of the drawHist, reset to previous settings
  if(plotDrawn) par(op)
  #saved the colnames of the matrix from the origin of the calculated differences: med, IQR
  colnames(r)<-nam
  return(r)
}


############################################################################
# aids function for boxplotParagLimits to drawn the histogram in every loop
###########################################################################
drawHistDiff <- function(differencenS, 
		r1i, nami)
{
  #formatted limits values with 2 digits
  r1if<- format(r1i, digits=2)
  #text to be wroten at the line
  line_txt <- paste( "critic value: ",r1if)
  hist(differencenS, breaks=12, main= nami, xlab="Absolut Intensities Differences\n(Samples-defaultSample)")
  #draw the limit line
  abline(v= r1i, col="red", lty=1)
  text(r1if, 6, line_txt, col ="red" )
 }

##############################################################################################
#Function to drawn the boxplotPara
#
##############################################################################################
boxplotParaDrawn <- function(boxpl.st,
		defaultS, qualityProblem,
		lim, nSample,
		plotAllBoxes=TRUE, verbose=TRUE)
{

  if (verbose) cat("\n\tPlot will be generated\n")

  if (verbose>1) cat("\tPlot ALL ? : ", plotAllBoxes,"\n")
  #When default is calculate for only one typDef or two (he median and mean)
  mNdef <- 0
  mDdef <- 0
  limmD <- 0
  limmN <- 0
  defaultSmD<- 0
  defaultSmN<- 0
  #when default has been calculated with median and mean
  if(length(defaultS) < 3){
    #save the median values of the defaultSample
    mDdef<- defaultS$median$stats[3,]
    mDHL <- defaultS$median$stats[2,]
    mDHU <- defaultS$median$stats[4,]
    defaultSmD<- defaultS$median
    #save the mean values of the defaultSample
    mNdef <- defaultS$mean$stats[3,]
    mNHL <- defaultS$mean$stats[2,]
    mNHU <- defaultS$mean$stats[4,]
    defaultSmN<- defaultS$mean
    # the limits values to be considered
    limmD<- lim[1, c(1,2)]
    limmN<- lim[1, c(3,4)]
  }else{ # when defaultSample only has been calculated as median or mean
    if(defaultS$names =="median"){
       mDdef  <- defaultS$stats[3,]
       mDHL<- defaultS$stats[2,]
       mDHU<- defaultS$stats[4,]
       limmD<- lim[1,]
       defaultSmD<- defaultS

    }else{
       mNdef <-defaultS$stats[3,]
       mNHL <- defaultS$stats[2,]
       mNHU <- defaultS$stats[4,]
       limmN<- lim[1,]
       defaultSmN<- defaultS
   }
  }
  #Indexes list of the Quality Problem Sanples
  uQP<- unlist(qualityProblem)
  namesuQP<- names(uQP)
  #number of samples to be ploted
  toBoxpl <- NULL
  #number of plots necessary to plot all Samples
  nplot <- 1
  #when parameter plotAllBoxes =FALSE we need to drawn a boxpl.st whith only the Quality problem Samples
  if(!plotAllBoxes){
    #get only the bad quality Problem samples
    toBoxpl <- getTotalBoxToDrawn(boxpl.st, uQP)
    nbxp <- length(toBoxpl$stats[1,])
    tmpnames <- toBoxpl$names
    cols <- rep("lightblue",nbxp )
    cols[grep("[a-zA-Z]+[.]mdIQR[1-9]+",namesuQP)] <- "red"
    cols[grep("[a-zA-Z]+[.]md[1-9]+",namesuQP)] <- "orange"
    cols[grep("[a-zA-Z]+[.]IQR[1-9]+",namesuQP)] <- "pink"
    
    if (verbose) cat(nbxp, " Samples will be used to create the boxplot\n")
  }else{ # When all samples should be ploted
    toBoxpl <- boxpl.st
    nbxp<- length(boxpl.st$stats[1,])
    tmpnames <- toBoxpl$names
    #to gives the color and the names of the qualityProblem Sample
    nam <- rep("", nbxp)
    cols <- rep("lightblue", nbxp)
    #class 1: color=red, Samples as medIQR
    nam[uQP[grep("[a-zA-Z]+[.]mdIQR[1-9]+",namesuQP)]] <-uQP[grep("[a-zA-Z]+[.]mdIQR[1-9]+",namesuQP)]
    cols[uQP[grep("[a-zA-Z]+[.]mdIQR[1-9]+",namesuQP)]] <- "red"
    #class 2: color=orange, Samples as med
    cols[uQP[grep("[a-zA-Z]+[.]md[1-9]+",namesuQP)]] <- "orange"
    nam[uQP[grep("[a-zA-Z]+[.]md[1-9]+",namesuQP)]] <- uQP[grep("[a-zA-Z]+[.]md[1-9]+",namesuQP)]
    #class 3: color= pink, Samples as IQR
    cols[uQP[(grep("[a-zA-Z]+[.]IQR[1-9]+",namesuQP))]] <- "pink"
    nam[uQP[grep("[a-zA-Z]+[.]IQR[1-9]+",namesuQP)]] <- uQP[grep("[a-zA-Z]+[.]IQR[1-9]+",namesuQP)]
    #remplace the names of the boxplot$name with only the names of the critical Samples
    toBoxpl$names <- nam
  }
 #to test if  are necessary more as one boxplot(nplot)
 if(nbxp > nSample){
   necBoxplot <- getNumberPlots(nbxp, nSample)
   nplot<- necBoxplot[1]
   nSample <- necBoxplot[2]
  }else{
   nSample = nbxp
   nplot <- 1
  }
 if (verbose) cat( paste(nSample, " Samples will be used to create ", nplot, "figure with boxplots\n"))
 #the boxplot are separated in median and mean to drawn the defaultSample

 #when the media for the DefaultSample has been calculated
 if(!mDdef == 0){
   nS <- 1
   nE <- nSample
   file_txt <- paste( "boxplot_", nS, "-", nE, "Samples",sep="")
   main_txt <- paste("Median", file_txt, sep="")
   #paper =size of paper (a4r : A4 landscape) , onefile= all generated boxplots are saved as one File
  #loop to drawn the boxplots
  for (n in 1:nplot){
      bxp.plot <- boxplotSplit(toBoxpl,nS, nE) 
      drawnBxp(bxp.plot, defaultSmD, cols, main_txt, mDdef, mDHL, mDHU, "median","darkviolet")

      #nS and nE values for the next loop
      nS = nS+nSample
      nE = nE +nSample
      #check the End sample at set it as last sample when greader than number of total samples to be ploted
      if(nE > nbxp) nE = nbxp

   }
   
  }

 #when the media for the DefaultSample has been calculated
 if(!mNdef == 0){
  nS <- 1
  nE <- nSample
  file_txt <- paste( "boxplot_", nS, "-", nE, "Samples",sep="")
  main_txt <- paste("Mean", file_txt, sep="")
  #loop to drawn the boxplots
  for (n in 1:nplot){
   bxp.plot <- boxplotSplit(toBoxpl,nS, nE) 
   drawnBxp(bxp.plot, defaultSmN, cols, main_txt, mNdef, mNHL, mNHU, "mean","green")
   #nS and nE values for the next loop
   nS = nS+nSample
   nE = nE +nSample
   #check the End sample at set it as last sample when greader than number of total samples to be ploted
   if(nE > nbxp) nE = nbxp
   }
  }

 #set back the correct names for the boxplot$names
 toBoxpl$names <- tmpnames

}


##############################################################################################
#Function to split automatically the Samples in plots with a total from n Samples when total Samples > 250
#parameter:  StatsObject = boxplot.stats, nSamplesS = start box, nSamplesE = Ende box
##############################################################################################
boxplotSplit <- function(StatsObject,
		nSamplesS, nSamplesE)
{
  #create a new stats list with the samples to be ploted at once
  bxp.plot<- vector("list", 6)
  names(bxp.plot)<-c("stats", "n","conf", "out", "group", "names")
    bxp.plot$stats <- matrix(StatsObject$stats[,c(nSamplesS:nSamplesE)],ncol=nSamplesE, nrow=5)
    bxp.plot$n <- StatsObject$n[c(nSamplesS:nSamplesE)]
    bxp.plot$conf <-  matrix(StatsObject$conf[,c(nSamplesS:nSamplesE)], ncol=nSamplesE, nrow=2)

    if(length(StatsObject$out) == 0) bxp.plot$out <- StatsObject$out[c(nSamplesS:nSamplesE)]
    bxp.plot$out <- NULL
    if(length(StatsObject$group) == 0)   bxp.plot$group <- StatsObject$group[c(nSamplesS:nSamplesE)]
    bxp.plot$group <- NULL
    bxp.plot$names <- StatsObject$names[c(nSamplesS:nSamplesE)]

    return(bxp.plot)
}


################################################################################
# Function to get out  only the problematic quality box, only them should be ploted
# when parameter plotAllBoxes= FALSE
####################################################################################
getTotalBoxToDrawn <- function(boxpl.st, 
		uqp)
{
  n <- length(uqp)
  bxp.plot<- vector("list", 6)
  names(bxp.plot)<-c("stats", "n","conf", "out", "group", "names")
  bxp.plot$stats <- boxpl.st$stats[,uqp]
  bxp.plot$n <- boxpl.st$n[c(1:n)]
  bxp.plot$conf <- boxpl.st$conf[,uqp]

  bxp.plot$out <- 0
  if (!is.null(boxpl.st$out) & length(boxpl.st$out)>0 ){ bxp.plot$out <- boxpl.st$out[,uqp]}

  bxp.plot$group <- 0
  if(!is.null(boxpl.st$group) & length(boxpl.st$group)>0) bxp.plot$group <- boxpl.st$group[,uqp]
  tempName<- uqp[[1]]
  for(i in 2:n){
  tempName<- c(tempName, uqp[[i]])
  }
  bxp.plot$names <- tempName
  return(bxp.plot)
}

#######################################################################################
# Function to dranw the boxplot with default Sample and limit lines as different colors
######################################################################################
drawnBxp <- function(bxp.plot, 
		defaultS, cols, 
		main_txt, mdef, 
		mHL,mHU, 
		typ, colTyp ){
 
 #boxplot for all samples
 bxp(bxp.plot, boxfill=cols, main=main_txt, lwd=0.5, medcol="white", varwidth=TRUE, xlab="Probe Arrays from the E-GEOD-7123 Data", ylab="Unprocesed log scale probe intensities", las=2)
 #boxplot for the defaultSample
 bxp(defaultS, add=TRUE, at=0.8, boxfill="blue", axes=FALSE, show.names=FALSE)
 #limit lines from the defaultSample
 abline(h= mdef, lty=1, lwd=2, col=colTyp)
 abline(h= mHL, lty=2, lwd=2, col=colTyp)
 abline(h= mHU, lty=2, lwd=2, col=colTyp)
 #legend
 legend(x="topright", inset=0.005, bty="0", bg="gray90", cex=0.8, title = "Color explanation", legend= c("Bad Quality Sample", "Bad Quality only Media", "Bad Quality only IQR", "Good Quality Sample", "calculated default Sample", typ), fill = c("red", "orange", "pink", "lightblue", "blue", colTyp))

}


#############################################################################
# function to calculate how many plots at the same File (loops) are necessary
# to plot all samples. Parameter:
# nSample:numeric, as parameter in the function boxplotPara.Indicates how many samples
# should be ploted at once.
############################################################################
getNumberPlots <- function(nbxp, 
		nSample)
{
  #calculate the number of boxplots
  if(nSample > 0 & nSample < nbxp) {
    nplot <- ceiling(nbxp / nSample)
    lastPlot <- nbxp %% nSample
    if(lastPlot > 0 & lastPlot <= 175 & nplot >= 1){
      nplot<- nplot-1
      nSample <- nSample + ceiling(lastPlot/nplot)
    }
   }else{
     nplot <- 1
     nSample <- nbxp
   }
  #to start variable whith the first and the last sample to be ploted
  return(c(nplot, nSample))

}

#############################################################################
# function to calculate the median and IQR from all Samples
# calculated the limits and get the bad quality Samples
# at the last drawn the boxplot with these values
#
#############################################################################
boxplotParagMedIQR <- function(IQRmedMatrix, 
		plot=TRUE, verbose=TRUE)
{
 #calculate IQR
   iqrSam<- abs(IQRmedMatrix[,1] - IQRmedMatrix[,2])

 # IQR critical values
    iqr.st <- boxplot(iqrSam, plot=FALSE)
    iqrCriticl <- iqr.st$stats[1]
    iqrCriticu <- iqr.st$stats[5]

 #Calculated critical values for median from all Samples
    med.st <- boxplot(IQRmedMatrix[,3], plot=FALSE)
    medCriticl <- med.st$stats[1]
    medCriticu <- med.st$stats[5]
    medMedian <- med.st$stats[3]

 #first boxplot with the medians from all samples
    iqrnam <- "IQR_Samples"
    mednam <- "Med_Samples"
    if (plot){ 
    op<- par(mfrow=c(1,2))
    #First boxplot: median
    bxp(med.st, main="Medians Boxplot from all Samples", ylab="log2(Median) from Intensities")
    abline(h= medCriticl, lty=2, lwd=1, col="red")
    abline(h= medCriticu, lty=2, lwd=1, col="red")

    #Second boxplot with the IQR from all Samples
    bxp(iqr.st, main="IQRs Boxplot from all Samples", ylab="log2(IQR) from Intensities")
    abline(h= iqrCriticl, lty=2, lwd=1, col="red")
    abline(h= iqrCriticu, lty=2, lwd=1, col="red")
    #End of the bxp, reset to previous settings
    par(op)
    }
    #return a matrix with the limits to calculate the Quality Problem Samples
    boxplotCrT <- matrix(c(iqrCriticl, iqrCriticu, medCriticl, medCriticu),nrow=1, ncol=4)
    colnames(boxplotCrT) <- c("IQR_HL", "IQR_HU", "med_HL", "med_HU" )
    rownames(boxplotCrT) <- c("limit values")
    return(boxplotCrT)

}

############################################################################
#
# Necessary functions to check the Parameter of the boxplotPara function:
#  typDef, percent, nSample
############################################################################

########################################### #########################################
# function to check the parameter typDef: how the default Sample should be calculated
####################################################################################
 boxplotParaCheckTypDef<- function(typDCheck)
 {
    cht<- TRUE
   if(length(typDCheck) == 1){
       if(!is.character(typDCheck) | !(typDCheck == "median" | typDCheck == "mean"))
            cht<- FALSE
       else
          cht <- TRUE
   }else if(length(typDCheck) == 2) {
            if(boxplotParaCheckTypDef(typDCheck[1]))
                boxplotParaCheckTypDef(typDCheck[2])
            else cht <- FALSE
   }else
       cht<-FALSE

}

###########################################
# function to check the parameter percent
##########################################
boxplotParaCheckPercent <- function(perc)
{

   if(!is.numeric(perc))
     chp <- FALSE
   else{

      if(trunc(perc) == 0 & (perc <= 0 | perc >=1 ))
        chp<- FALSE

      else if(trunc(perc)> 0)
        chp<- FALSE

      else
        chp <- TRUE
  }

    return(perc)

}

###########################################
# function to check the parameter nSample
##########################################
boxplotParaChecknSamples <- function(nS, 
		lnAffyB)
{

  if(!is.numeric(nS) | nS < 1 | nS >200){ 
     chn<- FALSE            
  }else if ((nS > 0) &  (nS/trunc(nS)> 1)){  
     chn<- FALSE    
  }else{
    chn <-TRUE
  }  
 }
 
 
###########################################
# function to check the parameter percent
##########################################
boxplotParaCheckPercentnSample <-  function(perc, 
		lobject)
{
      if(trunc(perc) == 0 & perc*lobject < 0.5){
        if(lobject <50) chps <- 2
        else chps<- FALSE
      } else if (trunc(perc) > 0 & perc/100*lobject < 0.5){
        if(lobject <50) chps <- 2
        else chps <- FALSE
      }else
        chps <- TRUE

}

