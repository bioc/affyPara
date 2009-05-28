# file: qa.R
# 
# Summarization for parallel quality assessment
#
# History
# 28.05.2009 : Version 0.1 - added to package
#
# Copyright (C) 2009 : Esmeralda Vicedo <e.vicedo@gmx.net>, Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de> 
###############################################################################


summaryM1M2Para <- function(method1, method2, 
		level, verbose=FALSE)
{
	#Check for affy
	require(affyPara)
	#check the mehode if boxplot or  Maplot
	method1Sum<- .getSumResults(method1,level, 1,verbose)
	method2Sum<- .getSumResults(method2,level, 2, verbose)
	
	methodSum<- cbind(method1Sum,method2Sum[,c(2:(level))])
	return (methodSum)
	
}

#########################################################################
.getSumResults<- function(method,level,N,verbose){
	SumMatrix <- NULL
	if(verbose) print(paste("Method", method , "Level ", level))  
	if(is.element("results_boxP",names(method))){
		if(verbose) print("boxplotPara")                          
		results <- method$results_boxP
		namResults <- colnames(results)
		level1 <- NULL
		level2 <- NULL
		level3 <- NULL
		
		if(level==3){
			level1 <- namResults[grep("[a-zA-Z]*.mdIQR", namResults)]
			level2 <- namResults[grep("[a-zA-Z]*[.]md$", namResults)] 
			level3 <- namResults[grep("[a-zA-Z]*[.]IQR", namResults)]     
		}else if(level==2){
			level1 <- "crit.mdIQR"
			level2 <- "crit.md"
			
		}else{
			level1 <- "crit.mdIQR"  
		}
		
		if(verbose) print(paste("level1", level1 , "Level2 ", level2, "level3", level3, "N", N))
		SumMatrix <- .createSumMatrix(results,level,level1, level2,level3,N, verbose)
		
	}else{
		if(verbose) print("MAplotPara")    
		results <-  method$results_MAP
		# possible names of the colums if method is  MAplotPara
		level1 <- "firstLevel"
		level2 <- c("s-loess", "s-sigma", "loess-sigma")
		level3 <- c("s", "loess", "sigma")
		if(verbose) print(paste("level1", level1 , "Level2 ", level2, "level3", level3, "N", N))
		SumMatrix <- .createSumMatrix(results,level,level1, level2,level3,N, verbose)
	}
	
	return (SumMatrix)
}

#####################################################################
.createSumMatrix<- function(methodResults,level,level1,level2, level3, N, verbose){ 
	# number of necessary Colums for matrix
	nCol<-(level)
	
	nSamples<-dim(methodResults)[1]
	if(verbose) cat("number of samples to be analysed", nSamples, "\n")
	
	#matrix wiht levels for both methods
	sumBoxMA <- matrix(c(rep(0, nCol*nSamples)), nrow=nSamples, ncol=nCol)
	colnames(sumBoxMA) <- c( paste("L1_M",N,sep=""), paste("L2_M",N,sep=""), paste("L3_M",N,sep=""))
	
	
	indexNameFirstL <- which(is.element(colnames(methodResults), level1))
	indexNameSecondL <- which(is.element(colnames(methodResults), level2))
	indexNameThirdL <- which(is.element(colnames(methodResults), level3))
	
	
	# to fill the matrix with the correct values for the differents values
	if(length(indexNameFirstL) > 0) indexFL <- .indexValueL(indexNameFirstL,methodResults,1)
	else indexFL <- 0
	if(length(indexNameSecondL) > 0) indexSL <- .indexValueL(indexNameSecondL,methodResults,1)
	else indexSL <- 0
	if(length(indexNameThirdL) > 0) indexTL <- .indexValueL(indexNameThirdL, methodResults,1)
	else indexTL <- 0 
	# fill the SumBoxMatrix with the correct values
	#sumBoxMA[,1]<- methodResults[,"sampleNames"]
	rownames(sumBoxMA)<- methodResults[,"sampleNames"]
	sumBoxMA[indexFL,1]<-1
	sumBoxMA[indexSL,2]<-1
	sumBoxMA[indexTL,3]<-1
	
	return(sumBoxMA)
	
	
}

################################################################
.indexValueL <- function(indexN, results, value){
	tmpIndex<- NULL
	for(i in 1: length(indexN))    
		tmpIndex <- sort(unique(c(tmpIndex,which(results[,indexN[i]]==value))))
} 