# combine.R
# 
# Functions to combine objects (in lists) after distributed computing
#
# History
# 22.08.2007 : Version 0.1
# 06.09.2007 : Version 0.2 - combineAffyBatch -> merge.AffyBatches : much faster
# 27.11.2007 : Version 0.3 - merge.AffyBatches renamed to mergeAffyBatches
# 06.12.2007 : Version 0.4 - code cleaning and dokumentation
# 10.12.2007 : Version 0.5 - #Remove NA's from list added
# 14.12.2007 : Version 0.6 - Input from combineMatrices changed to matrices.list
#
# Copyright (C) 2008 : Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de>
###############################################################################

mergeAffyBatches <- function (abatch.list,
    description = NULL, notes = character(0))
{
    #number auf probes in affyBatch
    adim <- dim(intensity(abatch.list[[1]]))[1]
	
	#Remove NA's from list (more slaves in universum then used)
	omit <- lapply(abatch.list, function(x){
				if (class(x)=="AffyBatch"){
					return(TRUE)
				}else return(FALSE)
			})
	abatch.list <- abatch.list[unlist(omit)]
	
    #check affyBatches
    if ( any(unlist(lapply(abatch.list,nrow)) != nrow(abatch.list[[1]])) 
      || any(unlist(lapply(abatch.list,ncol)) != ncol(abatch.list[[1]])) )
   		   stop("cannot merge chips of different sizes !")
   	if ( any(unlist(lapply(abatch.list,cdfName)) != cdfName(abatch.list[[1]])))
        warning("cdfName mismatch (using the cdfName of x)!")
    #Check description
    if (is.null(description)) {
        description <- new("MIAME")
        description@title <- "Created from merging several AffyBatches. No description was supplied. The description of the two original AffyBatches was not kept."
    }
    #Create phenodata
    phenodata <- phenoData(abatch.list[[1]])
    pData(phenodata) <- data.frame(sample=unlist(lapply(abatch.list,pData)),row.names=unlist(lapply(abatch.list,names<-function(data){
    										#row.names(pData(data))
    										sampleNames(data)
    										})))
    #Check notes
    if (length(notes) == 0){
      notes(description) <- list(paste("Merge of", length(abatch.list),"AffyBatches with notes:",lapply(abatch.list,data<-function(data){
    				notes(experimentData(data))
    				})))
    }
    else
      notes(description) <- notes
  
  	#Generat Exprs-Matrix
	exprs <- matrix( unlist(lapply(abatch.list,intensity)), nrow=ncol(abatch.list[[1]])*nrow(abatch.list[[1]]) )
  
    #Return ONE AffyBatch
    return(new("AffyBatch", exprs = exprs, 
        phenoData = phenodata, experimentData = description, 
        cdfName = cdfName(abatch.list[[1]]), nrow = nrow(abatch.list[[1]]), ncol = ncol(abatch.list[[1]]), 
        annotation = abatch.list[[1]]@annotation))
}

combineMatrices <- function(matrix.list, verbose=FALSE)
{
	if (verbose) cat("Rebuild matrix ")
  	#Dimension of matrix
  	nrow <- nrow(matrix.list[[1]])
  	if (is.null(nrow)) nrow <- length(matrix.list[[1]])

  	#Rebuilding
  	return( matrix(unlist(matrix.list),nrow=nrow) )
	if (verbose) cat(" DONE\n")
}
