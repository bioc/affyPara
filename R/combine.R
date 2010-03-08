# combine.R
# 
# Functions to combine objects (in lists) after distributed computing
#
# History
# 29.04.2009 : ... old stuff removed ...
# 10.12.2007 : Version 0.5 - #Remove NA's from list added
# 14.12.2007 : Version 0.6 - Input from combineMatrices changed to matrices.list
# 07.07.2008 : Version 0.7 - mergeAffyBatches improved for great Matrices
# 27.08.2008 : Version 0.8 - mergeAffyBatches bugfix in annotations
# 23.03.2009 : Version 0.9 - Option verbose set to getOption("verbose") and added . to 
#							  names of internatl functions
# 29.04.2009 : Version 0.10 - mergeAffyBatch improved for memory amd time
#
# Copyright (C) 2008 - 2010 : Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de>
###############################################################################

############################################
# Function to merge AffyBatches
mergeAffyBatches <- function (abatch.list,
    description = NULL, notes = character(0))
{
	if( length(abatch.list) == 1){
		#if there is only on AffyBatch in the list, no merge required!
		return(abatch.list[[1]])
	}else{
		#Remove non AffyBatch objects (NA's) from list (if more slaves in universum then used)
		omit <- lapply(abatch.list, function(x){
					if (class(x)=="AffyBatch"){
						return(TRUE)
					}else return(FALSE)
				})
		abatch.list <- abatch.list[unlist(omit)]

	    #check affyBatches
		anrow <- nrow(abatch.list[[1]])
		ancol <- ncol(abatch.list[[1]])
		acdfName <- cdfName(abatch.list[[1]])
	    if ( any(unlist(lapply(abatch.list,nrow)) != anrow) 
	      || any(unlist(lapply(abatch.list,ncol)) != ancol) )
	   		   stop("cannot merge chips of different sizes !")
	   	if ( any(unlist(lapply(abatch.list,cdfName)) != acdfName) )
	        warning("cdfName mismatch (using the cdfName of x)!")
		
		#data from affyBatch
		ab <- abatch.list[[1]]
		abatch.list[[1]]<-NA
		gc(verbose=FALSE)
		phenodata <- phenoData(ab)
		exprs <- intensity(ab)
		
	    #Check description
	    if (is.null(description)) {
	        description <- new("MIAME")
	        description@title <- "Created from merging several AffyBatches. No description was supplied. The description of the two original AffyBatches was not kept."
	    }

	    #Check notes
	    if (length(notes) == 0){
	      notes(description) <- paste("Merge of", length(abatch.list),"AffyBatches.")
	    } else
			notes(description) <- notes


		#Generat Exprs-Matrix and phenodata
		for( k in 2:length(abatch.list)){
			pData(phenodata) <- rbind(pData(phenodata), pData(abatch.list[[k]]))
			exprs <- cbind(exprs,intensity(abatch.list[[k]]))
			abatch.list[[k]]<-NA
			gc(verbose=FALSE)
		}
		
	    #Return ONE AffyBatch
		exprs(ab) <- exprs
		phenoData(ab) <- phenodata
		experimentData(ab) <- description
	    return(ab)
	}
}

############################################
# Function to combine Matrices
combineMatrices <- function(matrix.list, verbose=getOption("verbose"))
{
	if (verbose) cat("Rebuild matrix ")
  	#Dimension of matrix
  	nrow <- nrow(matrix.list[[1]])
  	if (is.null(nrow)) nrow <- length(matrix.list[[1]])

  	#Rebuilding
  	return( matrix(unlist(matrix.list),nrow=nrow) )
	if (verbose) cat(" DONE\n")
}
