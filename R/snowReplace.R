# snowReplace.R
# 
# functions for replacing snow functions to remove cluster object
# 
# History:
# 17.12.2008 : Version 0.1
#
#  Copyright (C) 2008 : Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de>
###############################################################################

makeCluster <- function(...)
{
	if( exists("cl", envir=.affyParaInternalEnv) ){
		stop("In one R master instance you can only run ONE cluster!\nThere is already running another snow cluster!")	
	}else{
		#makeCluster
		cl <- snow::makeCluster(...)
		#save cl in internal environment
		assign("cl", cl, envir=.affyParaInternalEnv)
	}
}
		
		
stopCluster <- function(cl)
{
	#Get cluster object
	if(missing(cl)){
		if( !exists("cl", envir=.affyParaInternalEnv) )
			stop("There is no running snow cluster!")
		else
			cl <- .affyParaInternalEnv$cl
	}
	checkCluster(cl)
	#stopCluster
	snow::stopCluster(cl)
	#remove cl variable
	if( exists("cl", envir=.affyParaInternalEnv) )
		rm("cl", envir=.affyParaInternalEnv)
}
