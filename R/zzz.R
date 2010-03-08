# zzz.R
# 
# functions for initialize and clean up.
#
# History:
# 28.11.2007 : Version 0.1
# 29.02.2008 : Version 0.2 : .onLoad added
# 17.12.2008 : Version 0.3 : affyParaInternalEnv added
# 24.03.2009 : Version 0.4 - Summarization optimized and . added to internal functions
#
# Copyright (C) 2008 - 2010 : Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de>
###############################################################################

.onLoad <- function (lib, pkg){
	require(affy)
	require(snow)
	
	##a place to store some variables that need to be accessed
	.affyParaInternalEnv <- new.env(parent=emptyenv())
	assign(".affyParaInternalEnv", .affyParaInternalEnv, envir=topenv(parent.frame()))
	
	## update summary methods
	upDate.express.summary.stat.methods(c(express.summary.stat.methods(), 'medianpolish_orig', 'liwong_orig', 'farms_orig', 'playerout_orig'))
}

.onAttach <- function(lib, pkg) {
	if(interactive() && .Platform$OS.type == "windows" && .Platform$GUI == "Rgui"){
		addVigs2WinMenu("affyPara")
	}
}

.onUnload <- function(libpath ) {
	library.dynam.unload("affyPara", libpath)
}
