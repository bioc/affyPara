# zzz.R
# 
# functions for initialize and clean up.
#
# History:
# 28.11.2007 : Version 0.1
# 29.02.2008 : Versiom 0.2 : .onLoad added
#
# Copyright (C) 2008 : Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de>
###############################################################################

.onLoad <- function (libname, pkgname){
	require(affy)
	require(snow)
}

.onAttach <- function(lib, pkg) {
	if(interactive() && .Platform$OS.type == "windows" && .Platform$GUI == "Rgui"){
		addVigs2WinMenu("affyPara")
	}
}

.onUnload <- function( libpath ) {
	library.dynam.unload("affyPara", libpath)
}
