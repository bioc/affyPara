# rmaPara.R
#
# Parallelization of the RMA Preprocessing Function
#
# History
# 22.02.2008 : Version 0.1
# 27.03.2008 : Version 0.2 - object.type as input removed
# 17.12.2008 : Version 0.3 - cluster object gets default parameter
# 18.12.2008 : Version 0.4 - cluster object gets default parameter: .affyParaInternalEnv$cl
# 20.03.2009 : Version 0.5 - Bug Fix in proproPara call
# 23.03.2009 : Version 0.6 - Option verbose set to getOption("verbose") and added . to names of internatl functions
# 16.07.2009 : Version 0.7 - summary.method added to function parameters "medianpolish_orig"
#
# Copyright (C) 2009 : Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de>
###############################################################################

rmaPara <- function(object, ids = NULL,
				phenoData = new("AnnotatedDataFrame"), cdfname = NULL,
				cluster, verbose=getOption("verbose"), summary.method="medianpolish_orig")
{
	#Get cluster object form default environment
	if(missing(cluster))
		cluster <- .affyParaInternalEnv$cl
	
	preproPara(object,
				bgcorrect=TRUE, bgcorrect.method="rma",
				normalize=TRUE, normalize.method="quantiles", normalize.param=list(type="pmonly"), 
				pmcorrect.method="pmonly", 
				summary.method=summary.method,
				ids = ids,
				phenoData = phenoData, cdfname = cdfname,
				cluster)
}