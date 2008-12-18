# rmaPara.R
#
# Parallelization of the RMA Preprocessing Function
#
# History
# 22.02.2008 : Version 0.1
# 27.03.2008 : Version 0.2 - object.type as input removed
# 17.12.2008 : Version 0.3 - cluster object gets default parameter
# 18.12.2008 : Version 0.4 - cluster object gets default parameter: .affyParaInternalEnv$cl
#
# Copyright (C) 2008 : Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de>
###############################################################################

rmaPara <- function(object, ids = NULL,
				phenoData = new("AnnotatedDataFrame"), cdfname = NULL,
				cluster, verbose=FALSE)
{
	#Get cluster object form default environment
	if(missing(cluster))
		cluster <- .affyParaInternalEnv$cl
	
	preproPara(cluster,
				object,
				bgcorrect=TRUE, bgcorrect.method="rma",
				normalize=TRUE, normalize.method="quantiles", normalize.param=list(type="pmonly"), 
				pmcorrect.method="pmonly", 
				summary.method="medianpolish",
				ids = ids,
				phenoData = phenoData, cdfname = cdfname,
				verbose=verbose)
}