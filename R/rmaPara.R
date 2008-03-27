# rmaPara.R
#
# Parallelization of the RMA Preprocessing Function
#
# History
# 22.02.2008 : Version 0.1
# 27.03.2008 : Version 0.2 - object.type as input removed
#
# Copyright (C) 2008 : Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de>
###############################################################################

rmaPara <- function(cluster,
				object,
				ids = NULL,
				phenoData = new("AnnotatedDataFrame"), cdfname = NULL,
				verbose=FALSE)
{
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