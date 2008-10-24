# vsnPara.R
#
# Parallelization of the VSN function
#
# History
# 10.06.2008 : Version 0.1 - first ideas
# xx.09.2008 : Version 0.2 - implementation
# 14.10.2008 : Version 0.3 - first running version
# 21.10.2008 : Version 0.4 - everything vectorized and strata removed 
# 24.10.2008 : Version 0.5 - subsamples no work correct
#
# TODO: vsnrma
# TODO: reference
#
# Sending AffyBatch form master to slave an back is very time consuming. Sending a list
# of CEL files from master to slave, creating the AffyBatch and do normalization is faster.
# Using the right combination "size of AffyBatch on slaves" - "number of slaves" the parallelized
# version is more than ten times faster as the serial version.
#
# Copyright (C) 2008 : Markus Schmidberger <schmidb@ibe.med.uni-muenchen.de>
###############################################################################

#we parallelize vsn2, because 'vsn' has been superseded by 'vsn2'.
justvsnPara <- function(cluster,
	object,
	phenoData = new("AnnotatedDataFrame"), cdfname = NULL,
	subsample, ...,
	verbose=TRUE) 
{	
    ########
    # Checks
    ########
	#Check for affy amd snow
	require(affy)
	require(snow)
	
	#Check cluster and generate number.parts
	checkCluster(cluster)
	number.parts <- length(cluster)
	
	#Check object type
	object.type <- getObjectType(object) 
	
	#Check size of partitions
	parts <- checkPartSize(object, number.parts)
	number.parts <- parts$number.parts
	object.length <- parts$object.length
	
	####################
	#Partition of object
	###################
	if (verbose) cat("Partition of object ")
		t0 <- proc.time();
		if (object.type == "AffyBatch"){
			object.list <- splitAffyBatch(object, number.parts)
			samples.names <- sampleNames(object)
		} else if( object.type == "CELfileVec" ){
			object.list <- splitFileVector(object, number.parts)
			samples.names <- gsub("^/?([^/]*/)*", "", unlist(object), extended = TRUE)
		} else if( object.type == "partCELfileList" ){
			object.list <- object
			object <- unlist(object)
			object.length <- length(object)
			samples.names <- gsub("^/?([^/]*/)*", "", unlist(object), extended = TRUE)
		}				
		t1 <- proc.time();
	if (verbose) cat(round(t1[3]-t0[3],3),"sec DONE\n")			
	
	#Info-Output for Distribution
	if (verbose) cat("Object Distribution:", unlist(lapply(object.list,length)), "\n") 
	
	#################################
	#Initialize AffyBatches at slaves
	##################################
	if (verbose) cat("Initialize AffyBatches at slaves ")
		t0 <- proc.time();
		#initialize affybatch
		dimAB <- clusterApply(cluster, object.list, initAffyBatchSF, object.type, rm.all=FALSE)
		#initialize intensity matrix and remove AB
		check <- clusterCall(cluster, setIntMatSF, rm.AB=FALSE)
		#calculate data distribution
		dist<-unlist(lapply(object.list,length))
		dist_list<-list()
		for(i in 1:length(dist)){
			if(i==1){
				dist_list[[1]]=1:dist[1]
				end <- dist_list[[i]][length(dist_list[[i]])]
			}else{
				dist_list[[i]]=(1+end):(end+dist[i])
				end <- dist_list[[i]][length(dist_list[[i]])]	
			}	
		}
		clusterApply(cluster, dist_list, function(x)assign("dist", value=x, envir= .GlobalEnv) )
	t1 <- proc.time();
	if (verbose) cat(round(t1[3]-t0[3],3),"sec DONE\n")
	
	#################################
	#Normalization depending on types
	#################################
	dimAB <- unlist(dimAB)
	ncol <- sum (dimAB[seq(2,length(dimAB),2)], na.rm=TRUE)
	nrow <- dimAB[1]
	dimAB <- c(nrow,ncol)
	
	vsn2Para(cluster, dimAB, subsample, ... , verbose=verbose)
	
	##############################
	#Combine / Rebuild affyBatches
	##############################
	if (verbose) cat("Rebuild AffyBatch ")
	t0 <- proc.time();
		AffyBatch.list.norm <- clusterCall(cluster, getAffyBatchSF)
		AffyBatch <- mergeAffyBatches(AffyBatch.list.norm)
	t1 <- proc.time();
	if (verbose) cat(round(t1[3]-t0[3],3),"sec DONE\n")
	
	#################
	#Return AffyBatch
	#################
	return(AffyBatch[,samples.names])
}

####################################################################################
## vsn2 paralleliziert
####################################################################################
vsn2Para <- function(cluster,
		dimAB, reference,
		subsample, ...,
		verbose=TRUE)
{
	nr <- dimAB[1]
	#TODO wollen wir die 30000 wirklich?
	if(missing(subsample)){
		if( nr > 30000L )
			subsample = 30000L
		else
			subsample = 0L
	}
	vsnMatrixPara(cluster,
			dimAB, reference,  
			subsample=subsample,
			defaultpar = list(factr=5e7, pgtol=2e-4, maxit=60000L, trace=0L, cvg.niter=4L, cvg.eps=0), 
			..., verbose=verbose)
}

####################################################################################
## vsnMatrix paralleliziert
####################################################################################
vsnMatrixPara <- function(cluster,
		dimAB,
		reference,
		lts.quantile = 0.9,
		subsample    = 0L,
		verbose      = interactive(),
		returnData   = TRUE,
		pstart,
		minDataPointsPerStratum = 42L,
		optimpar   = list(),
		defaultpar = list(factr=5e7, pgtol=2e-4, maxit=60000L, trace=0L, cvg.niter=7L, cvg.eps=0))
{
	nr <- dimAB[1]
	nc <- dimAB[2]
	
	if(missing(pstart))
		pstart = vsn:::pstartHeuristic(matrix(ncol=nc), list(all=seq_len(nr)))
	
	if(missing(reference)) {
		if(nc<=1L)
			stop("'x' needs to have 2 or more columns if no 'reference' is specified.")
		reference = new("vsn")
	}# else {
	#	if(nrow(reference)!=nrow(x))
	#		stop("'nrow(reference)' must be equal to 'nrow(x)'.")
	#	if(nrow(reference)!=length(reference@mu))
	#		stop(sprintf("The slot 'reference@mu' has length %d, but expected is n=%d", nrow(reference)))
	#}
	
	if(!(is.list(optimpar)&&all(names(optimpar)%in%vsn:::optimparNames)&&!any(duplicated(names(optimpar)))))
		stop(paste("Names of elements of 'optimpar' must be from: ",
						paste("'", vsn:::optimparNames, "'", collapse=", ", sep=""), ".", sep=""))
	opar = defaultpar
	opar[names(optimpar)] = optimpar
	
	v = new("vsnInputPara",
			dimAB  = dimAB,
			pstart = pstart,
			reference = reference,
			lts.quantile = lts.quantile,
			optimpar = opar,
			subsample = subsample,
			verbose   = verbose,
			ordered   = FALSE)
	
	## Print welcome message
	if (verbose)
		cat("vsn: ", nr, " x ", nc, " matrix. ", sep="")
	
	## Here, the actual work is done.
	res = vsnSamplePara(cluster, v, verbose)
	
	## hoffset: compute from average scale factors between arrays, but separately for the different strata.
	## coefficients 'cof' from are taken from the reference, if available.
	cof = coefficients( if(nrow(reference)==0L) res else reference )
	res@hoffset = log2(2*scalingFactorTransformation(rowMeans(cof[,,2L,drop=FALSE])))
	
	## calculate the data matrix transformed according to 'coefficients'
	vsn2trsfPara(cluster, coefficients(res), hoffset = res@hoffset, wh=FALSE, reset=TRUE)

	stopifnot(validObject(res))

	if(verbose) {
		cat("Please use 'meanSdPlot' to verify the fit.\n")
	}
	return(res)
}

####################################################################################
## vsnSample parallelized: allows for (balanced) subsetting / subsample
####################################################################################
vsnSamplePara <- function(cluster,
		v, verbose=TRUE)
{
	wh = NULL
  
  	if(v@subsample>0L) {
    	wh = sample(nrow(v), size=v@subsample)	
  	}

	if(!is.null(wh)){
		# reduce x to subsample dim / save subsample id
		check <- clusterCall(cluster, function(wh){
					assign("wh", value=wh, envir= .GlobalEnv)
					return(TRUE)
				}, wh)
		
		tmpDim <- v@dimAB[1] 
		v@dimAB[1] <- length(wh) 
		#call next vsn step
		res = vsnStrataPara(cluster, v, verbose)
		## put back the results from subsampling
		v@dimAB[1] <- tmpDim
		newmu = rep(NA_real_, nrow(v))
		newmu[wh] = res@mu
		res@mu = newmu
	} else {
		res = vsnStrataPara(cluster, v, verbose)
	}
	
	return(res)
}

######################################################################################
## vsnStrata parallelized: if necessary,
## reorder the rows of the matrix so that each stratum sits in a contiguous block
####################################################################################
vsnStrataPara <- function(cluster,
		v, verbose=TRUE)
{
  v@ordered = TRUE

  res =  vsnLTSPara(cluster, v, verbose)
  
  return(res)
} 

################################################################################
## vsnLTS : parallelized robust modification of the ML estimator
################################################################################
vsnLTSPara <- function(cluster,
		v, verbose)
{
	## for calculating a convergence criterion: earlier result
	oldhy   = Inf
	## the number of iterations that have already met the convergence criterion
	cvgcCnt = 0L
	
	for(iter in seq_len(v@optimpar$cvg.niter)) {
		if(v@verbose)
			vsn:::progress(iter-1L, v@optimpar$cvg.niter)
		
		sv  = if(iter==1L) v else{
					#v[whsel, ]
					check <- clusterCall(cluster, function(whsel){
								assign("whsel", value=whsel, envir= .GlobalEnv)
								return(TRUE)
							}, whsel)
					sv <- v
				}
		rsv = vsnMLPara(cluster, sv, verbose)
		
		## if LTS.quantile is 1, then the following stuff is not necessary
		if (vsn:::isSmall(v@lts.quantile-1))
			break
		
		## apply to all data
		vsn2trsfPara(cluster, coefficients(rsv), wh=TRUE)
		v@pstart = coefficients(rsv)
		
		## Calculate residuals
		if(length(v@reference@mu)>0L) {
			## with reference:
			hmean = v@reference@mu
		} else {
			## without reference:
			hmean <- rowMeansPara("hy",ncol(v))
			## double check with what was computed in vsnML:
			if(iter==1L) {
				if(rsv@lbfgsb==0L) stopifnot(vsn:::isSmall(rsv@mu-hmean))
			} else {
				if(rsv@lbfgsb==0L) stopifnot(vsn:::isSmall(rsv@mu-hmean[whsel]))
				## and create rsv@mu with mu of the right length (NA for the 'outliers' not in whsel)
				tmp = rep(NA_real_, nrow(v))
				tmp[whsel] = rsv@mu
				rsv@mu = tmp
			}
		}
		
		## row variances
		rvar <- rowVPara("hy", hmean)
		
		## select those data points whose rvar is within the quantile; do this separately
		## within each stratum, and also within strata defined by hmean
		## (see the SAGMB 2003 paper for details)
		nrslice = 5
		slice   = ceiling(rank(hmean, na.last=TRUE)/length(hmean)*nrslice)
		slice = factor(slice)
		grmed   = tapply(rvar, slice, quantile, probs=v@lts.quantile, na.rm=TRUE)
		if(any(is.na(grmed)))
			stop(sprintf("Too many data points are NA (%d of %d), not enough data for fitting, am giving up.",
							sum(is.na(sv@x)), length(sv@x)))
		
		meds    = grmed[as.character(slice)]
		whsel   = which(rvar <= meds)

        #convergence check
		#TODO convCheck bisher nicht implementiert
	}
	
	if(v@verbose){
		vsn:::progress(1L, 1L)
		cat("\n")
	}
	
	return(rsv)
}

################################################################################
# vsn2_trsf parallelized: for parallelized transformation
################################################################################
vsn2trsfPara<- function(cluster,
		par, hoffset=NULL, 
		wh=TRUE, reset=FALSE)
{
	check <- clusterCall(cluster, vsn2trsfParaSF, par, hoffset, wh, reset)			
}

vsn2trsfParaSF <- function(par, hoffset, wh, reset)
{
	# vectorized
	if (exists("dist", envir = .GlobalEnv) &&
			exists("x", envir = .GlobalEnv)) {
		
		dist <- get("dist", envir = .GlobalEnv)
		x <- get("x", envir = .GlobalEnv)
		if (exists("wh", envir = .GlobalEnv) && wh==TRUE){
			wh <- get("wh", envir = .GlobalEnv)
			x <- x[wh,]
			x <- as.matrix(x)
		}
		
		# calculate parameter
		a <- par[1:(length(par)/2)] #first part
		b <- par[(1+(length(par)/2)):length(par)] #second part
		a <- a[dist]
		b <- b[dist]
		
		nr <- dim(x)[1]
		nc <- dim(x)[2]
		
		#calculate transformation
		hy = asinh( matrix(exp(b), nrow=nr, ncol=nc, byrow=TRUE) * x + matrix(a, nrow=nr, ncol=nc, byrow=TRUE) ) 
		if(!is.null(hoffset)){
			hy <- hy / log(2) - hoffset						
		}
		
		#reset AffyBatch
		if(reset){
			if (exists("AffyBatch", envir = .GlobalEnv)){
				require(affy)
				AffyBatch <- get("AffyBatch", envir = .GlobalEnv)
				intensity(AffyBatch) <- hy
				assign("AffyBatch", value=AffyBatch, envir= .GlobalEnv)
			}else
				return(FALSE)
		}else
			assign("hy", value=hy, envir= .GlobalEnv)
		return(TRUE)
	}else
		return(FALSE)
}

################################################################################
## vsnML parallelized estimator
################################################################################
vsnMLPara <- function(cluster,
	v, verbose=TRUE) 
{
	p = as.vector(v@pstart)
	
	o = vsn2_optimPara(cluster, v@dimAB, p, v@optimpar, verbose)#istrat
	
	rv = new("vsn", coefficients=o$coefficients, 
			mu=o$mu, sigsq=o$sigsq, lbfgsb=o$fail, hoffset=rep(NA_real_,dim(o$coefficients)[1]))
	
	if (o$fail!=0L) {
		nrp = if(length(p)<=6L) length(p) else 6L
		msg = paste("** This is a diagnostic message on the performance of the optimizer,\n",
				"** it need not indicate a relevant problem:\n",
				sprintf("fail=%d\npstart[1:%d]=", o$fail, nrp),
				paste(signif(p[1:nrp], 4L), collapse=", "),
				sprintf("\n  coef[1:%d]=", nrp),
				paste(signif(coefficients(rv)[1:nrp], 4L), collapse=", "), "\n", sep="")
		## warning(msg)
	}
	
	return(rv)
}

################################################################################
# vsn2_optim parallelized: does optimazation
################################################################################
vsn2_optimPara <- function(cluster,
		SdimAB, Spar, 
		Soptimpar, verbose=TRUE)
{
	#set up Parameter
	lmm <- 5
	fail <- 0
	nREPORT <- 10
	
	factr <- Soptimpar$factr
	pgtol <- Soptimpar$pgtol
	maxit <- Soptimpar$maxit
	trace <- Soptimpar$trace
	
	fcount <- 0
	grcount <- 0
	
	px <- setupPara(cluster, Spar, SdimAB )
	cpar <- Spar
	
	lower <- vector(length=px$npar)
	upper <- vector(length=px$npar)

	lower[1:(length(cpar)/2)]<- -Inf
	upper[1:(length(cpar)/2)]<- Inf
	lower[(1+(length(cpar)/2)):length(cpar)]<- -100
	upper[(1+(length(cpar)/2)):length(cpar)]<- 100
	
	#Optimize
	erg <- optim(cpar, logikPara, grad_loglikPara, px, cluster,verbose,
			method = "L-BFGS-B", lower = lower, upper = upper,
			control = list(trace=trace, REPORT=nREPORT, lmm=lmm, pgtol=pgtol, factr=factr, maxit=maxit),
			hessian = FALSE) 
	px <- get("px", envir = .GlobalEnv)
	
	#set up output
	vfail <- erg$convergence
	dimcoef <- vector(length=3)
	dimcoef[1] <- px$npar / (px$ncol*2)
	dimcoef[2] <- px$ncol
	dimcoef[3] <- 2
	
	res <- list(fail=vfail, coefficients=array(erg$par, dimcoef), sigsq=px$sigsq, mu=px$mu)
	return(res)
}

##############################################################################
# This setup function is used by vsn2_optimPara
# It sets up workspaces px                   
##############################################################################
setupPara <- function(cluster,
		Spar, SdimAB)
{
	#set PX at Master
	px<-list()
	px$npar<-length(Spar)
	px$nrow<-SdimAB[1]
	px$ncol<-SdimAB[2]
	px$sigsq <- NA
	px$mu <- vector(length=px$nrow)
	
	#set PX at Nodes
	check <- clusterCall(cluster, setupParaSF, Spar)#, Sstrat)
	
	return(px)
}

setupParaSF <- function(Spar)
{
	#set PX at Nodes
	if (exists("x", envir = .GlobalEnv)) {
		px <- list()
		px$npar <- length(Spar)
		x <- get("x", envir = .GlobalEnv)
		if (exists("wh", envir = .GlobalEnv)){
			wh <- get("wh", envir = .GlobalEnv)
			x <- x[wh,]
			x <- as.matrix(x)
		}	
		if (exists("whsel", envir = .GlobalEnv)){
			whsel <- get("whsel", envir = .GlobalEnv)
			x <- x[whsel,]
			x <- as.matrix(x)
		}	
		px$y <- x 
		px$nrow <- dim(x)[1]
		px$ncol <- dim(x)[2]
		px$ntot <- length(which(is.na(x)))
	
		px$mu <- vector(length=px$nrow)
		px$profiling <-TRUE

		px$ly <- matrix(nrow=px$nrow, ncol=px$ncol)
		px$asly <- matrix(nrow=px$nrow, ncol=px$ncol)
		px$resid <- matrix(nrow=px$nrow, ncol=px$ncol)
		px$ma <- matrix(nrow=px$nrow, ncol=px$ncol)
		
		assign("px", value=px, envir= .GlobalEnv)
				
		return(TRUE)
	}else
		return(NA)
}

################################################################################
# parallel calculation of the function to be optimized:
#a are the offsets, b the factors. 
################################################################################
logikPara <- function(par,
		px, cluster, 
		verbose=TRUE)
{
	#######
	# 1st sweep through the data: compute Y_ki, h(y_ki), A_ki, B_ki 
	######
	jac <- clusterCall(cluster, logikParaSF1, par)
	jac1 <- sum( unlist(jac)[seq(1,length(unlist(jac)),3)] )
	jac2 <- sum( unlist(jac)[seq(2,length(unlist(jac)),3)] )
	nt <- sum( unlist(jac)[seq(3,length(unlist(jac)),3)] )
	jacobian = jac1 * 0.5 - jac2
	
	#######
	# 2nd sweep through the data: compute r_ki                      
	#######
	px$mu <-rowMeansPara("px", px$ncol, slot="asly")
		
	#vectorized
	ssq_list <- clusterCall(cluster, logikParaSF2, px$mu)
	ssq <- sum( unlist(ssq_list) )
	sigsq <- ssq/nt
	px$sigsq <- sigsq
	residuals <- nt/2
	
	scale = nt*0.5 * log(2* pi * sigsq)
	ll = scale + residuals + jacobian
	
	#output for debuging
	if(verbose>1){
		cat("\n\tll:",ll)
		cat(" | Scale: ", scale)
		cat(" | Res: ", residuals)
		cat(" | Jac1:", jac1)
		cat(" | Jac2:", jac2)
		cat(" | nt:", nt)
		cat(" | sigsq:", sigsq)
		cat(" | ssq/nt", ssq/nt)
	}
	cat("\n\tll:",ll)
	
	px <<- px
	return(ll)	
}

logikParaSF1 <- function(par)
{
	#vectorized!
	if (exists("px", envir = .GlobalEnv) &&
		exists("dist", envir = .GlobalEnv) ) {
			dist <- get("dist", envir = .GlobalEnv)
			px <- get("px", envir = .GlobalEnv)
			
			#calculate parameter
			a <- par[1:(length(par)/2)] # first part
			b <- par[(1+(length(par)/2)):length(par)] # second part
			a <- a[dist]
			b <- b[dist]
	
			# calculate values
			z <- px$y
			z <- matrix(exp(b), nrow=px$nrow, ncol=px$ncol, byrow=TRUE) * z +  matrix(a, nrow=px$nrow, ncol=px$ncol, byrow=TRUE)
			px$ly <- z
			px$asly <- asinh(z)
			px$ma <- 1/sqrt(1+z*z)
			jac1 <- sum(colSums(log(1+z*z), na.rm=TRUE))
			ni <- colSums(!is.na(z), na.rm=TRUE)
			jac2 <- sum(ni * b)
			nt <- sum(ni)	
			
			assign("px", value=px, envir= .GlobalEnv)
			return(list(jac1=jac1,jac2=jac2,nt=nt))
	}else
		return(NA)
}

logikParaSF2 <- function(mu)
{
	#vectorized
	if (exists("px", envir = .GlobalEnv)) {
		px <- get("px", envir = .GlobalEnv)
		px$mu <- mu		
		z <- px$asly - mu
		ssq <- sum(rowSums(z*z, na.rm=TRUE))
		px$resid=z
		assign("px", value=px, envir= .GlobalEnv)
		return(ssq)
	}else
		return(NA)
}

################################################################################
# parallel gradient calculation
################################################################################
grad_loglikPara <- function(par, px, cluster, verbose)
{
	px <- get("px", env=.GlobalEnv)
	rfac <- 1/px$sigsq
	
	#vectorized gradient calculation at nodes
	gr_list <- clusterCall(cluster, grad_loglikParaSF, rfac, par)
	gr1 <- vector(length=0)
	gr2 <- vector(length=0)
	for(i in 1:length(gr_list)){
		gr1 <- c(gr1,gr_list[[i]]$gr1)
		gr2 <- c(gr2,gr_list[[i]]$gr2)
	}	
	gr<- c(gr1,gr2)

	#output for debugging
	if(verbose>1){
		cat("\n\tpar:",par)
		cat("\n\tgr:",gr,"\n")
	}

	return(gr)
}

grad_loglikParaSF <- function(rfac, par)
{
	if (exists("px", envir = .GlobalEnv) &&
			exists("dist", envir = .GlobalEnv) ) {
		px <- get("px", envir = .GlobalEnv)
		dist <- get("dist", envir = .GlobalEnv)
		
		#calculate parameter 
		b <- par[(1+(length(par)/2)):length(par)] # second part
		b <- b[dist]
		
		#calculate gradients
		z = px$resid*rfac + px$ma * px$ly
		sa = colSums(z*px$ma, na.rm=TRUE)
		sb = colSums(z * px$ma * px$y, na.rm=TRUE)
		nj = colSums(!is.na(z))
		gr1 = sa
		gr2 = exp(b) * (sb - nj/exp(b))
		
		assign("px", value=px, envir= .GlobalEnv)
		return(list(gr1=gr1,gr2=gr2))
	}else
		return(NA)
}

################################################################################
## Class vsnInputPara for parallel vsn
################################################################################
setClass("vsnInputPara",
		representation(
				dimAB  = "vector",     ## The dimension of the n*d data matrix
				reference = "vsn", ## A result from a previous fit (for reference normalization)          
				ordered = "logical",  ## Are the levels consecutive in "strata"?               
				lts.quantile = "numeric",
				subsample = "integer",
				verbose   = "logical",
				pstart    = "array",     ## Start parameters: see comment on slot 'coefficients'
				## in definition of class 'vsn'
				optimpar  = "list"),     ## See below: optimparnames
		prototype = list(
				dimAB = vector(length=2),
				reference = new("vsn"),
				ordered = FALSE,
				lts.quantile = 1,
				subsample = 0L,
				verbose = TRUE,
				pstart = array(as.numeric(NA), dim=c(1L,0L,2L)),
				optimpar = list(factr=5e7, pgtol=2e-4, 
						maxit=60000L, trace=0L, cvg.niter=7L, cvg.eps=0)),
)
		#validity = vsn:::validityVsnInput)

setMethod("nrow", signature("vsnInputPara"), function(x) x@dimAB[1] )
setMethod("ncol", signature("vsnInputPara"), function(x) x@dimAB[2])
setMethod("dim",  signature("vsnInputPara"), function(x) x@dimAB)
