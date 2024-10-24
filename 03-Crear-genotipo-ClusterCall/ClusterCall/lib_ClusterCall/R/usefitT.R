
usefitT <- function(theta, params=NULL, n.core=1){
	# if(!requireNamespace("fitTetra", quietly=TRUE)&requireNamespace("reshape", quietly=TRUE)){
		# stop("the R packages fitTetra and reshape must be installed to use this function")
	# }

	file <- "tmp"
	if(is.null(params)){
		params=list(sd.threshold= 0.1, p.threshold=0.9, call.threshold=0.6, peak.threshold=0.85, try.HW=TRUE, dip.filter=1)
	}else{
		if(!all(names(params) %in% c("sd.threshold", "p.threshold", "call.threshold", "peak.threshold", "try.HW", "dip.filter"))){
			stop("one or more parameters in the list params is not supported")
		}
	}

	dir.create("fitTetra_running")
	setwd("fitTetra_running")


	if((n.core > 1)&requireNamespace("parallel",quietly=TRUE)) {
		it <- split(1:dim(theta)[1], factor(cut(1:dim(theta)[1], n.core, labels=FALSE)))

		fz <- function(ix, theta, file, params){
  		lapply(array(ix), function(x){ .runfitT_write(theta=theta[x,,drop=FALSE], file=file, params=params)})
		}
		
		parallel::mclapply(it, fz, theta=theta, file=file, params=params)  
			
	}else{
		.runfitT_write(theta=theta, file=file, params=params) 
	}
	
	
	.runfitT_read(theta, file)
	scores <- as.matrix(read.csv(paste0(file, "_fitTonTheta.csv"), check.names=F, row.names=1))
	
	setwd("..")
	unlink("fitTetra_running", recursive=T, force=T)
	return(.make.bipop(scores, theta))
	
}

