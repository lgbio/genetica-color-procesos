
CC.bipop <- function (pop, parent1, parent2, max.missing=0.05, min.sep=0.05, max.range=0.25, hom.theta=0.15, min.posterior=0.7, min.r=0.2, mono.thresh=0.95,chi2p.step=1,n.core=1){

	fz <- function(ix,theta,r,max.missing,max.range,min.posterior,hom.theta,min.sep,min.r,mono.thresh,chi2p.step) {
		ans <- apply(array(ix), 1, function(i){.biparental.call(theta=theta[i,],r=r[i,],max.missing=max.missing,max.range=max.range,min.posterior=min.posterior,hom.theta=hom.theta,min.sep=min.sep,min.r=min.r,mono.thresh=mono.thresh,chi2p.step=chi2p.step)})
		return(ans)
	}
	
	if(!inherits(pop,"pop")){
		stop("must provide object of class 'pop'; if error checking was used reading the data in, use popname[[1]]")
	}
	
	if(!(parent1 %in% colnames(pop@theta))){
		stop("parent1 must be included in the 'pop' ")
	}
	if(!(parent2%in% colnames(pop@theta))){
		stop("parent2 must be included in the 'pop' ")
	}
	
	par <- c(parent1, parent2)
	in.theta <- pop@theta
	order <- c(par, colnames(in.theta)[-which(colnames(in.theta) %in% par)])
	theta <- in.theta[,order]
	
	if(sum(dim(pop@r))>0){
		in.r <- pop@r
		r <- in.r[,order]
	} else {
		r <- NULL
	}
	
	m <- dim(theta)[1]
	n <- dim(theta)[2]

	if ((n.core > 1)&requireNamespace("parallel",quietly=TRUE)) {
		it <- split(1:m,factor(cut(1:m,n.core,labels=FALSE)))
		ans <- unlist(parallel::mclapply(it,fz,theta=theta,r=r,max.missing=max.missing,max.range=max.range,min.posterior=min.posterior,hom.theta=hom.theta,min.sep=min.sep,min.r=min.r,mono.thresh=mono.thresh,chi2p.step=chi2p.step,mc.cores=n.core),recursive=F)
	}else{
		ans <- fz(1:m,theta=theta,r=r,max.missing=max.missing,max.range=max.range,min.posterior=min.posterior,hom.theta=hom.theta,min.sep=min.sep,min.r=min.r,mono.thresh=mono.thresh,chi2p.step=chi2p.step)				
	}
	geno <- t(sapply(ans,function(x){x$geno}))
	colnames(geno) <- colnames(theta)
	rownames(geno) <- rownames(theta)
	chi2p <- sapply(ans,function(x){x$chi2p})
	called <- sapply(ans,function(x){x$called})
	mkr.data <- data.frame(marker=rownames(theta), called=called, negLog10p= -log10(chi2p), stringsAsFactors=F)
	
	if(is.null(r)){
		return(new("bipop", theta=theta[, !duplicated(colnames(theta))], geno=geno[, !duplicated(colnames(geno))], info=mkr.data, parent1=parent1, parent2=parent2))
	}else{
		return(new("bipop", theta=theta[, !duplicated(colnames(theta))], geno=geno[, !duplicated(colnames(geno))], r=r[, !duplicated(colnames(r))], info=mkr.data, parent1=parent1, parent2=parent2))	
	}
}








