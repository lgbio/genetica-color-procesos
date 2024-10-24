
read.pop <- function(theta.file, r.file=NULL, error.checking=FALSE, method="complete", thresh=5, max.missing.theta=0.2, ...){
	
	theta <- as.matrix(read.csv(file=theta.file, header=T, as.is=T, check.names=F, row.names=1))
	thresh=thresh # dendrogram height threshold for calling individuals identical
	
	if(dim(theta)[2]!=length(unique(colnames(theta)))){

			dup <- colnames(theta)[which(duplicated(colnames(theta))==TRUE)]
			if(length(dup) ==1){
			stop(paste0("all samples must have unique sample names- ", dup, " appears more than once"))
			}else{
			stop(paste0("all samples must have unique sample names- ", paste(dup, collapse=" , "), " appear more than once") )
			}

	}
	
	if(!missing(r.file)){
		r <- as.matrix(read.csv(file=r.file, header=T, as.is=T, check.names=F, row.names=1))
	
		if(dim(theta)[[1]] != dim(r)[[1]] || dim(theta)[[2]] != dim(r)[[2]]){
			stop("theta and r matrices must have the same dimensions")
		}
	
		if(colnames(theta) != colnames(r) || rownames(theta) != rownames(r)){
			stop("theta and r matrices differ in content")
		}
		
		
	}
	
	
	# remove samples with >max.missing.theta values NA
	ix <- which(apply(theta, 2, function(x) sum(is.na(x))/dim(theta)[1]) > max.missing.theta)
	if(length(ix) >0){
		options(warn=1)
		warning("samples removed for > max.missing.theta proportion of missing theta values: ")
		print(names(ix))
		theta <- theta[, -ix]
		if(!missing(r.file)){
			r <- r[, -ix]
		}
	}
	
	if(error.checking==FALSE){
		
		if(!missing(r.file)){
			return(new("pop", theta=theta, r=r))
		}else{
			return(new("pop", theta=theta))
		}
		
	}else{
		plot(hclust(dist(t(theta)), method=method), cex=.75, ...)
		abline(h=thresh, col="red", lty=3)
		
		result <- .find.matching.sets(theta, thresh, method)
		
		if(!missing(r.file)){
			return(list(pop=new("pop", theta=theta, r=r), equivalent=result, removed=names(ix)))
		}else{
			return(list(pop=new("pop", theta=theta), equivalent=result, removed=names(ix)))
		}
		
	}

}