
alter.pop <- function(pop, merge=NULL, remove=NULL, keep=NULL  ){
	
	if(length(pop)==1){
		if(!inherits(pop,"pop")){
			stop("must provide object of class 'pop' ")
		}
		
		theta <- pop@theta
		
		if(sum(dim(pop@r))>0){
			r <- pop@r
		}else{
			r <- matrix(rep(NA, dim(theta)[1]*dim(theta)[2]), dim(theta)[1], dim(theta)[2])
			colnames(r)=colnames(theta); rownames(r)=rownames(r)
		}			
		
		if(class(pop)=="anypop"| class(pop)=="bipop"){
			geno <- pop@geno
				
		}else{
			geno <- matrix(rep(NA, dim(theta)[1]*dim(theta)[2]), dim(theta)[1], dim(theta)[2])
			colnames(geno)=colnames(theta);rownames(geno)=rownames(theta)
			}
		
		
	}
	
	# combine pops
	if(length(pop) >1){
				
		f <- length(pop)
		
		markers <- unique(unlist(lapply(pop, function(p) rownames(p@theta))))
		gid <- unlist(lapply(pop, function(p) colnames(p@theta)))
		m <- length(markers)
		n <- length(gid)
		
		theta <- matrix(rep(NA, m*n), m, n)
		rownames(theta) <- markers; colnames(theta) <- gid
		
		geno <- theta
		r <- theta
			
		# combine data frames
		for(i in 1:f){
			theta[rownames(pop[[i]]@theta),colnames(pop[[i]]@theta)] <- pop[[i]]@theta
			
			if(sum(dim(pop[[i]]@r))>0){
				r[rownames(pop[[i]]@r), colnames(pop[[i]]@r)] <- pop[[i]]@r
			}
			
			if(class(pop[[i]])=="anypop"| class(pop[[i]])=="bipop"){
				geno[rownames(pop[[i]]@geno), colnames(pop[[i]]@geno)] <- pop[[i]]@geno
			}
			
		}
	
		
		if( length(gid) > length(unique(gid))){
			# consensus for theta is median, for geno {1, 1, NA}=1, {1, NA, NA}=1, {1, 1, 2}=1, {1, 2, 3}= NA
				
			tab <- table(gid)
			merge2 <- lapply(names(tab)[which(tab >1)], function(x) which(gid==x)) # ! tab sorts to alphabetical order
			
			drop <- unlist(merge2)
			#warning("Samples with identical sample names printed to the console.")
			options(warn=1)
			print("Non-unique samples merged: ")
			print(unique(gid[drop]))
		
			mtheta <- sapply(1:length(merge2), function(x) {.consensus(theta, merge2[[x]])})
			colnames(mtheta) <- sort(unique(gid[drop])); rownames(mtheta) <- markers
			theta <- cbind(theta[ ,-drop], mtheta)	
				
			mr <- sapply(1:length(merge2), function(x) {.consensus(r, merge2[[x]])})
			colnames(mr) <- sort(unique(gid[drop])); rownames(mr) <- markers
			r <- cbind(r[ ,-drop], mr)
		
			if(class(pop[[i]])=="anypop"| class(pop[[i]])=="bipop"){
				mgeno <- sapply(1:length(merge2), function(x) {.consensusgeno(geno, merge2[[x]])})
				colnames(mgeno) <- sort(unique(gid[drop])); rownames(mgeno) <- markers
				geno <- cbind(geno[ ,-drop], mgeno)
					
			}	
		
		}
	}
	
	# merge samples
	if(!is.null(merge)){
		
		if(!is.list(merge)){
			stop(" 'merge' must be a list")
		}
		
		if(length(names(merge))!= length(merge)){
			stop("the 'merge' list must have names for the consensus samples to create")
		}
		
		discard <- unname(unlist(merge))
		
		if(length(setdiff(discard, colnames(theta)))>0){
			stop("one or more samples indicated for merge are not present in the data")
		}
		
		mtheta2 <- sapply(1:length(merge), function(x) {.consensus(theta, merge[[x]])})
		colnames(mtheta2) <- names(merge)
		theta <- cbind(theta[, -which(colnames(theta) %in% discard)], mtheta2)
		
		if(sum(dim(pop@r))>0){
			mr2 <- sapply(1:length(merge), function(x) {.consensus(r, merge[[x]])})
			colnames(mr2) <- names(merge)
			r <- cbind(r[, -which(colnames(r) %in% discard)], mr2)
		}
		
		if(!all(is.na(geno))){
			
			mgeno <- sapply(1:length(merge), function(x) {.consensusgeno(geno, merge[[x]])})
			colnames(mgeno) <- names(merge)
			geno <- cbind(geno[, -which(colnames(geno) %in% discard)], mgeno)
	
			
		}

	}
	
	
	# keep / remove
	if(!is.null(keep) & !is.null(remove)){
		stop("only keep or remove may be indicated")
	}
	
	if(!is.null(keep)){
		
		if(length(setdiff(keep, colnames(theta))) >0){
			stop("one or more samples indicated to keep are not in the object")
		}
		
		theta <- theta[ ,keep]
		
		r <- r[ ,keep]

		
		if(!is.null(geno)){
			geno <- geno[, keep]
		}
		
	}
	
	if(!is.null(remove)){
		
		if(length(setdiff(remove, colnames(theta))) >0){
			stop("one or more samples indicated to remove are not in the object")
		}

		theta <- theta[ , -which(colnames(theta) %in% remove)]

		r <- r[ ,-which(colnames(r) %in% remove)]

		
		if(!is.null(geno)){
			geno <- geno[, -which(colnames(geno) %in% remove)]
		}

		
	}	
	
	
	
	
	
	
	
	if(length(pop) >1 ){
		if(all(is.na(geno))){
			if(!all(is.na(r))){
				return(new("pop", theta=theta, r=r))
			}else{
				return(new("pop", theta=theta))
			}	
		}else{
			bipops <- unlist(lapply(pop, function(x) inherits(x, "bipop")))
			pops2 <- lapply(pop, function(x) which(colnames(theta) %in% colnames(x@theta)))
			pops2[bipops] <- lapply(pop[bipops], function(x) c(which(colnames(theta) %in% c(x@parent1, x@parent2)), 
				which(colnames(theta) %in% setdiff(colnames(x@theta), c(x@parent1, x@parent2)))))
			
			names(pops2)[bipops] <- lapply(pop[bipops], function(x) 
			if(!is.na(x@parent1)){
				paste0(substr(x@parent1, 1, 1), "x", substr(x@parent2, 1, 1))	
			}else{
				paste0("diverse")
			}
			)
			names(pops2)[!bipops] <- rep("diverse", length(!bipops))
			
			warning("concordance is not calculated for the returned object of class anypop")
			return(new("anypop", theta=theta, geno=geno, r=r, pops=pops2))
		}	
	}else{
		
		if(all(is.na(geno))){
			if(!all(is.na(r))){
				return( new("pop", theta=theta, r=r))
			}else{
				return( new("pop", theta=theta))
			}
		}else{
			if(class(pop)=="bipop"){
				if(!is.null(merge) | !is.null(remove) | !is.null(keep)){
					warning("altering the sample composition of the bipop object may render some stored marker statistics inaccurate")
					
					if(length(intersect(c(pop@parent1, pop@parent2), colnames(theta)))==2){
						if(!all(is.na(r))){
							return(new("bipop", geno=geno, theta=theta, r=r, info=pop@info, parent1=pop@parent1, parent2=pop@parent2))
						}else{
							return(new("bipop", geno=geno, theta=theta, info=pop@info, parent1=pop@parent1, parent2=pop@parent2))
						}
					}else{
						warning("parental samples have been altered, an object of class anypop is returned")
						if(!all(is.na(r))){
							return(new("anypop", geno=geno, theta=theta, r=r, pops=list(1:dim(geno)[2]) ))
						}else{
							return(new("anypop", geno=geno, theta=theta, pops=list(1:dim(geno)[2]) ))
						}
					}		
					
				}
			}	
			
			if(class(pop)=="anypop"){	
					warning("altering the sample composition of the anypop object may render the concordance statistic inaccurate")			

					if(!all(is.na(r))){
						return(new("anypop", geno=geno, theta=theta, r=r, concordance=pop@concordance, pops=list(1:dim(geno)[2])))
					}else{
						return(new("anypop", geno=geno, theta=theta, concordance=pop@concordance, pops=list(1:dim(geno)[2])))
					}	
			}else{
				if(!all(is.na(r))){
					return(new("pop", theta=theta, r=r))
				}else{
					return(new("pop", theta=theta))
				}	
			}	
		
		}
	
			
	}	

}











