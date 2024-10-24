
CC.anypop <- function(train, predict=NULL,n.core=1, min.train=NULL, max.missing=0.05,max.range=0.25,min.posterior=0.7,min.sep=0.05,hom.theta=0.15, impute=TRUE){
	
		for (f in 1:length(train)){
			if( !inherits(train[[f]], "bipop"))
			stop("training populations must be of class bipop")
		}
				
		if(!is.null(min.train)){
			if(!(2 <= min.train & min.train <= length(train)))
			stop("min.train must be greater than 2 and fewer than the number of training populations")		
		}
		
					
		fz <- function(ix,TP,p.theta,TP.to.pred) {
			if (is.null(min.train)) {
				ans <- apply(array(ix),1,function(i){
					return(.theta2geno(train=data.frame(theta=TP$theta[i,],geno=TP$geno[i,]),predict=p.theta[i,],max.missing=max.missing,max.range=max.range,min.posterior=min.posterior,min.sep=min.sep,hom.theta=hom.theta,impute=impute))
					})
			} else {
				ans <- apply(array(ix),1,function(i){
					iv <- which(TP.to.pred[i,]==FALSE)
					train <- data.frame(theta=TP$theta[i,iv],geno=TP$geno[i,iv])
					iv <- which(TP.to.pred[i,]==TRUE)
					if (length(iv) > 0){
						predict <- c(p.theta[i,],TP$theta[i,iv])
					} else {
						predict <- p.theta[i,]
					}
					ans2 <- .theta2geno(train=train,predict=predict,max.missing=max.missing, max.range=max.range, min.posterior=min.posterior,min.sep=min.sep,hom.theta=hom.theta,impute=impute)
					return(ans2)
					})
			}
			return(ans)
		}
		
		n.fam <- length(train)
		# calls for markers called in at least min.train families and polymorphic in at least one
		present <- unique(unlist(lapply(train,function(x){x@info$marker})) )
		called <- unlist(lapply(train,function(x){x@info$marker[x@info$called]}))
		poly <- unlist(lapply(train,function(x){x@info$marker[!is.na(x@info$negLog10p)]}))
		marks <- unique(poly)
		tab1 <- table(called)
		
		if (!is.null(min.train)){
			train.m <- intersect(marks, names(tab1)[which(tab1 >= min.train)])
		}else{
			# require them to be called in all where present and poly in at least one
			tab2 <- table(unlist(lapply(train,function(x){x@info$marker})) )
			train.m <- marks[tab1[marks]==tab2[marks]]	
		}
		
		if(!is.null(predict)){
			train.m <- intersect(train.m, rownames(predict@theta))
		}
		
		m <- length(train.m)

		
		# read.pop and CC.bipop only return unique sample names
		# generate internal matrix with only unique sample names
		# in the case that an individual occurs more than once (in train & predict or as parent of >1 bipop), it'll get merged with its replicates by alter.pop
		options(warn=-1)		
		tmp <- alter.pop(train) #merges the pops
		options(warn=0)
		TP <- list(theta=tmp@theta[train.m,],geno=tmp@geno[train.m,])
		gid <- colnames(TP$theta)
		n <- length(gid)
								
		# separate values from TP to train vs predict when min.train!=NULL
		if(!is.null(min.train)){
			TP.to.pred <- matrix(FALSE, m, n) # index of those to move to pred
			rownames(TP.to.pred) <- train.m; colnames(TP.to.pred) <- gid
			i.pred <- lapply(train, function(x) intersect(train.m, x@info$marker[ which(x@info$called==FALSE)])) # those not called
			for(f in 1:length(train)){
				TP.to.pred[i.pred[[f]], which( gid%in% colnames(train[[f]]@theta))] <- TRUE
			}
		} else {
			TP.to.pred <- NULL
		}

	if(is.null(predict)){
		p.theta <- NULL
		nn <-0
	} else {
		P <- setdiff(colnames(predict@theta), gid)
		nn <- length(P)
		p.theta <- matrix(NA,m,nn)
		rownames(p.theta) <- train.m
		colnames(p.theta) <- P
		
		if(!nn >0){
			stop("all samples in the prediction set have already been assigned genotypes in the training data")
		}
		
		mark <- intersect(train.m, rownames(predict@theta))
		p.theta[mark, ] <- predict@theta[mark, colnames(p.theta)]	
	}
	
	if((n.core>1) & requireNamespace("parallel", quietly=TRUE)){
		it <- split(1:m, factor(cut(1:m, n.core, labels=FALSE)))
		ANS <- unlist(parallel::mclapply(it, fz, TP=TP, p.theta=p.theta,TP.to.pred=TP.to.pred, mc.cores=n.core),recursive=FALSE)
	}else{
		ANS <- fz(1:m,TP=TP,p.theta=p.theta, TP.to.pred=TP.to.pred)			
	}
	
	acc <- sapply(ANS, function(x){x$acc})
	names(acc) <- train.m

	if(is.null(min.train))	{
		names(ANS) <- train.m	
		geno <- t(sapply(ANS,function(x){x$data$geno}))
		colnames(geno) <- c(colnames(TP$theta),colnames(p.theta))
	}else{
		geno <- matrix(rep(NA, m*(n+nn)), m, n+nn )
		colnames(geno) <- c(colnames(TP$theta),colnames(p.theta))
		rownames(geno) <- train.m
		
		for (a in 1:length(ANS)){
			geno[train.m[a], rownames(ANS[[a]]$data)] <- ANS[[a]]$data$geno
		}
				
	}	
	

	
	# expand to dimensions of all provided data
	if(!is.null(predict)){
		present2 <- unique(c(present, rownames(predict@theta)))
	}else{
		present2 <- present	
	}
	
	theta2 <- matrix(NA, length(present2), n+nn); rownames(theta2) <- present2; colnames(theta2) <- c(gid, colnames(p.theta))
	geno2 <- theta2
		
	for(i in 1:n.fam){
			mark <- intersect(present2, rownames(train[[i]]@theta))
			theta2[mark, colnames(train[[i]]@theta)] <- train[[i]]@theta[mark, ]
		}
	if(!is.null(predict)){
		mark <- intersect(present2, rownames(predict@theta))
		theta2[mark, colnames(predict@theta) ] <- predict@theta[mark, ]
	}	
		
	geno2[rownames(geno) , colnames(geno)] <- geno[rownames(geno), colnames(geno)]		



	pops <- as.list(rep(NA, length(train)))
	pops <-lapply(train, function(X) c( which(colnames(geno2) %in% c(X@parent1, X@parent2)), which(colnames(geno2) %in% setdiff(colnames(X@theta), c(X@parent1, X@parent2)))))
	
	names(pops) <- lapply(train, function(x) 
		if(!is.na(x@parent1)){
			if(x@parent1 != x@parent2) {
				paste0(substr(x@parent1, 1, 1), "x", substr(x@parent2, 1, 1))
			}else{
				paste0(substr(x@parent1, 1, 1), "_selfed")
			}
		}else{
			paste0("diverse")	
		})
		
	if(!is.null(predict)){
		N <- c(n.fam+1)
		pops[[N]] <- which(colnames(geno2) %in% colnames(predict@theta))
		names(pops)[N] <- list("prediction")
	}

	concordance <- rep(NA, length(present2)); names(concordance) <- present2
	concordance[names(acc)] <- acc
	
	return(new("anypop", concordance = concordance, theta=theta2, geno=geno2, pops=pops))
}









