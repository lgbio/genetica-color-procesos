
.find.matching.sets <- function(theta, thresh, method){
	
	ans <- hclust(dist(t(theta)), method=method)
	assign <- cutree(ans, h=thresh)
	ind <- unique(assign)
	
	pick <- unname(which(table(assign)>1))
	
	sets <- lapply(1:max(ind), function(x) {names(assign)[which(assign==x)]})
	names(sets) <- ind
	
	matches <- sets[pick]
		
	return(matches)
}	

.consensus <- function(theta, indivs){
	cons<-apply(theta[, indivs], 1, median, na.rm=T)
	return(cons)
	}
	
	
	fx <- function(g) {
		if (length(unique(g[!is.na(g)]))==1){
			return(unique(g[!is.na(g)]))
		}else{
			return(NA)	
		}
	}
	
.consensusgeno <- function(geno, indivs){
	Mode <- function(x){
		ux <- unique(x[!is.na(x)])
		if(max(tabulate(match(x[!is.na(x)], ux))) > 1){
			ux[which.max(tabulate(match(x[!is.na(x)], ux)))]
		}else{
			if(length(unique(x[!is.na(x)]))==1){
				unique(x[!is.na(x)])
			}else{
				return(NA)
			}
			
		}
	}
	cons <- apply(geno[, indivs], 1, function(x) Mode(x))
	return(cons)
}

.seg.p <- function(y) {
	if (!all(!is.na(y[1:2]))) {return(NA)}
	options(warn=-1)
	tab <- table(factor(y[-(1:2)],levels=0:4))
	n <- sum(tab)
	p1 <- min(y[1:2])
	p2 <- max(y[1:2])
	if (p1==0) {
		if (sum(tab["3"],tab["4"])>0) {return(NA)}  #3 and 4 not possible when p1==0 	
	}
	if (p2==4) {
		if (sum(tab["0"],tab["1"])>0) {return(NA)}  #0 and 1 not possible when p2==4 	
	}
	if (p1==0) {
		if ((p2==0)|(p2==4)) {return(NA)}
		xp <- switch(p2,
			data.frame(x=c(tab["0"],tab["1"]),p=c(1,1)), #p2==1, DR can produce 2
			data.frame(x=c(tab["0"],tab["1"],tab["2"]),p=c(1,4,1)), #p2==2
			data.frame(x=c(tab["1"],tab["2"]),p=c(1,1)) #p2==3, DR can produce 0
		)
	}
	if (p1==1) {
		xp <- switch(p2,
			data.frame(x=c(tab["0"],tab["1"],tab["2"]),p=c(1,2,1)), #p2==1, DR can produce 3 or 4			
			data.frame(x=c(tab["0"],tab["1"],tab["2"],tab["3"]),p=c(1,5,5,1)), #p2==2, DR can produce 4 
			data.frame(x=c(tab["1"],tab["2"],tab["3"]),p=c(1,2,1)), #p2==3, DR can produce 0 or 4
			data.frame(x=c(tab["2"],tab["3"]),p=c(1,1))  #p2==4, DR can produce 4
		)
	}		
	if (p1==2) {
		xp <- switch(p2-1,	
			data.frame(x=as.vector(tab),p=c(1,8,18,8,1)), #p2==2
			data.frame(x=c(tab["1"],tab["2"],tab["3"],tab["4"]),p=c(1,5,5,1)), #p2==3, DR can produce 0
			data.frame(x=c(tab["2"],tab["3"],tab["4"]),p=c(1,4,1))  #p2==4
		)
	}		
	if (p1==3) { 
		xp <- switch(p2-2,	
			data.frame(x=c(tab["2"],tab["3"],tab["4"]),p=c(1,2,1)), #p2==3, DR can produce 0 or 1
			data.frame(x=c(tab["3"],tab["4"]),p=c(1,1))  #p2==4, DR can produce 2
		)
	}		
	if (p1==4) {return(NA)}
	return(chisq.test(xp$x,p=xp$p/sum(xp$p))$p.value)
	options(warn=0)	
}

.bayes.check <- function(theta,cluster,models,n.clust) {
	theta <- ifelse(theta < 1e-4,1e-4,theta)
	theta <- ifelse(theta > 1 - 1e-4,1-1e-4,theta)

	ans <- apply(array(theta),1,function(x){
		sapply(models,function(y){
			return(switch(y$type,
					dt((x-y$estimate[1])/y$estimate[2],df=y$estimate[3])/y$estimate[2],
					dlnorm(x,meanlog=y$estimate[1],sdlog=y$estimate[2]),
					dlnorm(1-x,meanlog=y$estimate[1],sdlog=y$estimate[2]),
					dnorm(x,mean=y$estimate[1],sd=y$estimate[2]),
					0)
			)
		})			
	})

	posterior <- t(ans)/apply(ans,2,sum)
	max.posterior <- apply(posterior,1,max,na.rm=T)
	return(ifelse(apply(posterior,1,which.max)==cluster,max.posterior,0))
}

.fit.dist <- function(x,center,hom.theta) {
	n <- length(x)
	x <- ifelse(x < 1e-4,1e-4,x)
	x <- ifelse(x > 1 - 1e-4,1 - 1e-4,x)
	if ((n==1)|(sd(x) < 1e-4)) {
		return(list(type=5))
	}
	options(warn=-1)
	ans <- vector("list",length=4)
	ans <- try(fitdistr(x,densfun="t"),silent=TRUE)
	if (class(ans)!="try-error") {
		return(list(type=1,estimate=ans$estimate))
	} 
	if (center <= hom.theta) {
		ans <- try(fitdistr(x,densfun="lognormal"),silent=TRUE)
		if (class(ans)!="try-error") {return(list(type=2,estimate=ans$estimate))}
	}
	if (center >= 1 - hom.theta) {
		ans <- try(fitdistr(1-x,densfun="lognormal"),silent=TRUE)
		if (class(ans)!="try-error") {return(list(type=3,estimate=ans$estimate))}
	}
	ans <- try(fitdistr(x,densfun="normal"),silent=TRUE)
	if (class(ans)!="try-error") {return(list(type=4,estimate=ans$estimate))}
	options(warn=0)
	return(NULL)	
}

.hom.check <- function(dosage,centers,hom.theta) {
	if (((min(dosage)==0)&&(centers[dosage==0]>hom.theta))|((max(dosage)==4)&&(centers[dosage==4]< 1 - hom.theta))) {
		return(NULL)
	} 
	if (((min(dosage) > 0)&&(min(centers) <= hom.theta))|((max(dosage) < 4)&&(max(centers) >= 1-hom.theta))) {
		return(NULL)
	}			
	return(dosage)
}


.get.dosage <- function(n.clust,centers,hom.theta) {
		if (n.clust==2) {
			dosage <- if (centers[1] <= hom.theta) {.hom.check(0:1,centers,hom.theta)} else {.hom.check(3:4,centers,hom.theta)}
		}
		if (n.clust==3) {
			if (centers[1] <= hom.theta) {
				dosage <- .hom.check(0:2,centers,hom.theta)
			} else {
				if (centers[n.clust] >= 1 - hom.theta) {
					dosage <- .hom.check(2:4,centers,hom.theta)
				} else {
					dosage <- .hom.check(1:3,centers,hom.theta)
				}
			}
		}
		if (n.clust==4) {
			dosage <- if (centers[1] <= hom.theta) {.hom.check(0:3,centers,hom.theta)} else {.hom.check(1:4,centers,hom.theta)}
		}
		if (n.clust==5) {
			dosage <- .hom.check(0:4,centers,hom.theta)
		}
		if (!is.null(dosage)) {
			names(dosage) <- names(centers)
		}
		return(dosage)
}

.biparental.call <- function(theta,r=NULL,min.r=0.2,min.posterior=0.7,max.missing=0.05,min.sep=0.05,max.range=0.25,mono.thresh=0.95,hom.theta=0.15,chi2p.step=1) {
	#first two elements of theta (and r) must be the parents
	cluster.method <- "complete"  #almost no difference in concordance using complete vs. average method, no clear winner
	if (!is.null(r)) {
		ix <- which(!is.na(r)&(r >= min.r)&!is.na(theta))
	} else {
		ix <- which(!is.na(theta))
	}
	n <- length(theta)
	geno <- integer(n)*NA
	called <- FALSE
	n2 <- length(ix)
	best.chi2p <- NA

	if (is.element(1,ix)&is.element(2,ix)&(n2/n >= 1-max.missing)) {
		theta2 <- theta[ix]		
		goo <- hclust(dist(theta2),method=cluster.method)
		valid.clust <- integer(0)
		bayes.ans <- vector("list",5)
		for (n.clust in 2:5) {
			goo2 <- cutree(goo,n.clust)
			tab <- table(goo2)
			spread <- tapply(theta2,factor(goo2),function(x){diff(range(x))})
			centers <- tapply(theta2,factor(goo2),median)		
			if ((min(diff(sort(centers))) >= min.sep)&(max(spread) <= max.range)) {
				models <- vector("list",length=n.clust)
				proceed <- TRUE
				for (k in 1:n.clust) {
					x <- theta2[goo2==k]
					tmp <- .fit.dist(x,centers[k],hom.theta)
					if (!is.null(tmp)){models[[k]] <- tmp} else {proceed <- FALSE}
				}
				if (proceed) {
					type5 <- which(sapply(models,function(x){return(x$type==5)}))
					NA.bayes <- which(.bayes.check(theta2,goo2,models,n.clust) < min.posterior & !is.element(goo2,type5))
					if ((length(NA.bayes)/n <= max.missing)&!is.element(1,NA.bayes)&!is.element(2,NA.bayes)) {
						valid.clust <- c(valid.clust,n.clust)
						bayes.ans[[n.clust]] <- NA.bayes
					}			
				}
			}
		}
		if (length(valid.clust)==0) {
			if (diff(range(theta2)) <= max.range) {
				max.clust <- 1
			} else {
				max.clust <- 0
			}
		} else {
			max.clust <- max(valid.clust)
		}
		if (max.clust > 0) {
			goo2 <- cutree(goo,max.clust)
			centers <- sort(tapply(theta2[names(goo2)],factor(goo2),median))
			prop <- table(goo2)/n2
		
			if (max(prop) > mono.thresh) {
				in.big.clust <- which(goo2==as.integer(names(prop)[which.max(prop)]))
				big.clust.center <- centers[names(prop)[which.max(prop)]]
				if (diff(range(theta2)) <=  max.range){	
					if (big.clust.center <= hom.theta) {
						geno[ix] <- 0
						called <- TRUE
					}
					if (1-big.clust.center <= hom.theta) {
						geno[ix] <- 4
						called <- TRUE
					} 						
				}else{
					# could be 0 x 4 pattern
					p1 <- theta2[1]
					p2 <- theta2[2]
					if (p1 <= hom.theta) {
						if (1-p2 <= hom.theta) {
							geno[ix] <- 2
							geno[1] <- 0
							geno[2] <- 4
							called <- TRUE
						}
					}
					if (p2 <= hom.theta) {
						if (1-p1 <= hom.theta) {
							geno[ix] <- 2
							geno[1] <- 4
							geno[2] <- 0
							called <- TRUE
						}
					}							
					if(!called){
						# assign theta in big.clust to homozygotes, leave the rest unassigned
						b <- which(names(theta) %in% names(goo2)[in.big.clust])					
						if (big.clust.center <= hom.theta){
					 		geno[b] <- 0
				 			called<-TRUE
				 		}else{
					 		if(1-big.clust.center<=hom.theta){
					 			geno[b] <-4
		 						called <- TRUE
							}
						}							
					}
				}
			} else {  
				for (n.clust in sort(valid.clust,decreasing=TRUE)) {  #go from highest to lowest
					goo2 <- cutree(goo,n.clust)
					centers <- sort(tapply(theta2[names(goo2)],factor(goo2),median))
					dosage <- .get.dosage(n.clust,centers,hom.theta)
					if (!is.null(dosage)) {
						x <- as.integer(dosage[as.character(goo2)])
						chi2p <- .seg.p(x)
						if (!is.na(chi2p)&&(is.na(best.chi2p)|(log10(chi2p) > log10(best.chi2p)+chi2p.step))) {# chi2p.step is penalty for fewer clusters
							best.chi2p <- chi2p
							if (length(bayes.ans[[n.clust]]) > 0) {x[bayes.ans[[n.clust]]] <- NA}
							geno[ix] <- x
							called <- TRUE
						}
					}
				}
			}
		}
	}	
	return(list(geno=geno,chi2p=best.chi2p,called=called))
}


.theta2geno <- function(train,predict=NULL,max.missing=0.05,min.sep=0.05,max.range=0.25,min.posterior=0.7,hom.theta=0.15,impute=TRUE) {
	#train must have rownames, as well as predict
	#train is data frame with two columns: theta, geno
	#predict is vector of theta values to predict
	#apply max.missing to TP and pred.set separately

	if (is.null(predict)) {
		data <- train
	} else {
		data <- data.frame(theta=c(train$theta,predict),geno=c(train$geno,1*rep(NA,length(predict))))
		rownames(data) <- c(rownames(train), names(predict))
	}
	
	not.NA <- which(!is.na(data$theta))
	data2 <- data[not.NA,]
	n <- length(not.NA)
	set1 <- which(rownames(data2) %in% rownames(train))
	if (!is.null(predict)){set2 <- which(rownames(data2) %in% names(predict))}
	t.set <- which(!is.na(data2$geno))
	
	# one cluster
	if (diff(range(data2$theta)) <= max.range) {
		tab <- table(data2$geno)
		best.geno <- rep(as.integer(names(tab)[which.max(tab)]),n)
		best.acc <- length(which(best.geno[t.set]==data2$geno[t.set]))/length(t.set)
		best.k <- 1
		best.spread <- diff(range(data2$theta))
	} else {
		best.acc <- NA
		best.geno <- rep(NA,n)
		best.spread <- NA
	}
			
	for (cluster.method in c("complete","ward.D2","average")) {  #Including multiple methods really does produce more curated markers!
		goo <- hclust(dist(data2$theta),method=cluster.method)
		n.clust <- 5
		proceed <- TRUE
		while ((n.clust >= 2)&(proceed)) {
			geno <- vector("integer",n)*NA		
			goo2 <- cutree(goo,n.clust)	
			spread <- tapply(data2$theta,factor(goo2),function(x){diff(range(x))})
			centers <- tapply(data2$theta,factor(goo2),median)	
			if ((min(diff(sort(centers))) >= min.sep)&(max(spread) <= max.range)) {
				models <- vector("list",length=n.clust)
				proceed2 <- TRUE
				for (k in 1:n.clust) {
					x <- data2$theta[goo2==k]
					tmp <- .fit.dist(x,centers[k],hom.theta)
					if (!is.null(tmp)){models[[k]] <- tmp} else {proceed2 <- FALSE}
				}
				if (proceed2) {
					type5 <- which(sapply(models,function(x){return(x$type==5)}))
					NA.bayes <- .bayes.check(data2$theta,goo2,models,n.clust) < min.posterior & !is.element(goo2,type5)
					# if does not meet compact criterion, remove from training set
					t.set <- intersect(t.set,which(NA.bayes==FALSE))
					tab <- table(goo2[t.set], data2$geno[t.set])
					cluster.dosage <- as.integer(colnames(tab)[apply(tab,1,which.max)]) 
					if (length(which(duplicated(cluster.dosage)))==0) {
						for (j in 1:length(cluster.dosage)) {geno[goo2==j] <- cluster.dosage[j]}				
						missing.clust <- setdiff(1:n.clust,as.integer(rownames(tab)))
						if (length(missing.clust) > 0) {
							#at least one of the clusters has no training set samples
							#but meets min.sep, max.range criteria, do not allow any smaller n.clust to avoid
							#spurious lumping of prediction set clusters
							proceed <- FALSE
							if (length(cluster.dosage)==4) {  #make guess about homozygous cluster
								if ((centers[missing.clust] <= hom.theta)&!(0 %in% cluster.dosage)) {
									geno[goo2==missing.clust] <- 0 
								}
								if ((centers[missing.clust] >= 1 - hom.theta)&!(4 %in% cluster.dosage)) {
									geno[goo2==missing.clust] <- 4 
								}
							}
						}
						geno2 <- geno
						if (length(which(NA.bayes)) > 0) {geno[NA.bayes] <- NA}
						set1.NA <- length(which(is.na(geno[set1])))/length(set1)
						if (!is.null(predict)){
							set2.NA <- length(which(is.na(geno[set2])))/length(set2)
						} else {
							set2.NA <- 0
						}
						if (max(set1.NA,set2.NA) <= max.missing) {
							acc <- length(which(geno[t.set]==data2$geno[t.set]))/length(t.set)
							if (is.na(best.acc)|(acc > best.acc)|((acc == best.acc)&(max(spread) < best.spread))) {
								best.acc <- acc
								best.k <- n.clust
								if (impute) {best.geno <- geno2} else {best.geno <- geno}
								best.spread <- max(spread)
							}
						}
					}
				}
			}
			n.clust <- n.clust - 1	
		}
	}
	data$geno[not.NA] <- best.geno
	return(list(acc=best.acc,data=data))
}








.runfitT_write <- function(theta=theta, file=file, params=params){
  # theta = matrix of theta values
  # rownames are marker names, colnames are sample names
  # name = project name, names output files
  
	if(requireNamespace("reshape", quietly=TRUE)){
			fitform <- reshape::melt(theta) 
			colnames(fitform) <- c("MarkerName", "SampleName", "ratio")
	}

	gid <- colnames(theta)
	mkrs <- rownames(theta)
	geno <- matrix(rep(NA, length(gid)*length(mkrs)),
	               nrow=length(mkrs), ncol=length(gid))
		colnames(geno) <- gid
		rownames(geno) <- mkrs

	

	for (m in 1:length(mkrs)){
		dat <- fitform[which(fitform$MarkerName==mkrs[m]), ]
		dat$MarkerName <- as.character(dat$MarkerName)
		dat$SampleName <- as.character(dat$SampleName)
		ans <- try(
			if(requireNamespace("fitTetra", quietly=TRUE)){
				fitTetra::saveMarkerModels(data=dat, 
					sd.threshold=params[["sd.threshold"]], 
					p.threshold=params[["p.threshold"]], 
					call.threshold=params[["call.threshold"]], 
					peak.threshold=params[["peak.threshold"]], 
					try.HW=params[["try.HW"]], 
					dip.filter=params[["dip.filter"]], 
					plot="none",     
					scorefile=paste0(mkrs[m],file, "_S.dat"), 
					modelfile=paste0(mkrs[m], file, "_M.dat"))} , silent=T
			
			)			

		}
	}



.runfitT_read <- function(theta, file){
	
  gid <- colnames(theta)
  geno <- matrix(NA, ncol=dim(theta)[2], nrow=dim(theta)[1])
  	for (i in 1:dim(theta)[1]){
		options(warn=1)
	    suppressWarnings(ans <- try(read.table ( paste0( rownames(theta)[i], 
	                                    file, "_S.dat"), header=T, sep="\t"),
	                                    silent=T))
	      if(class(ans)=="data.frame"){
	      	# warnings suppressed for markers that were not called, so no file to read
			  	suppressWarnings(score <- read.table(paste0(rownames(theta)[i],
			  	                           file, "_S.dat"), header=T, sep="\t"))
			  	if(requireNamespace("reshape", quietly=TRUE)){                        
			  		score.mat <- reshape::cast(markername~sample, data=ans, value="geno")	
			  
			  	}	
			  	score.mat2 <- score.mat[ ,-1]
			  	rownames(score.mat2) <- score.mat$markername
					
		    	geno[i, ] <- as.matrix(score.mat2[, gid])
	      }	
	    rownames(geno) <- rownames(theta)
	    colnames(geno) <- colnames(theta)
	
	  }
	write.csv(geno, file=paste0(file, "_fitTonTheta.csv"), row.names=T)
}


.make.bipop <- function(scores, theta){
	# scores = genotype calls from fitTetra
	# theta = a matrix of theta values
  	# columns of theta include all individuals in scores	
  
		gid <- colnames(scores)
		mark <- rownames(scores)
		scores <- as.matrix(scores)
		
		# poly if more than one dose is assigned at the marker
		poly <- apply(scores, 1, function(x) {length(unique(x[!is.na(x)]))>1}) 
		# if poly==TRUE, chi2p is set to 1, else NA (for CC.anypop to use it)
		chi2p <- ifelse(poly, 1, NA) 
		called <- ifelse(poly, TRUE, FALSE)
		mkr.data <- data.frame(marker=as.character(rownames(scores)), 
		                       called=called, negLog10p=chi2p)
		# fitTetra only returns called markers
		mkr.data$called[which(mkr.data$marker %in% rownames(scores))]<- TRUE 
		
		bipop <- new("bipop", theta=theta, geno=as.matrix(scores), 
		             info=mkr.data, parent1=as.character(NA), parent2=as.character(NA))
	
	return(bipop)

}

