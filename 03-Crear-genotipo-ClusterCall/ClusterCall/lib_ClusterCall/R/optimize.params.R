
optimize.params <- function(bipops, anypop=NULL, max.missing=0.05, max.range=0.25, hom.theta=0.15, min.posterior=0.7, min.sep=0.05, min.r=0.2, mono.thresh=0.95, chi2p.step=1,min.train=NULL, n.core=1){
	pops <- names(bipops)
	if (is.null(min.train)) {
		tests <- expand.grid(max.missing=max.missing, max.range=max.range, hom.theta=hom.theta, min.posterior=min.posterior, min.sep=min.sep,min.r=min.r, mono.thresh=mono.thresh, chi2p.step=chi2p.step,KEEP.OUT.ATTRS=F)
		tests$min.train <- rep(NA,nrow(tests))	
	} else {
		tests <- expand.grid(max.missing=max.missing, max.range=max.range, hom.theta=hom.theta, min.posterior=min.posterior, min.sep=min.sep,min.r=min.r, mono.thresh=mono.thresh,chi2p.step=chi2p.step, min.train=min.train,KEEP.OUT.ATTRS=F)		
	}
	
	t <- dim(tests)[1]
		
	CC <- as.list(rep(NA, t)) 
	bipop.results <- vector("list", t)
	
	for (i in 1:t){
		bipop.results[[i]] <-lapply(bipops, function(x) CC.bipop(new("pop",theta=x@theta,r=x@r), parent1=x@parent1, parent2=x@parent2, max.missing=tests$max.missing[i], max.range=tests$max.range[i], hom.theta=tests$hom.theta[i], min.posterior=tests$min.posterior[i], min.sep=tests$min.sep[i], min.r=tests$min.r[i],mono.thresh=tests$mono.thresh[i], chi2p.step=tests$chi2p.step[i],n.core=n.core))  			
		CC[[i]] <- CC.anypop(train=bipop.results[[i]], n.core=n.core, max.missing=tests$max.missing[i], max.range=tests$max.range[i], min.posterior=tests$min.posterior[i], hom.theta=tests$hom.theta[i],min.sep=tests$min.sep[i], min.train=if(is.na(tests$min.train[i])){NULL}else{tests$min.train[i]})
	}
return(list(tests=tests, anypop.results=CC, bipop.results=bipop.results))	
	
}