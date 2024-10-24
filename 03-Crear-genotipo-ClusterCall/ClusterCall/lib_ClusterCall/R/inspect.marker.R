
inspect.marker <- function(called, marker, ...){
		
		if(length(called)==1){
			if(sum(inherits(called, "bipop"), inherits(called, "anypop")) !=1)
				stop("called must be an object of class bipop or anypop")
					
		}else{
			
			for (f in 1:length(called )){
				if(sum (inherits(called[[f]], "bipop"), inherits(called[[f]], "anypop")) !=1)
					stop(paste0("population ", f, " must be of class bipop, or anypop"))
			}			
			called <- alter.pop(called)

		}
			
		if(is.character(marker)){
			if(length(intersect(marker, rownames(called@geno)))!=1){
				stop("the indicated marker was not called")
			}	
		}	
		
		if(!is.character(marker)){
			marker <- rownames(called@geno)[marker]
		}
		
		if(inherits(called, "bipop")){
			
			plot(called@theta[marker,], pch=as.character(called@geno[marker,]), xlab="Sample", ylab="Theta", ylim=c(0,1), ...)
			na <- which(is.na(called@geno[marker,]))
			points(na, called@theta[marker,na], pch=20, col="red")
			
			if(!is.na(called@parent1)){

				if(called@parent1 != called@parent2){
					points(called@theta[marker, called@parent1], pch="O", col="red", cex=1.5)
					points(2, called@theta[marker,  called@parent2], pch="O", col="blue", cex=1.5)
				}else{
					points(called@theta[marker, called@parent1], pch="O", col="purple", cex=1.5)
				}
			}
		}else{
			# for anypop, reconstruct the matrix with duplicates, report accuracy and add vertical breaks between populations
			n <- length(unlist(called@pops))
			gid <- colnames(called@geno)[unlist(called@pops)]
			
			theta <- rep(NA, n); names(theta) <- gid
			theta <- called@theta[marker,gid]
			
			geno <- rep(NA, n); names(geno) <- gid
			geno <- called@geno[marker, gid]

			plot(theta, pch=as.character(geno), xlab="Sample", ylab="Theta", ylim=c(0,1), ...)
			na <- which(is.na(geno))
			points(na, theta[na], pch=20, col="red")
				
			title(sub=paste0("concordance: ", round(called@concordance[marker], 2)), col.sub="blue")	 
						
			line <- unlist(lapply(called@pops, function(x) length(x)))[1:length(called@pops)]
			line <- line[1:length(line)-1]
			
			if(length(line) >1){
				for (l in 2:length(line)){
					line[l] <- line[l-1] +line[l]
				}
			}
			abline(v=c(line+.5), lty=2)
			
			mtext(names(called@pops), side=3, at=lapply(called@pops, function(x) median(x)))
			
			
			 # circle parents
			 for(p in 1:length(called@pops)){
				 if(names(called@pops[p]) != "prediction" ){
					 if(grepl("selfed", names(called@pops[p]))){
						 called@pops[p][[1]][1] <- "purple" 
					 }else{
						 called@pops[p][[1]][1:2] <- c("red", "blue")
					 }	
				 }	
			 }
			 purp <- which(unlist(called@pops)=="purple")
			 red <- which(unlist(called@pops)=="red")
			 blue <- which(unlist(called@pops)=="blue")
			 points(purp, theta[purp], pch="O", col="purple", cex=1.5)
			 points(red, theta[red], pch="O", col="red", cex=1.5)
			 points(blue, theta[blue], pch="O", col="blue", cex=1.5)

			
			}
			title(marker)
	
}