
CC.plot <- function(anypop, names=NULL, xlim=c(1,0),...){

	if(length(anypop)==1){
		
		if(length(anypop@concordance) < 1){
			stop("Concordance not found for this object. Use the CC.anypop function to calculate concordance.")
		}	
		
		acc <- anypop@concordance
		x <- sort(acc)
		ac.ct <- apply(array(x), 1, function(x, acc){length(which(acc>=x))}, acc)
		plot(x, ac.ct, type='l', ylab="Markers", xlab="Concordance", xlim=xlim, ...)
	}else{
		line <- c(1:length(anypop))
		color <- c(rep(c("black", "red", "blue"), length(anypop)))
		
		if(length(anypop[[1]]@concordance) <1){
			stop("Concordance not found for the first object. Use the CC.anypop function to calculate concordance.")
		}
		
		acc <- anypop[[1]]@concordance
		x <- sort(acc)
		ac.ct <- apply(array(x), 1, function(x, acc){length(which(acc>=x))}, acc)
		plot(x, ac.ct, type='l', ylab="Markers", xlab="Concordance", xlim=xlim, lty=line[1], ...)
		
		for (i in 2:length(anypop)){
			
			if(length(anypop[[i]]@concordance) <1){
				stop(paste0("Concordance not found for object ", i, " in the anypop list. Use the CC.anypop function to calculate concordance."))
			}
			
			acc <- anypop[[i]]@concordance
			x <- sort(acc)
			ac.ct <- apply(array(x), 1, function(x, acc){length(which(acc>=x))}, acc)
			points(x, ac.ct, type='l', ylab="Markers", xlab="Concordance", xlim=xlim, lty=line[i], col=color[i])
		
		}
		
		if(!is.null(names)){
			legend("bottomright", legend=names, lty=line, col=color)
		}
		
	}
	
}