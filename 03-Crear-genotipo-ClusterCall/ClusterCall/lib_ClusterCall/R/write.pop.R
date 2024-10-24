
write.pop <- function(result, file, thresh=NA){
	
	if (sum(inherits(result, "bipop"), inherits(result, "anypop")!=1)){
		stop("result must be a bipop or anypop object")
	}
	
	if (inherits(result, "bipop")){
		
		write.csv(result@geno, file=file, row.names=TRUE)
		
	}else{
		if(!is.na(thresh)){
			write.csv(result@geno[names(result@concordance)[which(result@concordance >= thresh)], ], file=file, row.names=TRUE)
		}else{
			write.csv(result@geno, file=file, row.names=TRUE)
		}
	}	
	
}