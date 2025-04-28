
read.pop <- function(theta.file, r.file=NULL, error.checking=FALSE, method="complete", thresh=5, max.missing.theta=0.2, ...){
	
    print (">>> "); print (theta.file); print (r.file)
	theta <- as.matrix(read.csv(file=theta.file, header=T, as.is=T, check.names=F, row.names=1))
	thresh=thresh # dendrogram height threshold for calling individuals identical
	
    #compare_csv_headers_rows (theta.file, r.file)

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
	
		if (!all(colnames(theta) == colnames(r)) || !all (rownames(theta) == rownames(r))){
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

compare_csv_headers_rows <- function(tFile, rFile) {
  # Read the input files
  theta <- as.matrix(read.csv(file=tFile, header=TRUE, as.is=TRUE, check.names=FALSE, row.names=1))
  r <- as.matrix(read.csv(file=rFile, header=TRUE, as.is=TRUE, check.names=FALSE, row.names=1))

  # Compare column names
  col_comparison <- data.frame(
    Theta_Columns = colnames(theta),
    R_Columns = c(colnames(r), rep(NA, max(0, length(colnames(theta)) - length(colnames(r)))))
  )
  write.csv(col_comparison, "column_comparison.csv", row.names=FALSE)

  # Compare row names
  row_comparison <- data.frame(
    Theta_Rows = rownames(theta),
    R_Rows = c(rownames(r), rep(NA, max(0, length(rownames(theta)) - length(rownames(r)))))
  )
  write.csv(row_comparison, "row_comparison.csv", row.names=FALSE)

  print("Comparison files created: column_comparison.csv and row_comparison.csv")
}

