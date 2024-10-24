
\name{usefitT}
\alias{usefitT}
\title{Interface fitTetra and ClusterCall}
\description{
Provides an interface to the fitTetra package	
}

\usage{
	usefitT(theta, params=NULL, n.core=1)
}
\arguments{
\item{theta}{ m x n matrix of theta values with marker names as rownames and sample names as column names}
\item{params}{optional list of parameters for the fitTetra function \code{saveMarkerModels}. See Details.}
\item{n.core}{number of cores to use for parallel processing}
}

\details{
This runs the \code{saveMarkerModels} function in R package fitTetra for each marker in \code{theta} and returns the result as a genotype matrix within a variable of class \code{\linkS4class{bipop}}, regardless of whether the dataset is actually a biparental population. This is done so that the fitTetra results can be combined with true bipop training sets to calculate concordance using \code{\link{CC.anypop}}. This workflow is illustrated in the vignette available at the URL shown below. The optional \code{params} argument needs to be a named list, e.g., \code{params = list(sd.threshold=0.05, p.threshold=0.90)}
}
\value{
An object of class \code{\linkS4class{bipop}}
}
\references{
Schmitz Carley et al. (2017) Automated tetraploid genotype calling by hierarchical clustering. Theoretical and Applied Genetics 130:717–726 DOI: 10.1007/s00122-016-2845-5. 
	
Voorrips et al. (2011)  Genotype calling in tetraploid species from bi-allelic marker data using mixture models. BMC Bioinformatics 12:172.	
}

\examples{
##see Vignette at http://potatobreeding.cals.wisc.edu/software
}