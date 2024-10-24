
\name{optimize.params}
\alias{optimize.params}
\title{Determine optimal parameters}
\usage{
optimize.params(bipops, anypop = NULL, max.missing = 0.05,
  max.range = 0.25, hom.theta = 0.15, min.posterior = 0.7,
  min.sep = 0.05, min.r = 0.2, mono.thresh = 0.95, chi2p.step = 1,
  min.train=NULL, n.core = 1)
}
\arguments{
\item{bipops}{List of objects of class \code{\linkS4class{bipop}}}

\item{anypop}{Object of class \code{\linkS4class{pop}} for which marker genotypes are to be predicted}

\item{max.missing}{Maximum proportion of samples assigned NA genotype}

\item{max.range}{Maximum allowed spread of a cluster}

\item{hom.theta}{Maximum distance between theta = 0 and the center of a cluster with genotype = 0, or between theta = 1 and a cluster with genotype = 4}

\item{min.posterior}{Minimum posterior membership probability of the sample to its assigned cluster}

\item{min.sep}{Minimum distance between medians of neighboring clusters}

\item{min.r}{Minimum threshold for r values}

\item{mono.thresh}{Minimum proportion of samples to belong to a single cluster for the marker to be assigned a single genotype}

\item{chi2p.step}{Incremental change in -log10(p) goodness of fit required to reduce the number of clusters}

\item{min.train}{Minimum number of bipops with genotype calls at a given marker for the marker to be predicted with \code{\link{CC.anypop}}. By default, the marker must be called in all bipops.}

\item{n.core}{Number of cores to use for multicore processing}
}
\value{
	A list of length three with:
	\itemize{ 
	\item{\code{tests}}: {A data frame of tests by parameters indicating the parameter values for each run}
	\item{\code{anypop.results}}: {A list of \code{\linkS4class{anypop}} results corresponding to the rows of \code{tests}}
	\item{\code{bipop.results}}: {A list of lists of \code{\linkS4class{bipop}} results (one per family in the nested list) corresponding to the rows of \code{tests}}
	}
}
\description{
Runs all possible combinations of supplied parameters with \code{\link{CC.bipop}} on each family supplied as an object of class \code{\linkS4class{bipop}} and then \code{\link{CC.anypop}} is used to determine concordance of the resulting calls.  
}
\details{
This function provides an organized approach to optimizing parameters for \code{\link{CC.bipop}} and \code{\link{CC.anypop}}. All possible combinations of indicated parameters are analyzed. 

To test multiple levels of any parameter (excluding \code{n.core}), pass a vector to the parameter. For instance, setting \code{hom.theta=c(0.1, 0.14)} will produce results for \code{hom.theta=0.1} and \code{hom.theta=0.14}. 
	
}
\seealso{
\code{\link{CC.bipop}}, \code{\link{CC.anypop}}

\code{\link{CC.plot}} to visually compare tests
}

\references{
Schmitz Carley et al. (2017) Automated tetraploid genotype calling by hierarchical clustering. Theoretical and Applied Genetics 130:717–726 DOI: 10.1007/s00122-016-2845-5. 
}

\examples{
##see Vignette at http://potatobreeding.cals.wisc.edu/software
}