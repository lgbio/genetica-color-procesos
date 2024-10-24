
\name{CC.bipop}
\alias{CC.bipop}
\title{ClusterCall for a biparental population}
\description{
Make genotype calls for a biparental family
}
\usage{
CC.bipop(pop, parent1, parent2, max.missing = 0.05, min.sep = 0.05, max.range = 0.25,
  hom.theta = 0.15, min.posterior = 0.7, min.r = 0.2, 
  mono.thresh = 0.95, chi2p.step = 1, n.core = 1)
}
\arguments{
\item{pop}{Object of class \code{\linkS4class{pop}} containing a biparental family and two parental samples}

\item{parent1}{Sample name of first parent}

\item{parent2}{Sample name of second parent}

\item{max.missing}{Maximum proportion of samples assigned NA genotype}

\item{max.range}{Maximum allowed spread of a cluster}

\item{hom.theta}{Maximum distance between theta = 0 and the center of a cluster with genotype = 0, or between theta = 1 and a cluster with genotype = 4}

\item{min.posterior}{Minimum posterior membership probability of the sample to its assigned cluster}

\item{min.sep}{Minimum distance between medians of neighboring clusters}

\item{min.r}{Minimum r value to assign a genotype, if r matrix was supplied}

\item{mono.thresh}{Minimum proportion of samples to belong to a single cluster for the marker to be assigned a single genotype (called monomorphic)}

\item{chi2p.step}{Penalty for reducing the number of clusters in terms of the -log10(p-value) for segregation test (see Details)}

\item{n.core}{Number of cores to use for multicore processing}
}

\value{
	An object of class \code{\linkS4class{bipop}}
}


\details{
Hierarchical clustering within a biparental family is used to group samples with similar theta values, and then allele dosages are assigned to clusters based on expected segregation ratios, given the parental genotypes under the proposed clustering solution. 

For each marker, the number of sample clusters is varied from two to five. If a partitioning satisfies the \code{max.range} and \code{min.sep} thresholds and the fraction of samples with cluster membership probability less than \code{min.posterior} is less than \code{max.missing}, the partitioning is passed to the genotype assignment module.

In the genotype assignment module, clusterings with five clusters are assigned genotypes of {0,1,2,3,4} in order of ascending median theta value. With fewer than five clusters, several genotype-to-cluster assignments can be made. \code{CC.bipop} cycles through those possibilities that are consistent with segregation in a biparental population, allowing for double reduction products at a low frequency. 

When several plausible genotype-to-cluster assignments exist, the best solution is reported based on goodness-of-fit (highest p-value) to the expected genotype ratios under the assumption of random bivalent formation (i.e. random chromosome segregation) in meiosis.

Markers exceeding thresholds for all partitionings into two to five clusters are not assigned genotypes (i.e. not called).	
}

\seealso{
\code{\link{CC.anypop}}: ClusterCall for arbitrary populations
}

\references{
Schmitz Carley et al. (2017) Automated tetraploid genotype calling by hierarchical clustering. Theoretical and Applied Genetics 130:717–726 DOI: 10.1007/s00122-016-2845-5. 
}

\examples{
##see Vignette at http://potatobreeding.cals.wisc.edu/software
}