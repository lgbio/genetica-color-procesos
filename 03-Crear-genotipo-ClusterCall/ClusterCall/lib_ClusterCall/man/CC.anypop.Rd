\name{CC.anypop}
\alias{CC.anypop}
\title{ClusterCall for any population}
\usage{
CC.anypop(train, predict = NULL, n.core = 1, min.train=NULL, max.missing = 0.05,
  max.range = 0.25, min.posterior = 0.7, min.sep = 0.05, hom.theta=0.15, impute = TRUE)
}
\arguments{
\item{train}{List of objects of class \code{\linkS4class{bipop}}}

\item{predict}{Object of class \code{\linkS4class{pop}}}

\item{n.core}{Number of cores to use for multicore processing}

\item{min.train}{Minimum number of bipops with genotype calls at a given marker for the marker to be predicted. When \code{min.train=NULL}, the marker must be called in each bipop in which it appears, and segregate in at least one of them.}

\item{max.missing}{Maximum proportion of samples assigned NA genotype}

\item{max.range}{Maximum allowed spread of a cluster}

\item{min.posterior}{Minimum posterior membership probability of the sample to its assigned cluster}

\item{min.sep}{Minimum distance between the medians of neighboring clusters}

\item{hom.theta}{Maximum distance between theta = 0 and the center of a cluster with genotype = 0, or between theta = 1 and a cluster with genotype = 4}

\item{impute}{When TRUE, samples that do not meet the min.posterior criterion are still assigned a genotype based on their cluster.}
}
\description{
Generates consensus marker genotypes across biparental families, optional prediction of genotypes for an arbitrary population.
}

\value{
	An object of class \code{\linkS4class{anypop}}
}

\details{
The consistency of genotype calls across two or more biparental populations is determined for markers called in all families (where present) and segregating in at least one family when min.train = NULL. Alternatively, min.train can be set to an integer value to specify the minimum number of families for which a marker must be called to be considered. 

All theta values are clustered as a single population by hierarchical clustering, and the genotype assignment for each cluster is the most prevalent genotype (mode). Concordance is defined as the proportion of the samples with genotye equal to the mode. The number of clusters and genotype-to-cluster assignment is chosen to maximize concordance.

The biparental families constitute a training population that can be used to predict genotypes for an arbitrary population (e.g. diversity or elite panel). A prediction can be made by clustering the training and prediction samples together, and the number of clusters is chosen to maximize concordance in the training set. Cluster range and separation thresholds are also applied.
}
\seealso{
\code{\link{CC.bipop}}: ClusterCall for biparental families
}

\references{
Schmitz Carley et al. (2017) Automated tetraploid genotype calling by hierarchical clustering. Theoretical and Applied Genetics 130:717–726 DOI: 10.1007/s00122-016-2845-5. 
}

\examples{
##see Vignette at http://potatobreeding.cals.wisc.edu/software
}