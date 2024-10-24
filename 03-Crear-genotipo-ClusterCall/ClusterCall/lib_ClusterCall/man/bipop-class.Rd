\docType{class}
\name{bipop-class}
\alias{bipop-class}

\title{Class \code{bipop}}
\description{
Extends the class \code{\linkS4class{pop}} by the addition of genotype calls for a biparental family, using the function \code{\link{CC.bipop}}.
}
\section{Slots}{
\describe{
	\item{\code{geno}}{Marker by sample genotype matrix}
	\item{\code{theta}}{Marker by sample matrix of theta values}
	\item{\code{info}}{Data frame of marker segregation information with columns
		\itemize{
		\item \code{marker}: {Marker name}
		\item \code{called}: {TRUE/FALSE whether the marker is called}
		\item \code{negLog10p}: -{log10(p-value) from the chi-squared goodness-of-fit test of the resulting genotypes to the expected genotype ratios, given the genotypes assigned to parental samples}
	}
	}
	\item{\code{parent1}}{Name of parent 1}
	\item{\code{parent2}}{Name of parent 2}
	\item{\code{r}}{Marker by sample matrix of r values (optional)}
}
}
\references{
Schmitz Carley et al. (2017) Automated tetraploid genotype calling by hierarchical clustering. Theoretical and Applied Genetics 130:717–726 DOI: 10.1007/s00122-016-2845-5. 
}
\examples{
showClass("bipop") }