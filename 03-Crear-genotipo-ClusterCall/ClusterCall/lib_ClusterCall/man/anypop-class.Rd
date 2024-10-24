\docType{class}
\name{anypop-class}
\alias{anypop-class}
\title{Class \code{anypop}}
\description{
Extends the class \code{\linkS4class{pop}} by the addition of genotype calls for a population, using the function \code{\link{CC.anypop}}. 
}
\section{Slots}{

\describe{
\item{\code{geno}}{Marker by sample genotype matrix}

\item{\code{concordance}}{The proportion of the population with genotype call equal to the mode of its cluster (vector of values for every marker). }

\item{\code{pops}}{List identifying the component populations}

\item{\code{theta}}{Marker by sample matrix of theta values}

\item{\code{r}}{Marker by sample matrix of r values (optional)}
}}

\references{
Schmitz Carley et al. (2017) Automated tetraploid genotype calling by hierarchical clustering. Theoretical and Applied Genetics 130:717–726 DOI: 10.1007/s00122-016-2845-5. 
}

\examples{
showClass("anypop")
}