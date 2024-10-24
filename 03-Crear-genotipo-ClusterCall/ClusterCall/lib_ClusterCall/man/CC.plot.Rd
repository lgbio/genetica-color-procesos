\name{CC.plot}
\alias{CC.plot}
\title{Plot concordance for population(s)}
\usage{
CC.plot(anypop, names = NULL, xlim = c(1, 0), ...)
}
\arguments{
\item{anypop}{One or more objects of class \code{\linkS4class{anypop}} }

\item{names}{Names for legend when anypop is a list of objects of class \code{\linkS4class{anypop}}, see details}

\item{xlim}{See plotting parameters in base R}

\item{...}{Additonal parameters passed to plotting function}
}
\description{
Plots marker count by concordance threshold for one or more \code{\linkS4class{anypop}} objects produced from \code{\link{CC.anypop}}. 
}

\details{
To plot multiple curves on the same figure, provide a list of objects of class \code{\linkS4class{anypop}}. For instance, \code{CC.plot(anypop=list(anypop1, anypop2, anypop3), names=c("Pop 1", "Pop 2", "Pop 3"))} will plot concordance curves for each of the objects of class \code{\linkS4class{anypop}} with legend labels as specified in the character vector \code{names}.
}

\references{
Schmitz Carley et al. (2017) Automated tetraploid genotype calling by hierarchical clustering. Theoretical and Applied Genetics 130:717–726 DOI: 10.1007/s00122-016-2845-5. 
}

\examples{
##see Vignette at http://potatobreeding.cals.wisc.edu/software
}