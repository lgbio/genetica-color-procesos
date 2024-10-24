\name{inspect.marker}
\alias{inspect.marker}
\title{Inspect marker genotypes}
\usage{
inspect.marker(called, marker, ...)
}
\arguments{
\item{called}{One or more objects of class \code{\linkS4class{bipop}} and/or \code{\linkS4class{anypop}}}

\item{marker}{Name or index number of marker to plot}

\item{...}{Additonal parameters to pass to plotting function}
}
\description{
Plots genotype calls for a designated marker from a \code{\linkS4class{bipop}} or \code{\linkS4class{anypop}} object.
}
\details{
Sample theta values are plotted against sample index (corresponding to sample column order) with sample genotypes as plot characters (i.e. point symbols; \code{pch}). Red points indicate samples with unassigned (NA) genotypes. 

To plot several populations in one plot, provide a list of \code{\linkS4class{bipop}} and/or \code{\linkS4class{anypop}} objects. For instance, \code{inspect.marker(list(pop1, pop2), 9)} will plot genotypes for populations 1 and 2 consecutively for the 9th marker in the combined data. The listed objects are combined with the function \code{\link{alter.pop}} prior to plotting. 

}

\references{
Schmitz Carley et al. (2017) Automated tetraploid genotype calling by hierarchical clustering. Theoretical and Applied Genetics 130:717–726 DOI: 10.1007/s00122-016-2845-5. 
}

\examples{
##see Vignette at http://potatobreeding.cals.wisc.edu/software
}