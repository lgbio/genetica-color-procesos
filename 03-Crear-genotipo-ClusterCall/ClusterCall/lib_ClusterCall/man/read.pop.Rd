
\name{read.pop}
\alias{read.pop}
\title{Read in population marker data}
\usage{
read.pop(theta.file, r.file = NULL, error.checking = FALSE,
  method = "complete", thresh = 5, max.missing.theta =0.2,  ...)
}
\arguments{
\item{theta.file}{Comma-separated values (csv) file of theta values (marker x samples)}

\item{r.file}{Comma-separated values (csv) file of theta values (marker x samples)}

\item{error.checking}{Whether to check for duplicate samples. If \code{error.checking=TRUE}, plots dendrogram of sample marker distance.}

\item{method}{Agglomeration method used in hierarchical clustering of samples, see \code{hclust}}

\item{thresh}{Height threshold for identifying replicate samples when \code{error.checking=TRUE}}

\item{max.missing.theta}{Maximum proportion of markers with NA theta values for each sample}

\item{...}{Graphical parameters to be passed to dendrogram plot when \code{error.checking=TRUE}}
}
\value{
When \code{error.checking=FALSE} An object of class \code{\linkS4class{pop}}.

When \code{error.checking=TRUE} A list with [[1]] an object of class \code{\linkS4class{pop}}, [[2]] a list of replicate samples, and [[3]] a list of samples that were removed (see \code{thresh} parameter)
}
\description{
Reads in population theta and r values. Option to identify replicate samples.
}
\details{
Returns an object of class \code{pop} that is input to other functions. Additonal plotting arguments can be passed to adjust the dendrogram plot. Samples with >max.missing.theta proportion of theta values missing do not appear in the resulting \code{pop} object. Identifies and reports duplicate/replicate samples and samples removed for >max.missing.theta proportion missing theta values when \code{error.checking=TRUE}. When \code{error.checking=TRUE}, the dendrogram is plotted and a red, dashed line indicates the \code{thresh} level.
}

\references{
Schmitz Carley et al. (2017) Automated tetraploid genotype calling by hierarchical clustering. Theoretical and Applied Genetics 130:717–726 DOI: 10.1007/s00122-016-2845-5. 
}
\examples{
##see Vignette at http://potatobreeding.cals.wisc.edu/software
}