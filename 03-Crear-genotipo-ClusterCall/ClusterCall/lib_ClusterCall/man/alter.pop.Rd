\name{alter.pop}
\alias{alter.pop}
\title{Combine population objects, merge samples, select samples}
\description{
Combine objects inheriting from class \code{\linkS4class{pop}}, remove or merge samples
}
\usage{
alter.pop(pop, merge = NULL, remove = NULL, keep = NULL)
}
\arguments{
\item{pop}{One or more objects inheriting from class \code{\linkS4class{pop}} (pass multiple objects as a list)}

\item{merge}{List of samples to merge}

\item{remove}{Character vector of samples to eliminate}

\item{keep}{Character vector of samples to retain}
}
\details{
A function for altering the sample contents of objects inheriting from class \code{\linkS4class{pop}}. The function combines populations (if applicable), then merges (i.e. creates consensus) samples, and finally selects samples to return using the list passed to \code{keep} or \code{remove}. 

When multiple objects inheriting from class \code{\linkS4class{pop}} are provided as a list, the populations are first combined to a single object (i.e. \code{pop=list(pop1, pop2)} combines 'pop1' and 'pop2').

A named list of vectors of samples to merge can be provided to the \code{merge} parameter. For example, \code{merge = list(Consensus=c("Rep1", "Rep2"))} creates a consensus individual from samples 'Rep1' and 'Rep2' with the new name 'Consensus'. Additionally, if more than one sample has the same name a consensus sample is automatically created. Consensus values for theta values are median values. For genotypes, the consensus value is the mode genotype or NA (no call) if calls are conflicting with no mode value.
 
Only \code{remove} or \code{keep} can be used in a call to the function; assigning a vector to both parameters will result in an error.

}
\value{
When the \code{pop} input is a list of objects of class \code{\linkS4class{pop}}, an object of class \code{\linkS4class{pop}} is returned. If one or more of the objects is of class \code{\linkS4class{bipop}} or class \code{\linkS4class{anypop}} containing marker calls, the returned object is of class \code{\linkS4class{anypop}}. 

When the input is a single object of the class \code{\linkS4class{pop}}, an object of class \code{\linkS4class{pop}} is returned. 

When the input is a single object of the class \code{\linkS4class{bipop}}, a \code{\linkS4class{bipop}} object will be returned, unless the parental samples (indicated in calling the function \code{\link{CC.bipop}}) are absent, in which case an object of class \code{\linkS4class{anypop}} is returned. Stored marker statistics may be inaccurate with the removal or merging of samples.

When the input is a single object of the class \code{\linkS4class{anypop}}, an object of class \code{\linkS4class{anypop}} is returned. 
	
}
\references{
Schmitz Carley et al. (2017) Automated tetraploid genotype calling by hierarchical clustering. Theoretical and Applied Genetics 130:717–726 DOI: 10.1007/s00122-016-2845-5. 
}

\examples{
##see Vignette at http://potatobreeding.cals.wisc.edu/software
}