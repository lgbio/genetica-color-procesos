
\name{write.pop}
\alias{write.pop}
\title{Write genotype calls to file}
\description{
Writes marker x sample genotype matrix to file.
}
\usage{
write.pop(result, file, thresh = NA)
}
\arguments{
\item{result}{\code{\linkS4class{anypop}} or \code{\linkS4class{bipop}} from which the genotype matrix will be written}

\item{file}{File name for comma-separated values (csv) file}

\item{thresh}{Concordance threshold above which markers are included (applies only to \code{\linkS4class{anypop}} objects)}
}
\details{
Writes a comma-separated values (csv) file of the marker matrix. A concordance threshold for selecting markers can be applied to \code{\linkS4class{anypop}} objects.
}

\references{
Schmitz Carley et al. (2017) Automated tetraploid genotype calling by hierarchical clustering. Theoretical and Applied Genetics 130:717–726 DOI: 10.1007/s00122-016-2845-5. 
}

\examples{
##see Vignette at http://potatobreeding.cals.wisc.edu/software
}