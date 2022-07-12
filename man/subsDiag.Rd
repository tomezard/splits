\name{subsDiag}
\alias{subsDiag}

\title{Apply two types of diagnostics to clustered data}
\description{Calculate diagnostics on the subspace identified by cluster analysis}

\usage{
subsDiag(X, ncl, clustMethod = "hc", nSim = 2000, sigLvl = 0.05, status = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{The Data.}
  \item{ncl}{The Structure of the data, obtained using a clustering statistic or some other hypothesis, e.g. \code{\link{ddwtGap}}.}
  \item{clustMethod}{Is the cluster definition obtained using hierarchical clustering "hc" or k-means "km".  See details in \code{\link{ddwtGap}} or on the dedicated help-pages.}
  \item{nSim}{The number of simulations used for Monte Carlo estimates of significance.}
  \item{sigLvl}{The significance level for the chi-squared testing whether observations are significantly, or otherwise, influential on the structure of the data.}
  \item{status}{Report the status of the functions?}
}

\details{
Model diagnostics assess the validity of particular assumptions.  Application of the model diagnostics requires at least two individuals within each well-separated group; the cluster identification algorithms can identify isolated individuals as whole groups.  Depending upon the circumstances, it might be reasonable to consider such individuals suspicious.  The diagnostics aimed to identify individuals that (a) were extreme in measurement and (b) affected significantly the definition of the data structure.

Brooks (1994) calculated the influence of each data point by jack-knifing, i.e. by comparing the dominant eigenvalues of the data with and without a focal observation.  A large difference in dominant eigenvalues implies that the focal observation exerts large influence in the sample, whose significance can be assessed using Monte Carlo estimates.  If variables are not normally distributed , reference data sets can be generated using singular-value decomposition.

Fung (1999) devised a method to identify extreme observations outside the expected range of a particular sample.

Note that data need not be questionable or unusual to exert large influence.
}

\value{
A list containing:
  \item{both}{The index of observations that are BOTH infleuntial and extreme.}
  \item{influence}{The index of infleuntial observations.}
  \item{distance}{The index of extreme observations.}
}
\author{Thomas H.G. Ezard \email{tomezard [at] gmail [dot] com}}
\references{
Brooks, S. P. 1994. Diagnostics for Principal Components: Influence Functions as Diagnostic Tools. The Statistician 43:483-494.
Ezard, T.H.G., Pearson, P.N. & Purvis, A. 2010.  Algorithmic Approaches to Delimit Species in Multidimensional Morphospace. BMC Evol. Biol. 10: 175, doi:10.1186/1471-2148-10-175.
Fung, W.-K. 1999. Outlier Diagnostics in Several Multivariate Samples. The Statistician 48:73-84.
}
\seealso{\code{\link{dimReduct}}, \code{\link{ddwtGap}}}
\examples{
##following the example in ddwtGap ....
data(iris)
subsDiag(as.matrix(iris[,1:4]), 3)

}
