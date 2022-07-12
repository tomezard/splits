\name{dimReduct}
\alias{dimReduct}
\title{Strip uninformative measurements}
\description{
This function is essentially a wrapper around \code{\link[nFactors]{nScree}} and \code{\link[paran]{paran}}, taking a matrix of measurements and returning the predictions from a Principal Components Analysis PCA.  The number of axes in the PCA predictions is determined from either the Kaiser-Guttman criterion \code{KG} or Horn's Parallel Analysis \code{PA} depending on which is specified.
}
\usage{
dimReduct(XX, how = "PA", scale = TRUE, iterations=0, centile=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{XX}{A matrix of measurements.  Each row corresponds to one observation.}
  \item{how}{PA or KG.  defaults to PA}
  \item{scale}{Scale variances, defaults to TRUE}
  \item{iterations}{The number of iterations used for PA.}
  \item{centile}{The centile }
}
\details{
PA is implemented using the excellent \code{\link[paran]{paran}} with default settings, but can be adapted as required, e.g. to apply the revisions to address Glorfeld's claims of over-retention of variables.  See \code{\link[paran]{paran}}.  KG is implemented from the standard deviations of the PCA model output.  

\code{\link{prcomp}} is used to fit the PCA model, which assumes linear measurements.  Transformation may be required prior to the PCA.

Multiple measurements on a focal individual are often highly correlated.  Similarities along one axis might mask differences along another and thus generate inflated estimates of similarity between individuals unless prior weighting is given to certain axes.  Rather than arbitrary assignment of a, say, 0.05 threshold, the number of axes that adequately describe each individual retained can be determined using PA or KG.  Both too many or too few axes of variation can generate bias.

}
\value{
  \item{sS}{A matrix of predicted values from the PCA, with observations as rows and model predictions as columns.}}

\author{Thomas H.G. Ezard \email{tomezard [at] gmail [dot] com}, but the hard work for PA by Alexis Dinno (see \code{\link[paran]{paran}})}

\references{
Dinno, A. in press. Implementing Horns Parallel Analysis for Principal Components Analysis and Factor Analysis. The Stata Journal.

Ezard, T.H.G., Pearson, P.N. & Purvis, A. 2010.  Algorithmic Approaches to Delimit Species in Multidimensional Morphospace. BMC Evol. Biol. 10: 175, doi:10.1186/1471-2148-10-175.

Glorfeld, L. W. 1995. An Improvement on Horn's Parallel Analysis Methodology for Selecting the Correct Number of Factors to Retain. Ed. Psych.Meas. 55, 377-393

Guttman, L. 1954. Some necessary conditions for common factor analysis. Psychometrika 19, 149-162.

Horn, J. L. 1965. A rationale and a test for the number of factors in factor analysis. Psychometrika 30, 179-185.

Kaiser, H. F. 1960. The application of electronic computer to factor analysis. Ed. Pysch. Meas. 20, 141-151.
}

\seealso{\code{\link{prcomp}} and \code{\link{princomp}},\code{\link[paran]{paran}}}

\examples{
data(iris)
dimReduct(iris[,-5])

}
