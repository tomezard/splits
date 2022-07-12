\name{test.tr}
\alias{test.tr}
\docType{data}
\title{A test tree for use with gmyc.}
\description{
A test tree for use with gmyc
}

\usage{data(test.tr)}
\format{
  The format is:
List of 4
 $ edge       : int [1:298, 1:2] 151 152 153 154 155 156 157 157 156 155 ...
 $ Nnode      : int 149
 $ tip.label  : chr [1:150] "spec1.5" "spec1.4" "spec1.3" "spec1.2" ...
 $ edge.length: num [1:298] 0.0454 0.35302 0.59264 0.00541 0.00243 ...
 
 - attr(*, "class")= chr "phylo"
}
\details{
An example tree for use with \code{\link{gmyc}}}

\author{Tomochika Fujisawa \email{tomochika.fujisawa05@imperial.ac.uk} and Tim Barraclough \email{t.barraclough@imperial.ac.uk}}

\seealso{\code{\link{gmyc}}}

\examples{
data(test.tr)
str(test.tr)

}
