\name{gmyc}
\alias{gmyc}
\alias{plot_gmyc}
\alias{summary.gmyc}
\alias{compare}
\alias{remove.terminal.zeros}
\alias{spec.list}

\title{Optimizes genetic clusters using the generalized mixed yule coalescent.}

\description{
This function optimizes either the single threshold (Pons et al. 2006; Fontaneto et al. 2007) or multiple threshold (Monaghan et al. Submitted) versions of the generalized mixed Yule coalescent. Note that the multiple threshold version is in development and can take a long time to run.  We recommend use of the single threshold version. Results are returned as a list object of class 'gmyc'.}

\usage{
gmyc(tr, method = "single", interval = c(0, 5), quiet = FALSE)
}

\arguments{
  \item{tr}{An ultrametric, dichotomous tree object in ape format.}
  \item{method}{Method of analysis, either "single" for single-threshold version or "multiple" for multiple-threshold version.}
  \item{interval}{Upper and lower limit of estimation of scaling parameters, e.g. c(0,10)}
  \item{quiet}{By default shows progress on console. Use quiet = FALSE to disable.}
}

\details{
The function optimizes the likelihood function described in Pons et al. 2006, which specifies the likelihood of branching intervals assuming: i) between species branching according to a Yule model or assuming evenly increasing or decreasing net speciation rates towards the present; ii) coalescence within species according to a neutral coalescent or with evenly increasing or decreasing rates of coalescent branching towards the ancestral node. The single threshold version uses a single threshold to specify the transition from between- to within-species branching. The multiple threshold version starts with a single threshold, but uses an iterative algorithm to test whether relaxation of the single threshold leads to a significant improvement in model fit (described in Monghan et al. submitted). After completing optimization, use plot to see the likelihood surface, a log-lineage through time plot and a tree with clusters indicated. Use summary to summarize the output statistics, including a log-likelihood ratio test of the fitted model against a null model of no distinct species clusters. Use compare(x1,x2) to perform a log likelihood ratio test of a multiple threshold and single threshold solutions for the same tree.
}

\note{
Important: Having the input tree in the correct format - ultrametric, fully dichotomous - is very important as the code assumes these properties. You can use r8s by Sanderson to derive ultrametric trees or simply fit a molecular clock model in PAUP. TreeEdit is useful for manipulating trees in terms of rooting and converting polytomies to dichotomies with zero branch lengths. Use \code{\link[ape]{read.tree}} or \code{\link[ape]{read.nexus}} to read trees into R using ape (>=1.9).}

\value{
gmyc returns an object of class "gmyc": a list with the follow elements
  \item{method }{method used for an analysis}
  \item{likelihood }{likelihood values for each optimization}
  \item{parameters }{estimated parameters for each optimization. (lambda1, lambda2, pp1, pp2)}
  \item{entity }{numbers of entities}
  \item{cluster }{numbers of clusters}
  \item{MRCA }{index of MRCA nodes, i.e. ancestral node of each delimited cluster}
  \item{threshold.time }{The optimized threshold times.}
  \item{tree }{The tree.}
  }

\references{
Fontaneto, D. E. Herniou, C. Boschetti, M. Caprioli, G. Melone, C. Ricci and T.G. Barraclough. 2007. Independently evolving species in asexual bdelloid rotifers. PLoS Biology. 5:e87


Monaghan, M.T., Wild, R., Elliot, M., Fujisawa, T., Balke, M., Inward, D.J.G., Lees, D.C., Ranaivosolo, R., Eggleton, P., Barraclough, T.G., Vogler, AP. Submitted. Accelerated species inventory on Madagascar using coalescent-based models of species delineation.


Pons, J, Barraclough, TG, Gomez-Zurita, J, Cardoso, A, Duran, DP, Hazell, S, Kamoun, S, Sumlin, WD, Vogler, AP,2006. Sequence-based species delimitation for the DNA taxonomy of undescribed insects. Systematic Biology. 55:595-609.

}

\author{Tomochika Fujisawa \email{tomochika.fujisawa05@imperial.ac.uk} and Tim Barraclough \email{t.barraclough@imperial.ac.uk}}

\section{Warning }{The multiple threshold version works, but is very experimental (highly sensitive to initial conditions) and takes a long time to run. Tomochika directs users to the web interface maintained by Jiajie Zhang and to his new Bitbucket repositories.}

\seealso{\code{\link{spec.list}}, \code{\link{test.tr}}}

\examples{
data(test.tr)
ultra_tree <- force.ultrametric(test.tr)

#*#* run the single threshold model
test1 <- gmyc(ultra_tree, method="single", interval=c(0, 10))
summary(test1)  ##show summary results

#*#* plot results
plot_gmyc(test1)

#*#* un-comment to run the multiple threshold model
#*#* see warning and suggestion anove
test2 <- gmyc(ultra_tree, method="multiple", interval= c(0, 10))
summary(test2)
plot_gmyc(test2)

#*#* compare multiple vs single threshold models, former a significant improvement here
compare(test1, test2) 
}

