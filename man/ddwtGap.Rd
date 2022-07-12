\name{ddwtGap}
\alias{ddwtGap}
\alias{getWbarK}
\alias{genRefdata}
  
\title{Identify the Number of Clusters in a Data Matrix.}
\description{
Calculate the number of clusters using adapted versions of the weighted-gap (wtGap, Tibshirani et al., 2001) and double-difference weighted-gap (DDwtGap, Yan & Ye, 2007) algorithms.  Two clustering methods are implemented to identify the 'true' structure of the data: kmeans and hierarchical clustering.  See \code{details} below.
}
\usage{
ddwtGap(XX, maxClust=11, nTruth = 20, nRep = 100, 
     status = TRUE, method = "hc", genRndm = "pc")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{XX}{A matrix of measurements.  Each row corresponds to one observation.}
\item{maxClust}{The maximum number of clusters to calculate the statistics for.  Note that the DD weighted gap only returns values between 2 and (maxClust-1) clusters.}
\item{nTruth}{Number of times to repeat procedures with different seeds used to calculate the so-called truth.}
\item{nRep}{The number of reference data sets used to calculate standard errors around the weighted-gap statistic.}
\item{status}{Should progress updates be returned?  Defaults to TRUE.}
\item{method}{Which clustering method should be employed?  Defaults to \code{hc} for hierarchical clustering; use \code{km} for kmeans.}
\item{genRndm}{How should the random reference data be generated?  Defaults to \code{pc} for principal components, which simulates data that approximates the distribution of the original data in XX; use "uni" for uniform random data across the range of observations in XX.}
}
\details{

Cluster analysis quantifies similarities and dissimilarities among observations. There are two broad-classes of cluster analysis: distance-type methods, which are based upon pairwise distances between observations (e.g. Tibshirani et al., 2001; Yan and Ye, 2007), and model-based methods, which describe the variance in the sample statistically (e.g. Fraley and Raftery, 2002; see \code{\link[mclust]{Mclust}}).  Both pivot on calculating within-group variability.  Tight groups have low within-group variability; loose groups have high within-group variability.  The optimal delimitation for a sample has low intra- and high intergroup variability.  The statistics implemented in this function are distance-based.

If k groups minimize within-group variability, k is considered the best approximation to the true structure of the sample (kSTAR).  Identification of kSTAR is a non-trivial problem however, since additional groups likely reduce within-group variability until each individual is isolated in its own group.  Tibshirani et al. (2001) developed the gap statistic to test null (H0: homogeneous, non-grouped data) versus alternative hypotheses (H1: heterogeneous data containing more than one group).  The "gap" is that between observed and expected summed distance (in Euclidean space) between individuals.  It can however overestimate kSTAR (Dudoit and Fridlyand, 2002).  To overcome this bias, Yan and Ye (2007) defined the (unbiased) weighted-gap statistic and double-difference weighted-gap statistic, with the logic being that kSTAR is identified if k is much better fit than (k-1) and not significantly worse than (k+1).  The methods are employed successively: the weighted-Gap assesses H0 vs H1 and the double-difference weighted-gap if and only if H0 is rejected (Yan and Ye 2007).

The significance of improvements in fit was assessed using Monte Carlo simulation, drawing values either uniformly over the total range of the data or using using singular-value decomposition.  Tibshirani et al. (2001) consider these as data-independent and data-dependent geometries, respectively.

Yan and Ye (2007) used one 'true' data set.  In an ideal world, allGhat below shows no variance.  Since different starting points for clustering algorithms can generate different clusters however, our approach repeats calculation \code{nTruth} times. 

Yan & Ye (2007) used \code{link{kmeans}}, which (assuming two-dimensions) is biased towards identification of circular groups (i.e. equal radii in all dimensions, Tibshirani and Walther, 2005).  In reality, species may ordinate themselves in non-circular morphospace.  Hierarchical methods are more appropriate in the case of (partially-) nested designs, such as individuals within species: broader groups contain smaller ones nearer to the tips of the dendrogram.  Similarity is typically measured using Euclidean distance.  To obtain k groups using hierarchical clustering \code{link{hclust}}, detail beneath the (k-1)th split is neglected.  Other methods are available, note especially the Bayesian approaches outlined by Fraley and Raftery (2002) in \code{\link[mclust]{Mclust}}, which, amongst others, allows cluster shapes to have a multitude of different shapes (e.g. egg-shaped clusters).}


\value{
A list containing:
  \item{mnGhatWG}{The mean number of well separated clusters based on wGap.}
  \item{allGhatWG}{all Ghats: the number of well separated clusters from each simulation for each truth based on wGap.}
  \item{mnGhatDD}{The mean number of well separated clusters based on DDwGap.}
  \item{allGhatDD}{all Ghats: the number of well separated clusters from each simulation/for each truth based on DDwGap.}
  \item{wGap}{A matrix of weighted gap statistics (Tibshirani et al., 2001) for each truth and number of clusters.}
  \item{seGap}{The standard errors associated with each wGap entry.}
  \item{DDwGap}{The double difference weighted gap statistic (Yan and Ye, 2007) for each truth and number of clusters.}
}

\references{
Dudoit, S., and J. Fridlyand. 2002. A prediction-based resampling method for estimating the number of clusters in a dataset. Gen. Biol. 3::research0036.1-0036.21.

Ezard, T.H.G., Pearson, P.N. & Purvis, A. 2010.  Algorithmic Approaches to Delimit Species in Multidimensional Morphospace. BMC Evol. Biol. 10: 175, doi:10.1186/1471-2148-10-175.

Fraley, C., and A. E. Raftery. 2002. Model-based clustering, discriminant analysis, and density estimation. J. Amer. Stat. Ass. 97:611-631.

Oh, M.-S., and A. E. Raftery. 2007. Model-Based Clustering with Dissimilarities: A Bayesian Approach. J. Comp. Graph. Stat. 16:559-585.

Tibshirani, R., G. Walther, and T. Hastie. 2001. Estimating the number of data clusters via the gap statistic. J. Roy Stat. Soc. B 63: 411-423.

Yan, M., and K. Ye. 2007. Determining the Number of Clusters Using the Weighted Gap statistic. Biometrics 63:1031-1037.
}
\author{Thomas H.G. Ezard \email{tomezard [at] gmail [dot] com}}
\seealso{\code{\link{dimReduct}}, \code{\link{subsDiag}}}
\examples{
#############################################################
##GENERATE RANDOM DATA
x1  <- rnorm(50,0.5,0.4)  ;  y1 <- rnorm(50,mean=0.5,sd=0.4)
x2  <- rnorm(50,1.5,0.4)  ;  y2 <- rnorm(50,mean=1.5,sd=0.4)
sp1 <- matrix(c(x1,x2,y1,y2),ncol=2,byrow=FALSE)
##PLOT DATA
par(mfrow=c(3,1))
plot(sp1,pch=c(rep(1,50),rep(16,50)))

##IDENTIFY CLUSTERS
eg1 <- ddwtGap(sp1)
##PLOT WEIGHTED-GAP AND DD WEIGHTED-GAP
plot(colMeans(eg1$wGap),pch=15,type='b',xlab="Number of Clusters", ylab="wtGap",
    ylim=extendrange(colMeans(eg1$wGap),f=0.2))

arrows(1:dim(eg1$wGap)[2],colMeans(eg1$wGap)-colMeans(eg1$seGap), 
    1:dim(eg1$wGap)[2], colMeans(eg1$wGap)+colMeans(eg1$seGap),length=0.05,angle=90,code=3)
plot(colMeans(eg1$DDwGap),pch=15,type='b',ylim=extendrange(colMeans(eg1$DDwGap),f=0.2),
    xlab="Number of Clusters", ylab="DDwtGap")

}
