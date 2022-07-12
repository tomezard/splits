`subsDiag` <-
function (X , ncl, clustMethod="hc", nSim = 2000, sigLvl=0.05, status=TRUE)
	{
	separatedXX <- mhXX <- vector("list", ncl)
	odd <- c()
	if (clustMethod=="hc") {Xc <- cutree(hclust(dist(X)),ncl)}	
	if (clustMethod=="km") {Xc <- kmeans(X,ncl)[[1]]}		
	for (i in 1:ncl) 
		{
		separatedXX[[i]] <- X[Xc==i,]
		if (dim(separatedXX[[i]])[1]==1) warning("you have >=1 cluster with only one element")
		mhXX[[i]] <- mahalanobis(separatedXX[[i]],colMeans(separatedXX[[i]]),cov(separatedXX[[i]]))
		odd <- c(odd, which(mhXX[[i]] > qchisq((1-sigLvl/2),df=dim(separatedXX[[i]])[2])))
		}

	cat("generating reference data sets\n")
	for(ii in 1:(nSim/dim(X)[1]))
		{
		Xdot<-X-matrix(colMeans(X),nrow=dim(X)[1],ncol=dim(X)[2],byrow=TRUE)
		Xdotdot <- Xdot %*% svd(Xdot)$v		#do svd, assume that X* = UDV'; i.e. U%*%D%*%t(V)
		Y<-matrix(0,nrow=dim(X)[1],ncol=dim(X)[2])
		for (i in 1:dim(Y)[2]) {Y[,i]<-runif(dim(Y)[1],min=min(Xdotdot[,i]),max=max(Xdotdot[,i]))}
		ifelse (ii==1, refX <- Y, refX <- rbind(refX,Y))
		}
	cat("jackknifing observed data\n")
	for (i in 1:dim(X)[1])
		{
		infI <- abs((prcomp(X, scale = TRUE)$sdev^2 - prcomp(X[-i,], scale = TRUE)$sdev^2)*(dim(X)[1]-1))
		ifelse (i==1, obsInf <- infI, obsInf <- rbind(obsInf,infI))
		}
	cat("jackknifing reference data\n")
	for (a in 1:dim(refX)[1])
		{
		infI <- abs((prcomp(refX, scale = TRUE)$sdev^2 - prcomp(refX[-a,], scale = TRUE)$sdev^2)*(dim(refX)[1]-1))
		ifelse (a==1, empInf <- infI, empInf <- rbind(empInf,infI))
       	if (status == TRUE)
       		{
			if (a%%(dim(refX)[1]/10) == 0 & dim(refX)[1] >= 10 & a >= 1) 
				{
   	            pct <- (a%/%(dim(refX)[1]/10)) * 10
       	        cat(pct, "%  ", sep = "")
           		}
           	if (a == dim(refX)[1]) {cat("\n\n")}
       		}
		}

	badd <- NULL
	for (i in 1:dim(X)[2]) {badd <- c(badd,which(sort(empInf[,i])[nSim*(1-sigLvl/2)] < obsInf[,i]))}
	toCheck <- list(both=intersect(unique(badd),unique(odd)),influence=unique(badd),distance=unique(odd))
	for (i in 1:3) {if (length(toCheck[[i]])==0) {toCheck[[i]] <- NA}}
	return (toCheck)
	}

