`ddwtGap` <-
function (XX, maxClust=11, nTruth=20, nRep=100, status=TRUE, method="hc", genRndm="pc")
	{
	getWbarK<-function (kcl,space=XX,method=method) 
		{
		Wbark<-0
		if(method=="km") {ifelse (kcl!=1, clustDefn<-kmeans(space,kcl)$cluster, clustDefn<-rep(1,length(space[,1])))}
		if(method=="hc") {ifelse (kcl!=1, clustDefn<-cutree(hclust(dist(space)),kcl), clustDefn<-rep(1,length(space[,1])))}	
		for (m in 1:kcl)
			{
			Dmat<-as.matrix(dist(subset(space,clustDefn==m)))
			nm<-length(subset(space,clustDefn==m))
			wbarM<-sum(Dmat)/(2*nm*(nm-1))
			Wbark<-Wbark+wbarM	
			}
		return(Wbark)
		}
	genRefdata<-function(XX,genRndm=genRndm,nrep=nrep) 
		{
		BB<-vector("list", nrep)			#generate nrep reference data sets.
		if (genRndm=="pc") {
			for (b in 1:nrep) {
				Xdot<-XX-matrix(colMeans(XX),nrow=dim(XX)[1],ncol=dim(XX)[2],byrow=TRUE)
				Xdotdot <- Xdot %*% svd(Xdot)$v		#do svd, assume that X* = UDV'; i.e. U%*%D%*%t(V)
	
				Y<-matrix(0,nrow=dim(XX)[1],ncol=dim(XX)[2])
				for (i in 1:dim(Y)[2]) {Y[,i]<-runif(dim(Y)[1],min=min(Xdotdot[,i]),max=max(Xdotdot[,i]))}
				BB[[b]]<-Y%*%t(svd(Xdot)$v)		#final reference data, i.e. Z
				}					#this is /pc;
			}
		if (genRndm=="uni") {
			for (b in 1:nrep) {
				Y<-matrix(0,nrow=dim(XX)[1],ncol=dim(XX)[2])
				for (i in 1:dim(Y)[2]) {Y[,i]<-runif(dim(Y)[1],min=min(XX[,i]),max=max(XX[,i]))}
				BB[[b]]<-Y
				}}
		return(BB)
		}
		
	AA1<-c();AA2<-c();AA3<-c()
	wGap  <- DDwGap <- seGap <- matrix(0,nrow=nTruth,ncol=maxClust)
	XX <- as.matrix(XX)

	for (a in 1:nTruth) 
		{									
												##THESE STEPS ARE DESCRIBED BY YAN & YE 2007
		Wbarks<-numeric(maxClust)					###STEP 1: compute the value of Wbark for all K
		for (k in 1:maxClust) {Wbarks[k]<-getWbarK(k,XX,method)}

		BB<-genRefdata(XX,genRndm,nRep)			###STEP 2a: generate nrep reference data sets
		
		WbarBK<-matrix(0,nrow=nRep,ncol=maxClust)	###STEP 2b: cluster each reference data set to obtain WbarBK
		for (k in 1:maxClust) for (b in 1:nRep) {WbarBK[b,k]<-getWbarK(k,BB[[b]],"hc")}
		
		gapBARnk<-colMeans(log(WbarBK))-log(Wbarks)
	
		lbar<-colMeans(log(WbarBK))					###STEP 3: calculate the standard deviation
		sdk<-colMeans((log(WbarBK)-matrix(colMeans(log(WbarBK)),ncol=maxClust,nrow=nRep,byrow=TRUE))^2)^.5
		sk<-sdk*sqrt(1+(1/nRep))

		Ghat   <- min(which(gapBARnk[1:(maxClust-1)]>=(gapBARnk[2:maxClust]-sk[2:maxClust])))
		DgpNK  <- gapBARnk[2:maxClust]-gapBARnk[1:(maxClust-1)]
		DDgpNK <- DgpNK[1:(maxClust-1)]-DgpNK[2:maxClust]

		wGap[a,] <- gapBARnk
		seGap[a,] <- sk
		DDwGap[a,] <- c(NA,DDgpNK)

        if (status == TRUE)
        	{
			if (a%%(nTruth/10) == 0 & nTruth >= 10 & a >= 1) 
				{
                pct <- (a%/%(nTruth/10)) * 10
                cat(pct, "%  ", sep = "")
            	}
            if (a == nTruth) {cat("\n\n")}
        	}

		AA1<-c(AA1,Ghat)
		DFFS<-(-(gapBARnk[1:(maxClust-1)]-(gapBARnk[2:maxClust]-sk[2:maxClust])))
		AA3 <- c(AA3,which(DFFS==max(DFFS))+1)
		}

	output <- list(mnGhatWG=mean(AA1),allGhatWG=AA1,mnGhatDD=mean(AA3),allGhatDD=AA3,wGap=wGap,seGap=seGap,DDwGap=DDwGap)
	return(output)
	}
	