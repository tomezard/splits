`compare` <-
function(test1,test2) {
	chi.sq<-abs(2*(max(test1$likelihood)-max(test2$likelihood)))
	dfs<- 3*(length(test2$threshold.time[[which.max(test2$likelihood)]])-1)
	cat("Comparison of single and multiple threshold GMYC\n")
	cat("\tChi.sq",chi.sq,"Degrees of freedom",dfs,"\n", sep=" ")
	cat("\tSignificance", 1-pchisq(chi.sq,dfs), sep=" ")
}

