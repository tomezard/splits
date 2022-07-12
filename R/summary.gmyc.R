`summary.gmyc` <-
function(object, second.peak=FALSE, ...) {
	#res = result of GMYC
	#display summary of GMYC; likelihood values, chi-square test, estimated parameters, etc...
	
		if (second.peak==TRUE) {
		tmp<-table(cummax(object$likelihood))
		lik.peaks<-names(tmp[tmp>20])
		peak<-which(object$likelihood==lik.peaks[(length(lik.peaks)-1)])}

	
	cat("Result of GMYC species delimitation\n")
	cat("\n\tmethod:\t", object[["method"]], sep="")
	cat("\n\tlikelihood of null model:\t", object$likelihood[1], sep="")
		if (second.peak==FALSE) {
			cat("\n\tmaximum likelihood of GMYC model:\t", max(object$likelihood), sep="")} else
			{cat("\n\tmaximum likelihood of GMYC model:\t", object$likelihood[peak], sep="")}
	
	#chisq test
	if (second.peak==FALSE) {
				LR <- 2*(max(object$likelihood)-object$likelihood[1])} else
				{LR <- 2*(object$likelihood[peak]-object$likelihood[1])}
	cat("\n\tlikelihood ratio:\t", LR, sep="")

	if (object[["method"]] == "single") {
		pvalue <- 1-pchisq(LR, 3)
	}  else if (object[["method"]] == "multiple" || object[["method"]] == "exhaustive") {
		pvalue <- 1 - pchisq(LR, 3 + length(object$threshold.time[[which.max(object$likelihood)]]) - 1)
	}
	
	cat("\n\tresult of LR test:\t", pvalue, ifelse(pvalue<0.001, "***", ifelse(pvalue<0.01, "**", ifelse(pvalue<0.05, "*", "n.s."))), sep="")
	
		if (second.peak==FALSE) {
	cat("\n\n\tnumber of ML clusters:\t", object$cluster[which.max(object$likelihood)], sep="")
		tmp<-object$cluster[object$likelihood>(max(object$likelihood)-2)]
		cat("\n\tconfidence interval:\t", paste(min(tmp),max(tmp),sep="-"), sep="")
	cat("\n\n\tnumber of ML entities:\t", object$entity[which.max(object$likelihood)], sep="")
		tmp<-object$entity[object$likelihood>(max(object$likelihood)-2)]
		cat("\n\tconfidence interval:\t", paste(min(tmp),max(tmp),sep="-"), sep="")

	if (object[["method"]] == "single") {	
		cat("\n\n\tthreshold time:\t", object$threshold.time[which.max(object$likelihood)], "\n", sep="")
	} else if (object[["method"]] == "multiple" || object[["method"]] == "exhaustive") {
		cat("\n\n\tthreshold time:\t", object$threshold.time[[which.max(object$likelihood)]], "\n", sep=" ")
	}	
	cat("\n")} else
	
	{cat("\n\n\tnumber of ML clusters:\t", object$cluster[peak], sep="")
	cat("\n\tnumber of ML entities:\t", object$entity[peak], sep="")
	if (object[["method"]] == "single") {	
		cat("\n\tthreshold time:\t", object$threshold.time[peak], "\n", sep="")
	} else if (object[["method"]] == "multiple" || object[["method"]] == "exhaustive") {
		cat("\n\tthreshold time:\t", object$threshold.time[[peak]], "\n", sep=" ")
	}	
	cat("\n")}
}

