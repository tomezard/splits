`gmyc` <-
function(tr, method = "single", interval = c(0, 5), quiet = FALSE) {
	#run GMYC
	#contains;
	#reading data
	#applying single threshold
	#applying multiple threshold
	
	#local.env <- globalenv()	#
	local.env <- environment()
	
	### function to read tree and data
	read.data<-function (z=1) {
		#tr<<-tt[[5]]
		#if (length(tt)>1) {tr<<-tt[[z]]} else {tr<<-tt}
	 	bt<- -branching.times(tr)
		#Set any ages after present day or below a certain limit to finite age
	 	bt[bt>-0.000001]<- -0.000001
	 	names(bt)<-NULL
	 	#sb<-sort(bt)
		
		#assing variables to parent environment
		#02/07/08...
		assign("bt", bt, envir=local.env)
		assign("sb", sort(bt), envir=local.env)
		
	 	#19/11/07
	 	##add number of tip and number of all nodes###
		#03/07/08
		##move to assign function
		assign("numnod", length(bt), envir=local.env)
		assign("numtip", length(tr$tip.label), envir=local.env)
		assign("numall", length(bt)+length(tr$tip.label), envir=local.env)
		
		internod<-sb[2:numnod]-sb[1:numnod-1]
	 	internod[numnod]<-0-sb[numnod]
		
		assign("internod", internod, envir=local.env)
		##lists of which nodes are nesting or nested within each node, and arrays with sisters and desc
			  
		#change parameters of sapply(nesting/nested.nodes) -1:numnod -> numtip+1:numall
		assign("nesting", sapply((numtip+1):numall,nesting.nodes), envir=local.env)
		assign("nested", sapply((numtip+1):numall,nest.nodes), envir=local.env)
			
		numnested<-array(NA,numnod)
			  
		numnesting<-array(NA,numnod)
		des<-matrix(NA,nrow=numnod,ncol=2)
		sis<-array(NA,numnod)
			  
		for (j in (1:numnod)) {
			numnested[j] <- length(nested[[j]])
			numnesting[j] <- length(nesting[[j]])
			des[j,]<-as.integer(tr$edge[,2])[as.integer(tr$edge[,1])== j+numtip]
			if (des[j,1] > numtip) { sis[des[j,1]-numtip] <- des[j,2]}	#??????????
			if (des[j,2] > numtip) { sis[des[j,2]-numtip] <- des[j,1]}	#??????????
		}
		
		assign("numnested", numnested, envir=local.env)
		assign("numnesting", numnesting, envir=local.env)
		
		des[des <= numtip] <- NA
		
		assign("des", des, envir=local.env)
		assign("sis", sis, envir=local.env)
		
		##ancestral nodes, column 1, of each node, column 2
			ancs<-cbind(tr$edge[pmatch((1:numnod+numtip),tr$edge[,2]),1],(1:numnod+numtip))
		##node ages for the nodes in ancs
			bt.ancs<-cbind(bt[ancs[,1]-numtip],bt[ancs[,2]-numtip])
		assign("bt.ancs", bt.ancs, envir=local.env)

				
	}

	## Nested nodes
	nest.nodes<-function(x,p=0) {
		#print(paste("nest.nodes() called with x=", as.character(x)))
		numtip <- length(tr$tip.label)
		   
		nods<-array(NA,0)
		desc<-as.integer(tr$edge[,2][tr$edge[,1]==x])
		 
		if (desc[1] > numtip) {
			nods <- c(nods, desc[1], nest.nodes(desc[1]))
		}
		if (desc[2] > numtip) {
			nods <- c(nods, desc[2], nest.nodes(desc[2]))
		}
		   
		if (length(nods)>0) {
			return (nods) 
		} else {
			return(NULL)
		}      
	}

	## Nesting nodes
	nesting.nodes<-function(x,p=0) {  
		#print(paste("nesting.nodes() called with x=", as.character(x)))
		numtip <- length(tr$tip.label)
		
		nod<-array(NA,0)
		  
		##change -1 -> numtip+1
		if (x >= numtip+2)  { 
			anc <- as.integer(tr$edge[,1][tr$edge[,2]==x])
		}  else  { 
			anc <- 1	#?????????????????????? 1 ?? 
		}	
		  
		if (anc	>= numtip+2) {
			nod <- c(nod, anc, nesting.nodes(anc)) 
		}
		if (anc == numtip+1) {
			nod <- c(nod, anc)
		}
		  
		if (length(nod)>0)  {
			return (nod)
		} else {
			return(NULL)
		}
	}
	  

	###Derived from SpeciesTestLikFunctions2005.3

	##replaced global variables with list.
	##02/07/08...
	l.prep<-function () {
		  #print ("l.prep() called")
		  assign("n", length(mrca.nodes), envir=local.env)
		  
		  i.mat<-matrix(0,ncol=numnod,nrow=(n+1))
		  s.nod<-matrix(0,ncol=numnod,nrow=(n+1))
		  
		  assign("nod", nod.type[order(bt)], envir=local.env)
		  
		  ##set of coalescent processes for each mrca node
		  for (i in (1:n)) {
			  s.nod[i,mrca.nodes[i]]<-2
			  
			  #change s.nod[i, -nested[[mrca.nodes[i]]]> s.nod[i, nested[[mrca.nodes[i]][-numtip]
			  if (!is.null(nested[[mrca.nodes[i]]])) { s.nod[i, nested[[mrca.nodes[i]]] - numtip] <- 1}
			  
			  s.nod[i,]<-s.nod[i,order(bt)]	##????
			  i.mat[i,][s.nod[i,]==2]<-2
			  i.mat[i,][s.nod[i,]==1]<-1
			  i.mat[i,]<-cumsum(i.mat[i,])
		  }
		  s.nod[s.nod==2]<-1	###???????????????
		  ##transform number of lineages to coalescent type
		  i.mat<-i.mat*(i.mat-1)
		  
		  
		  #one diversification process
		  s.nod[n+1,]<-nod==0
		  i.mat[n+1,nod==0]<-1
		  i.mat[n+1,nod==2]<--1
		  i.mat[n+1,]<-cumsum(i.mat[n+1,])+1
		  
		  assign("s.nod", s.nod, envir=local.env)
		  assign("i.mat", i.mat, envir=local.env)
		  
	}

	###LIKELIHOOD FUNCTIONS - NULL AND MINIMUM MODEL

	##1) Calc likelihood for null model, generalized Yule model with scaling power, p
	l.null<-function(p=1) {
		#print("l.null() called")
		i.div<-2:(numnod+1)
		i.div<-i.div*(i.div-1)
		
		lambda.div<-numnod/sum(internod*i.div^p)
		assign("lambda.div", lambda.div, envir=local.env)
		
		lik<-i.div^p*lambda.div*exp(-i.div^p*lambda.div*internod)
		
		return (sum(log(lik)))
	}

	#Calcul
	l.min2<-function (q=c(1, 1)) {
	  #print("l.min2() called")

	  ##prob of any event happening
		p<-c(rep(q[1],n),q[2])
		lambda<-sum(s.nod[1:n,])/sum(i.mat[1:n,]^p[-(n+1)]%*%internod)
		lambda<-c(rep(lambda,n),sum(s.nod[n+1,])/i.mat[n+1,]^p[n+1]%*%internod)
				
		b<-t(i.mat^p)%*%lambda
		lik<-b*exp(-b*internod)	
		lambda<-lambda
		
		assign("b", b, envir=local.env)
		assign("lik", lik, envir=local.env)
		assign("lambda", lambda, envir=local.env)
		
		return(sum(log(lik)))
	}

	gmyc.likelihood <- function() {	
		 ###APPLY SIMPLE THRESHOLD MODEL - ASSUME ALL CLUSTERS HAVE SAME BRANCHING PARAMETERS, LAMBDA AND P
		l.prep()
		
		l.results <- rep(NA, 7)
		#x<-optim(c(1, 1),l.min2,method = "BFGS",control=list(fnscale=-1))	#likelihood for GMYC model
		
		#x <- optim(temp.params, l.min2, method = "Nelder-Mead",control=list(fnscale=-1))
		
		x <- optim(c(1, 1), l.min2, method = "Nelder-Mead",control=list(fnscale=-1))
		
		#x <- optim(temp.params, l.min2, method = "L-BFGS", control=list(fnscale=-1))

		
		l.results[1:6]<-c(x$value, lambda[n+1], lambda[1], x$par[2],x$par[1],as.integer(sum(s.nod[n+1,])+1))
		l.results[7]<-n
		#results[j,]<-nod.type
		assign("temp.params",x$par,envir=local.env)
		#print(l.results)
		
		return (l.results)
	}
	
	#exhaustive search of MRCAs...03/12/08
	gmyc.exhaustive <- function() {
		
		#obtaining all combination of MRCAs
		re.combn <- function(node=numtip+1, tr=tr, se=numtip+1) {
			
			parents <- tr$edge[,1]
			children <- tr$edge[,2]

			temp.list <- list()
				
			ch <- children[which(parents==node)]
			se <- se[-(which(se==node))]
			se <- c(se, ch)
			
			temp.list <- c(temp.list, list(se))
			#print(se)
			
			if (any(ch > numtip)) {
				for (cc in ch) {
					temp.temp.list <- temp.list
					if (cc > numtip) {
						for (tl in temp.temp.list) {
							temp.list <- c(temp.list, re.combn(node=cc, tr=tr, se = tl))
						}
					}
				}
			} else {
				
			}
			return(temp.list)
		}
		
		remove.tips <- function(vec) {
			return(vec[vec > numtip])
		}

		MRCA <- lapply(re.combn(tr=tr), remove.tips)
		MRCA <- MRCA[-length(MRCA)]
		#print(MRCA)
		if (!quiet) {
			cat("exhaustive search", "\npossible combinations of MRCAs are", length(MRCA), "\n")
		}
		
		l.results<-c()#matrix(ncol=7, nrow=length(MRCA))
		
		l.mrca <- lapply(MRCA, "-", numtip)
		
		#results<-matrix(nrow=length(MRCA) ,ncol=numnod)

		
		
		assign("temp.params",c(1,1),envir=local.env)
		
		count <- 1
		for (nodes in MRCA) {
			mrca.nodes <- nodes - numtip
			
			assign("mrca.nodes", mrca.nodes, envir=local.env)
					
			#print(mrca.nodes)
			nod.type <- rep(0, numnod)
					
			#print(length(nod.type))
					
			for (j in 1:length(mrca.nodes)) {
				nod.type[mrca.nodes[j]] <- 2
				nod.type[nested[[mrca.nodes[j]]]-numtip] <- 1	#nested[[xxx]] - numtip
			}
			assign("nod.type", nod.type, envir=local.env)
			
			l.results <- rbind(l.results, gmyc.likelihood())
			if (!quiet) {
				#cat("testing nodes", mrca.nodes, "\n")
				cat(l.results[count,1], "\n")
				count <- count+1
			}
		}
		
		return (list(l.results, l.mrca))
		
	}
	
	
	gmyc.single <- function() {
		##SET UP RESULTS MATRICES
		l.results<-matrix(ncol=7,nrow=(nthresh))
		l.mrca <- list()
		results<-matrix(nrow=nthresh,ncol=numnod)

		###WORK OUT RESULTS FOR NULL MODEL - JUST ONE CLUSTER
		#lambda.div <- NULL
		#x<-optimise(l.null,c(0,5),maximum=1)		#likelihood for null model in gmyc.single()
		
		x <- optimise(l.null, interval=interval, maximum=1)
		
		l.results[1,c(6:7,1:2,4)]<-c(as.integer(1),as.integer(1),x$objective, lambda.div, x$maximum)
		
		ml <- 0
		
		if (!quiet) { cat("node\t", "T\t", "loglik\n") }

		stthresh <- 2
		while (sb[stthresh] == sb[1]) {
			stthresh <- stthresh + 1
		}

		assign("temp.params",c(1,1),envir=local.env)

		#stthresh <- 2
		for (j in (stthresh:nthresh)) {
			 
			 ##THIS FOR-LOOP SLIDES THROUGH THE TREE AND DEFINES EACH NODE AGE IN TURN AS THE THRESHOLD
			 ##FOR THE SWITCH BETWEEN CLADOGENESIS AND COALESCENCE WITHIN CLUSTERS
			threshy<-sb[j]
						
			tmp<-(bt.ancs[,1]<threshy)&(bt.ancs[,2]>=threshy)
			nod.type<-tmp+(bt>=threshy)
			mrca.nodes<-which(nod.type==2)
			
			#print(nod.type)
			
			if (nod.type[1]==1) nod.type[1] <- 2
			mrca.nodes <- which(nod.type==2)
			
			assign("mrca", mrca, envir=local.env)
			assign("mrca.nodes", mrca.nodes, envir=local.env)
			assign("nod.type", nod.type, envir=local.env)
			
			l.mrca[[j]] <- mrca.nodes
			
			l.results[j,] <- gmyc.likelihood() 
			
			if (!quiet) { cat (j, threshy, l.results[j,1], "\n") }

		}
		 #results[j,]<-nod.type
		return(list(l.results, l.mrca))
	}
	
	#optimizing multiple MRCAs
	gmyc.multi <- function() {
		assign("temp.params",c(1,1),envir=local.env)
		renew.mrca <- function(n) {	#optimization of MRCA nodes with dividing-fusing algoritm
			#obtain a new set of mrca nodes from the current set.
			#fusing a pair
			f.renew <- function(n) {
				parents <- tr$edge[,1]
				children <- tr$edge[,2]
				
				renewed <- list()
				
				if (length(n) > 1) {
					pair <- combn(n, 2)
					
					for (i in 1:length(n)) {
						parent.node <- parents[children==n[i]]
						sibling.nodes <- children[parents==parent.node]
						sibling <- sibling.nodes[sibling.nodes != n[i]]
						
						if (sibling <= numtip) {
							renewed <- c(renewed, list(c(n[-i], parent.node)))
						} else if (any(n == sibling)) {
							renewed <- c(renewed, list(c(n[-c(i, which(n==sibling))] , parent.node)))
						}
					}
				}
				return (unique(renewed))
			}

			#dividing a node 
			d.renew <- function(n) {
				parents <- tr$edge[,1]
				children <- tr$edge[,2]
				
				renewed <- list()
				if (length(n) > 1) {
					for (i in 1:length(n)) {
						child.nodes <- children[parents==n[i]]
						if (any(child.nodes > numtip)) {
							renewed <- c(renewed, list(c(n[-i], child.nodes[child.nodes > numtip])))
						} 
					}
				}
				return (renewed)
			}

		
			return (c(f.renew(n), d.renew(n)))
		}
		
		##rewritten to avoid redundant part of starting point setting
		##search for 3 part 11/09/08
		##revised 15/09/08
		select.start.re <- function(time, ml) {	#recurive division of starting point of d-f algoritm, avoiding trying parts with less improvement
			result <- c()
			result.mrca <- list()
			
			m <- mean(time)
			
			left <- time[time < m]
			right <- time[time >= m]	
			mid <- time[time >= mean(left) & time < mean(right)]
			
			part <- list(left, mid, right)
			
			start <- c(mean(left), m, mean(right))
			
			if (!any(is.na(start)) || length(unique(start)) == 1) {	#stop recursive division when length of starting points are less than 1 ...03/12/08
				improve <- c(FALSE , FALSE , FALSE)
				num.improve <- c(0, 0, 0)
				
				#print(start)
				
				for (i in 1:length(start)) {
					if (!quiet) {cat("start at", start[i], "\n")}
					
					temp <- renew.likelihood.from.thresh(start[i])
					lik <- temp[[1]]
					mrca <- temp[[2]]
					
					if (any(lik[,1] > ml)) {
						if (!quiet) {cat("improvement found\n")}
						result <- rbind(result, lik[which(lik[,1] > ml),])
						result.mrca <- c(result.mrca, mrca[which(lik[,1] > ml)])
						improve[i] <- TRUE
						num.improve[i] <- length(lik[lik[,1] > ml, 1])
					}
				}
				if (!quiet) {cat(num.improve, "\n")}
				
				if (!all(num.improve == 0)) {
					temp <- select.start.re(part[[which.max(num.improve)]], max(result[,1]))
					result <- rbind(result, temp[[1]])
					result.mrca <- c(result.mrca, temp[[2]])
				}
			}
			
			return(list(result, result.mrca))
			
		}
		
		
	
		
		##rewritten to avoid redundant part of starting point setting
		##search for 3 part 
		select.start <- function(time, ml) {	#recurive division of starting point of d-f algoritm
			result <- c()
			result.mrca <- list()
			
			m <- mean(time)
			
			left <- time[time < m]
			right <- time[time >= m]	
			mid <- time[time >= mean(left) & time < mean(right)]
			
			part <- list(left, mid, right)
			
			start <- c(mean(left), m, mean(right))
			#improve <- c(FALSE , FALSE , FALSE)
			
			for (i in 1:length(start)) {
			
				if (!quiet) {cat("start at", start[i], "\n") }
				
				temp <- renew.likelihood.from.thresh(start[i])
				lik <- temp[[1]]
				mrca <- temp[[2]]
				
				if (any(lik[,1] > ml)) {
					if (!quiet) { print("improvement found\n") }
					result <- rbind(result, lik[which(lik[,1] > ml),])
					result.mrca <- c(result.mrca, mrca[which(lik[,1] > ml)])
					
					if (length(part[[i]]) > 1) {
						if (!quiet) {cat("recursively trying part of", start[i], "\n")} 
						temp <- select.start(part[[i]], max(lik[,1]))
						result <- rbind(result, temp[[1]])
						result.mrca <- c(result.mrca, temp[[2]])
					}
				}
			}
			return(list(result, result.mrca))
			
		}
		
		renew.likelihood.from.thresh <- function(start) {	#find ML solution with renew.mrca and select.start
			l.results <- c()
			
			l.mrca <- list()	######add, 08/09/08
			
			#start <- mean(sb)
			#cat("start: ", start, "\n")
			max.lik <- NA
			mrca<-array(FALSE,numnod)
			mrca[2:numnod]<- (bt[as.integer(tr$edge[,1][tr$edge[,2]>numtip]) - numtip]<start)&(bt[as.integer(tr$edge[,2][tr$edge[,2]>numtip]) - numtip]>=start)	##???
						
			nod.type <- (bt>=start)+mrca

			if (nod.type[1]==1) nod.type[1] <- 2
			mrca.nodes <- which(nod.type==2)
			
			assign("mrca", mrca, envir=local.env)
			assign("mrca.nodes", mrca.nodes, envir=local.env)
			assign("nod.type", nod.type, envir=local.env)
			

			initial.mrca <- mrca.nodes
			while (TRUE) {
				
				found <- FALSE
				max.mrca <- initial.mrca
				for (nodes in renew.mrca(initial.mrca+numtip)) {
					mrca.nodes <- nodes - numtip
					if (mrca.nodes[[1]] != 1) {
						
						assign("mrca.nodes", mrca.nodes, envir=local.env)
						
						#cat(mrca.nodes, "\n")
						nod.type <- rep(0, numnod)
						
						#print(length(nod.type))
						
						for (j in 1:length(mrca.nodes)) {
							nod.type[mrca.nodes[j]] <- 2
							nod.type[nested[[mrca.nodes[j]]]-numtip] <- 1	#nested[[xxx]] - numtip
						}
						assign("nod.type", nod.type, envir=local.env)
					
						#print(c(length(nod.type), numnod))
						#print(nod.type)
						
						res <- gmyc.likelihood()
						#print(length(res))
						
						if (!is.na(max.lik)) {
							if (max.lik <= res[1]) {
								max.lik <- res[1]
								l.results <- rbind(l.results, res)
								max.mrca <- nodes - numtip
								
								l.mrca <- c(l.mrca, list(max.mrca))	######add, 08/09/08
								
								if (!quiet) {
									cat(max.lik, "\n")
								}
								#print("improved")
								found <- TRUE
							}
						} else {
							max.lik <- res[1]
							l.results <- rbind(l.results, res)
							max.mrca <- nodes - numtip
							
							l.mrca <- c(l.mrca, list(max.mrca))	######add, 08/09/08
							
							found <- TRUE
						}
					}
				}
				
				if (found) {
					initial.mrca <-  max.mrca
				} else { 
					if (!quiet) {cat("break\n")}
					break
				}
			}
			
			rownames(l.results) <- NULL
			
			return (list(l.results, l.mrca))
		
		}
		
		
		##SET UP RESULTS MATRICES
		l.results<-c()#matrix(ncol=7,nrow=(nthresh))
		#results<-matrix(nrow=nthresh,ncol=numnod)
		l.mrca <- list()
		
		
		###WORK OUT RESULTS FOR NULL MODEL - JUST ONE CLUSTER
		#lambda.div <- NULL
		#x<-optimise(l.null,c(0,5),maximum=1)	#likelihood of GMYC model in gmyc.multi()
		
		x <- optimise(l.null, interval=interval, maximum=1)
		
		l.results<- rbind(l.results, c(x$objective, lambda.div, NA, x$maximum, NA, as.integer(1),as.integer(1)))	#[1,c(6:7,1:2,4)]
		l.mrca <- c(l.mrca, list(numtip+1))
		
		if (!quiet) {
			cat("null likelihood\n")
			cat(l.results[1, 1], "\n")
		}
		
		###########################
		##get a set of renewed MRCAs for a threshold
		##10/07/08...
		##########################
		#start <- mean(sb)
		#l.results <- rbind(l.results, renew.likelihood.from.thresh(start))
		#l.results <- select.start1(sb, 0)
		
		
		
		#temp <- select.start(sb, l.results[1, 1])	#select.start(sb, 0) ????????
		
		if (!quiet) {
			cat("\nGMYC likelihood\n")
		}
		
		temp <- select.start.re(sb, l.results[1, 1])	
		
		l.results <- rbind(l.results, temp[[1]])
		l.mrca <- c(l.mrca, temp[[2]])
		
		return (list(l.results, l.mrca))
		
	}
	
	#######################
	##PROGRAM ENTRY POINT##
	#######################
	##check if the input tree is appropriate
	if (!is.ultrametric(tr)) {
		stop("Your ultrametric tree is not ultrametric, please check")
	}
	if (!is.binary.tree(tr)) {
		stop("Your input tree is not fully bifurcating, please resolve with zero branch lengths")
	}
		
	##SET UP DATA IN FORMAT USED LATER
	read.data()
		
	ntrees<-1

	##numnod <- -min() -> numnode<-max()-numtip
	numnod <- max(as.integer(tr$edge[,1])) - numtip

	##which is last node to try for threshold model
	nthresh<-numnod
	
	#run analysis
	if (method == "single" || method == "s") {
		method <- "single"
		temp <- gmyc.single()
		l.results <- temp[[1]]
		l.mrca <- temp[[2]]
		
	} else if (method == "multiple" || method == "m") {	
		method <- "multiple"
		temp <- gmyc.multi()
		
		l.results <- temp[[1]]
		l.mrca <- temp[[2]]
	} else if (method == "exhaustive" || method == "e") {
		method <- "exhaustive"
		temp <- gmyc.exhaustive()
		l.results <- temp[[1]]
		l.mrca <- temp[[2]]
	} else {
		stop("Invalid name of optimiztion method. Only single(s) or multiple(m) and exhaustive(e) are available.")
	}
	
	#colnames(l.results)<-c("likelihood","lambda1","lambda2","p1","p2","num entities","num clusters")
	cat("\n", date(), "\n", sep="")

	#plot.cluster1(tr, sb[which.max(l.results[,1])])
	cat("finish.\n")

	#making result data structure...
	result <- list()
	result[["method"]] <- method
	result[["likelihood"]] <- l.results[,1]
	result[["parameters"]] <- l.results[,2:5]
	colnames(result[["parameters"]]) <- c("lambda.div", "lambda.coal", "p.div", "p.coal")
	result[["entity"]] <- l.results[,6]
	result[["cluster"]] <- l.results[,7]
	
	if (method == "single") {
		result[["MRCA"]] <- l.mrca
		result[["threshold.time"]] <- sb
	# } else if (method == "exhaustive") {
		# result[["MRCA"]] <- l.mrca
		# result[["threshold.time"]] <- NA
	} else if (method == "multiple" || method == "exhaustive") {
		
		reduce.threshold <- function(mrcas) {	#function to obtain multiple threshold times from a set of mrcas
 			parent <- tr$edge[,1]
			child <- tr$edge[,2]
		
			thresh.group <- list()
			thresh.time <- c()
		
			mrcas <- mrcas + numtip
			k <- 1
			
			while (TRUE) {
				times <- bt[mrcas-numtip]
				thresh1.time <- min(times)
				thresh1.node <- mrcas[which.min(times)]
				
				mrcas <- mrcas[-which.min(times)]
				
				if (length(mrcas) == 0) { 
					thresh.time <- c(thresh.time, thresh1.time)	###??????????????????????? last MRCA ???
					thresh.group[[k]] <- thresh1.node
					break 
				}	
				
				member <- thresh1.node
				del <- c()
				for (i in 1:length(mrcas)) {
					#print(mrcas[i])
					#print(length(mrcas))
					par.nod <- parent[child==mrcas[i]]
					t.par <- bt[par.nod-numtip]
					#print(c(mrcas[i], par.nod, t.par, length(mrcas)))
					if (t.par < thresh1.time) {
							member <- c(member, mrcas[i])
							del <- c(del, i)
					}
					
				}
				thresh.time <- c(thresh.time, thresh1.time)
				thresh.group[[k]] <- member
				
				k <- k+1
				
				if (length(del) != 0) {	mrcas <- mrcas[-del]}
				
				if (length(mrcas) == 0) { break }	
			}
			
			return (thresh.time)
		}
		
		result[["MRCA"]] <- l.mrca
		result[["MRCA"]] <- l.mrca
		result[["threshold.time"]] <- lapply(l.mrca, reduce.threshold)
	
		#result[["threshold.time"]] <- l.mrca
	}
	#result[["threshold.time"]] <- sb
	result[["tree"]] <- tr
	
	class(result) <- "gmyc"
	
	return (result)
}

