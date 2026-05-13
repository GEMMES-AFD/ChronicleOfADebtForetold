calibrateModel<-function(sys,
												 dataParms,
												 dataVariables,
												 dataSeries,
												 dataWeight=NULL,
												 nbOfRuns=1,plotOptim=T,plotFileName=NA,resultFile=NULL,printEvol=TRUE,timeframe){
	allSols <- list()
	########################################
	## PARAMETERS CALIBRATION USING CMAES ##
	########################################
	
	#Getting lower/upper bounds for parms, plus which parms to optimize
	#Extract Min/Max Value for parameters, and if we should calibrate them or keep them fixed
	parmsLower <- dataParms$Min
	parmsUpper <- dataParms$Max
	parmsToOptimize <- as.logical(dataParms$'To Calibrate')
	parmsInit <- dataParms$`Init Value`
	names(parmsInit) <- names(parmsLower) <- names(parmsUpper) <- names(parmsToOptimize) <- dataParms$Name
	
	#Getting the series to be fitted
	targetNames <- dataVariables$Name[which(dataVariables$Fitted==1)]
	target <- as.data.frame(dataSeries[,targetNames])
	names(target) <- targetNames
	## dataMinDist is the data set used for optimization
	## WARNING ## elements in dataMinDist must be named ## WARNING ##
	dataMinDist <- lapply(1:ncol(target), function(i) target[,i])
	names(dataMinDist) <- targetNames
	samplingTimeMinDist <- lapply(dataMinDist, function(x) dataSeries$time)
	names(samplingTimeMinDist) <- names(dataMinDist)
	# pointsWeightMinDist <- NULL #same weight for all observations by default
	pointsWeightMinDist <- dataWeight
	## PARAMETERS FOR CMAES ##
	# Reduce sigma if the algorithm diverges, or increase it to reduce the probability of stopping in a local extrema
	# decrease lambda to reduce computation time, or increase it to reduce the probability of stopping in a local extrema
	# I advise not to change nIterMax
	# Below are default values considered as optimal for optimization of a "reasonably simple" function
	# but on complex problems it is often necessary increase lambda or to adjust sigma to the sensitivity of the system to parameters values
	lambda <- 4*floor(4+3*log(sum(parmsToOptimize))) # number of points per generation
	sigma <- 0.2 # standard deviation at initialization (i.e. initial dispersion of the points)
	nIterMax <- 200 + 50*(sum(parmsToOptimize)+3)**2/sqrt(lambda)       # max number of iterations before termination
	tolCMAES <- 1e-3                                                    # numerical tolerance
	nReDrawMax <- 100 #Do not touch unless you know what you are doing. You might want to increase it in case of a very poorly defined system (returns lots of 1e50)
	
	## CALL THE CPP FUNCTION TO PERFORM OPTIMIZATION
	allPars<-matrix(nrow=(length(parmsInit)+1),ncol=nbOfRuns)
	rownames(allPars)<-c(names(parmsInit),"fit")
	bestFit=Inf
	bestIndex=0

	for(i in 1:nbOfRuns){
	  myTime <- system.time(
	    sol <- cppCMAES(sys=sys, parmsToOptimize=parmsToOptimize, parmsInit=parmsInit,
									  dataMinDist=dataMinDist, samplingTimeMinDist=samplingTimeMinDist, pointsWeightMinDist=pointsWeightMinDist,
									  lambda=lambda, sigma=sigma, nIterMax=nIterMax, tol=tolCMAES, parmsLower=parmsLower, parmsUpper=parmsUpper, nReDrawMax=nReDrawMax,  solver="dopri")
	  )
	  print(myTime)
		# sol <- cppMinDist(sys=sys,
		# 									parms=dataParameters$`Init Value`,
		# 									parmsToOptimize = parmsToOptimize,
		# 									lowerParms = lowerParms,
		# 									upperParms = upperParms,
		# 									dataMinDist=dataMinDist,
		# 									samplingTimeMinDist = samplingTimeMinDist,
		# 									pointsWeight = pointsWeight,
		# 									lambda = lambda,
		# 									sigma = sigma,
		# 									nIterMax = nIterMax,
		# 									tol = tol,
		# 									standardizeParms=TRUE,
		# 									useParallel=TRUE)

		if(sol$dist<bestFit && sol$dist!=0){
			bestFit=sol$dist
			bestIndex=i
		}
		trueParms <- c(sol$parms)
		# print(names(sys$parms))
		names(trueParms) <- names(parmsInit)
		allPars[,i]=c(trueParms,sol$dist)
		resSysCalib <- cppSolve(sys,parms=trueParms[which(names(trueParms) %in% names(sys$parms))])
		# write.csv(resSysCalib,file=paste("solutions/solCon_",i,".csv",sep=""))
		allSols[[i+1]]=resSysCalib
		if(printEvol){
			cat("Distance at the estimated parameters: ", sol$dist,'\n')
			cat(paste(i," of ",nbOfRuns,'\n'))
		}
	}

	if(!is.null(resultFile))
		write.csv(allPars,file=paste(nbOfRuns,resultFile))

	if(plotOptim){

		trueParms <- allPars[,bestIndex]
		names(trueParms) <- names(parmsInit)
		resSysCalib <- cppSolve(sys,parms=trueParms[which(names(trueParms) %in% names(sys$parms))])

		indices<-c(0, seq(1,nrow(target)-1)*10)+1
		simul<-resSysCalib[indices,targetNames]
		if(is.null(dim(targetNames))){
			targetNames=as.vector(targetNames)
			simul<-as.data.frame(simul)
			names(simul)<-targetNames
		}

		nbPlots<-ceiling(length(targetNames)/3)
		for (nb in 1:nbPlots){
			if(!is.na(plotFileName)){
				thisFileName <- paste(strsplit(plotFileName,'\\.')[[1]][1],nb,".png",sep="")
				png(filename=thisFileName,width = 960,height=480)
			}

			remainingPlots=min(3,length(targetNames)-(nb-1)*3)
			par(mfrow=c(2,remainingPlots))
			for(i in 1:remainingPlots){
				index<-(nb-1)*3+i
				name=targetNames[index]
				if (name=="V") {
				  target[1, targetNames[index]] <- 1/trueParms["beta_HUC"]*target[1, "YP"]
				  for (i in 2:12) target[i, targetNames[index]] <- target[i-1, targetNames[index]] + sys$samplesExogVar$Vdot[i]
				}
				matplot(timeframe,cbind(target[,targetNames[index]],simul[targetNames[index]]),type='l',main=targetNames[index],ylab="",xlab="")
				legend("topleft",lty=1:2,col=1:2,legend=c("Observations","Simulations"),bty='n')
			}
			for(i in 1:remainingPlots){
				index<-(nb-1)*3+i
				name=targetNames[index]
				matplot(timeframe,100*cbind(simul[targetNames[index]]-target[,targetNames[index]])/target[,targetNames[index]],type='l',main=targetNames[index],ylab="% erreur",xlab="")
			}
			if(!is.na(plotFileName)){
				dev.off()
			}
		}
	}
	trajectoriesScenar <<- allSols
	return(allPars)
}