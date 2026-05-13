source("SourceCodeParallel/utilities.R")


initSysCpp <- function(fileName=NULL, tderiv=NULL, parms=NULL, times=NULL, 
                       y0=NULL, y0Text=NULL, intermediateVar=NULL,
                       samplesExogVar=NULL, 
                       eventTime=NULL, eventVar=NULL, useParallel=FALSE, verboseCMAES=0,
                       returnRK4="all", argsDopri=NULL, longIntForHCMSubdivision="long long unsigned int",
                       compileForPython=FALSE) {

  if(!is.null(y0) && !is.null(y0Text)) {
    warning("Both y0 and y0Text are defined, only y0 will be used.")
  }
  
  if (verboseCMAES!=0 && useParallel!=FALSE) {
    stop('Setting both verboseCMAES>0 and useParallel==TRUE will cause Rcpp to crash ("C Stack usage too close too the limit"). Either set verbose to zero or useParallel to FALSE.')
  }
  
  ## INITIALIZE SYS OBJECT
  # load sys from external file 
  if (!is.null(fileName)) {
    sys <- loadModel(fileName, samplesExogVar)
  } else {
    sys <- NULL
  }
  
  # Below tests might look weird but are correct, do not edit !
  if (is.null(eventTime) || (is.logical(eventTime) && eventTime==FALSE)) sys$eventTime <- NULL
  if (is.null(eventVar) || (is.logical(eventTime) && eventTime==FALSE)) sys$eventVar <- NULL
  if (!is.null(eventTime) && (is.logical(eventTime) && eventTime==TRUE)) eventTime <- NULL
  if (!is.null(eventVar) && (is.logical(eventTime) && eventTime==TRUE)) eventVar <- NULL
  if(!is.null(eventTime) && is.null(sys$eventTime)) eventTime <- NULL
  if(!is.null(eventVar) && is.null(sys$eventVar)) eventVar <- NULL
  
  if ((!is.null(samplesExogVar) && isFALSE(samplesExogVar))) { #Need for a specific treatment of exogenous variables initialization (as they are included in sys$intermediateVar by loadModel)
    if (length(sys$samplesExogVar)>0) {
      for (i in 1:length(sys$samplesExogVar)) {
        sys$intermediateVar <- sys$intermediateVar[-which(names(sys$intermediateVar)==names(sys$samplesExogVar)[i])]
      }
    }
  }
  
  #Complete argsDopri
  if(is.null(argsDopri)) argsDopri <- list()
  if (is.null(argsDopri$atol)) argsDopri$atol <- 1e-4                             #absolute tolerance for step validation
  if (is.null(argsDopri$rtol)) argsDopri$rtol <- 0                                #relative tolerance for step validation
  if (is.null(argsDopri$fac)) argsDopri$fac <- 0.85                               #reduction factor for step length update
  if (is.null(argsDopri$facMin)) argsDopri$facMin <- 0.1                          #min multiplicative factor for step length update
  if (is.null(argsDopri$facMax)) argsDopri$facMax <- 4                            #max multiplicative factor for step length update
  if (is.null(argsDopri$nStepMax)) argsDopri$nStepMax <- 100                      #maximum number of intermediate steps
  if (is.null(argsDopri$hInit)) argsDopri$hInit <- sys$times[2] - sys$times[1]    #initial step length
  if (is.null(argsDopri$hMin)) argsDopri$hMin <- argsDopri$hInit/100              #minimum step length
  if (is.null(argsDopri$hMax)) argsDopri$hMax <- sys$times[2] - sys$times[1]      #maximum step length
  if (is.null(argsDopri$nStiffMax)) argsDopri$nStiffMax <- 100                    #maximum step length
  if (is.null(argsDopri$nStiffSuccessiveMax)) argsDopri$nStiffSuccessiveMax <- 15 #maximum step length
  
  sys <- completeSys(allNames = c("tderiv", "parms",  "times", "y0",
                                  "y0Text", "intermediateVar", "eventTime",
                                  "eventVar",  "useParallel", "returnRK4", "verboseCMAES",
                                  "samplesExogVar", "longIntForHCMSubdivision", "compileForPython", "argsDopri"), 
                     sys, envir=environment())
  missingArguments(c("tderiv", "parms", "times", "returnRK4", "longIntForHCMSubdivision"), sys, envir=environment())

  sys$verboseCMAES <- ifelse(is.numeric(verboseCMAES) && (verboseCMAES %in% 0:3), verboseCMAES, 0)
  
  if (!is.null(sys$compileForPython)) {
    sys$useParallel <- 0 # No parallelization in python (because we only run the solver, not CMAES or any other parallelized algorithm)
  }
  sys$useParallel <- ifelse(is.null(sys$useParallel) || sys$useParallel=="no" || sys$useParallel==FALSE || sys$useParallel==0, 0,1)
  
  sys$returnRK4 <- ifelse(sys$returnRK4=="variable" || sys$returnRK4==0, 0, 
                          ifelse(sys$returnRK4=="derivative" || sys$returnRK4==1, 1, 
                                 ifelse(sys$returnRK4=="intermediateVar" || sys$returnRK4==2, 2, 3)))
  #Initialization
  sys$dim <- length(sys$tderiv)
  
  ## MODEL STRUCTURE VERIFICATION
  checkModeDefinition(sys)
  
  ## FORMATTING FOR C++
  initScipen <- getOption("scipen")
  options(scipen=999) #avoid scientific notation, that can cause issues with the different conversions from R to cpp
  

  #Define strings for the cpp function
  allEquations <- unloopEquations(tderiv = sys$tderiv, 
                                  intermediateVar = sys$intermediateVar)

  # Defines numerical values for y0 from equations in y0Text
  # if (!is.null(sys$y0Text)) sys$y0 <- computeInitialPosition(sys)
  if (!is.null(sys$y0Text)) sys$y0 <- sys$y0Text
  
  # Define cpp code defining the model
  sys$strFunc <- makeCppFunc(allEquations, sys$parms, sys$tderiv, sys$intermediateVar)
  sys$strEventTime <- makeEventTime(sys)
  sys$strEventVar <- makeEventVar(sys)
  
  # Define parameters required at compile time 
  editPreprocRCPP(sys)
  createFuncCpp(sys)
  
  if(sys$useParallel==1) {
    Sys.setenv(CXX_STD="CXX14", PKG_LIBS = "-fopenmp -larmadillo")
  } else {
    Sys.setenv(CXX_STD="CXX14", PKG_LIBS = "-larmadillo")
  }
  sourceCpp("SourceCodeParallel/cppCode/functionsForR.cpp", showOutput = FALSE)
  options(scipen=initScipen)
  sys
}

###########################################################################
## Simulate the previously compiled system (multiple solvers available)  ##
## solver can be "euler", "RK4Fixed" or "dopri" (more will come later)   ##
###########################################################################
cppSolve <- function(sys, y0=NULL, parms=NULL, times=NULL, samplesExogVar=NULL, 
                     solver="dopri", argsDopri=NULL, updateY0=FALSE) {
  sys <- completeSys(allNames=c("parms", "y0", "times", "samplesExogVar", "argsDopri"), sys=sys, envir=environment())
  if (is.null(sys$samplesExogVar)) sys$samplesExogVar <- list()
  missingArguments(c("parms", "times", "y0"), sys, envir=environment())
  if (updateY0==TRUE) computeInitialPosition(sys)
  if (is.null(solver) || !(solver %in% c("euler", "RK4Fixed", "dopri", "dopriStiff"))) {
    stop("Unknown solver. Set solver to \"euler\", \"RK4Fixed\" or \"dopri\" ")
  }

  #if(solver!="dopri" && !is.null(sys$argsDopri))
  #  warning("Arguments for dopri method are specified in argsDopri, but method used is not dopri; argsDopri will be ignored")
  if(is.null(sys$argsDopri)) sys$argsDopri <- list()
  if (is.null(sys$argsDopri$atol)) sys$argsDopri$atol <- 1e-4                            #absolute tolerance for step validation
  if (is.null(sys$argsDopri$rtol)) sys$argsDopri$rtol <- 0                               #relative tolerance for step validation
  if (is.null(sys$argsDopri$fac)) sys$argsDopri$fac <- 0.85                              #reduction factor for step length update
  if (is.null(sys$argsDopri$facMin)) sys$argsDopri$facMin <- 0.1                         #min multiplicative factor for step length update
  if (is.null(sys$argsDopri$facMax)) sys$argsDopri$facMax <- 4                           #max multiplicative factor for step length update
  if (is.null(sys$argsDopri$nStepMax)) sys$argsDopri$nStepMax <- 100                     #maximum number of intermediate steps
  if (is.null(sys$argsDopri$hInit)) sys$argsDopri$hInit <- sys$times[2] - sys$times[1]   #initial step length
  if (is.null(sys$argsDopri$hMin)) sys$argsDopri$hMin <- sys$argsDopri$hInit/100             #minimum step length
  if (is.null(sys$argsDopri$hMax)) sys$argsDopri$hMax <- 10*sys$argsDopri$hInit              #maximum step length
  if (is.null(sys$argsDopri$nStiffMax)) sys$argsDopri$nStiffMax <- 100                   #maximum step length
  if (is.null(sys$argsDopri$nStiffSuccessiveMax)) sys$argsDopri$nStiffSuccessiveMax <- 15 #maximum step length
  
  if (solver=="dopri" || solver=="dopriStiff") { #complete arguments for dopri

    
    if (sys$argsDopri$hMax>((sys$tEnd-sys$tInit)/(sys$nt-1))) {
      print(sys$argsDopri$hMax)
      warning("hMax is higher than the step length for output, reducing hMax to sys$times[2]-sys$times[1]")
      sys$argsDopri$hMax <- (sys$tEnd-sys$tInit)/(sys$nt)
    }
  }
  
  #DEFINE INPUTS FOR CPP FUNCTION (mostly turn matrices to 1dimensional vectors)
  ntCpp <- sys$nt
  tInitCpp <- sys$tInit
  tEndCpp <- sys$tEnd
  nVCpp <- sys$dim
  nIVCpp <- length(sys$intermediateVar)
  y0Cpp <- sys$y0
  parmsCpp <- sys$parms; names(parmsCpp) <- NULL
  
  formatedExogVar <- formatExogVar(sys)
  samplesExogVarCpp <- formatedExogVar[["samplesExogVarCpp"]]
  nSamplesVarExogVarCpp <-  formatedExogVar[["nSamplesVarExogVarCpp"]]
  nVarExogVarCpp <- formatedExogVar[["nVarExogVarCpp"]]
  argsDopriCpp <- sys$argsDopri
  if(solver=="dopri") {
    if(!exists("dopriForR"))
      stop("No Compiled C++ code found. Use the function cppMakeSys with the argument compile=TRUE to compile.")
    solution <- dopriForR(ntCpp, tInitCpp, tEndCpp, nVCpp, nIVCpp, y0Cpp, parmsCpp, samplesExogVarCpp, nSamplesVarExogVarCpp, nVarExogVarCpp, argsDopriCpp)
  } else if (solver=="RK4Fixed") {
    if(!exists("EulerForR"))
      stop("No Compiled C++ code found. Use the function cppMakeSys with the argument compile=TRUE to compile.")
    solution <- RK4FixedForR(ntCpp, tInitCpp, tEndCpp, nVCpp, nIVCpp, y0Cpp, parmsCpp, samplesExogVarCpp, nSamplesVarExogVarCpp, nVarExogVarCpp)
  } else if (solver=="euler") {
    if(!exists("EulerForR"))
      stop("No Compiled C++ code found. Use the function cppMakeSys with the argument compile=TRUE to compile.")
    solution <- EulerForR(ntCpp, tInitCpp, tEndCpp, nVCpp, nIVCpp, y0Cpp, parmsCpp, samplesExogVarCpp, nSamplesVarExogVarCpp, nVarExogVarCpp)  
  } else if (solver=="dopriStiff") {
    if(!exists("dopriStiffForR"))
      stop("No Compiled C++ code found. Use the function cppMakeSys with the argument compile=TRUE to compile.")
    solution <- dopriStiffForR(ntCpp, tInitCpp, tEndCpp, nVCpp, nIVCpp, y0Cpp, parmsCpp, samplesExogVarCpp, nSamplesVarExogVarCpp, nVarExogVarCpp, argsDopriCpp)
  }
  
  dim <- length(sys$tderiv); dimIv <- length(sys$intermediateVar)
  solution <- matrix(solution,
                     length(sys$times),
                     c(dim, 2*dim, dim+dimIv, 2*dim+dimIv)[sys$returnRK4+1],
                     byrow=TRUE)
  
  #add time vector and name elements
  solution <- cbind(sys$times, solution)
  if(sys$returnRK4==3) {
    colnames(solution) <- c("time", names(sys$tderiv), paste0(names(sys$tderiv),"Dot"), names(sys$intermediateVar))
  } else if(sys$returnRK4==2) {
    colnames(solution) <- c("time", names(sys$tderiv), paste0(names(sys$tderiv),"Dot"))
  } else if(sys$returnRK4==1) {
    colnames(solution) <- c("time", names(sys$tderiv), names(sys$intermediateVar))
  } else {
    colnames(solution) <- c("time",names(sys$tderiv))
  }
  return(solution)
}

###########################################################
## Find ONE zero of a system using Newton-Raphton Method ##
###########################################################
cppNR <- function(sys, parms=NULL, y0=NULL, times=NULL, samplesExogVar=NULL,
                  epsilonJacobian=1e-6, tolNR=1e-6, nIterMaxNR=100, updateY0=FALSE) {
  sys <- completeSys(c("parms", "times", "y0", "samplesExogVar"),
                     sys, envir=environment())
  if(!is.null(sys$samplesExogVar) ) {
    stop("Newton-Raphson algorithm cannot be used on a system with exogenous variables")
  }
  if (is.null(sys$samplesExogVar)) sys$samplesExogVar <- list()
  missingArguments(c("parms", "times", "y0"), sys, envir=environment())
  if (updateY0==TRUE) computeInitialPosition(sys)
  
  #DEFINE INPUTS FOR CPP FUNCTION (mostly turn matrices to 1dimensional vectors)
  ntCpp <- sys$nt
  tInitCpp <- sys$tInit
  tEndCpp <- sys$tEnd
  nVCpp <- length(sys$tderiv)
  nIVCpp <- length(sys$intermediateVar)
  y0Cpp <- sys$y0
  parmsCpp <- sys$parms; names(parmsCpp) <- NULL
  formatedExogVar <- formatExogVar(sys)
  samplesExogVarCpp <- formatedExogVar[["samplesExogVarCpp"]]
  nSamplesVarExogVarCpp <-  formatedExogVar[["nSamplesVarExogVarCpp"]]
  nVarExogVarCpp <- formatedExogVar[["nVarExogVarCpp"]]
  
  
  #Make sure cpp code has been compiled
  if(!exists("NRForR")) #(note: NR stands for Newton-Raphson)
    stop("No Compiled C++ code found. Use the function cppMakeSys with the argument compile=TRUE to compile.")
  
  #Run Cpp code 
  NRForR(ntCpp, tInitCpp, tEndCpp, nVCpp, nIVCpp, y0Cpp, parmsCpp, 
         samplesExogVarCpp, nSamplesVarExogVarCpp, nVarExogVarCpp, 
         epsilonJacobian, tolNR, nIterMaxNR)
}
  
  
##############################
## Run Sensitivity Analysis ##
##############################
cppSA <- function(sys, allParmsSA, fullTrajectory=TRUE, times=NULL, parms=NULL, y0=NULL, solver="dopri", argsDopri=NULL, updateY0=FALSE) {
  sys <- completeSys(c("parms", "times", "y0", "argsDopri"), sys, envir=environment())
  if (is.null(sys$samplesExogVar)) sys$samplesExogVar <- list()
  missingArguments(c("parms", "times", "y0"), sys, envir=environment())
  if (updateY0==TRUE) computeInitialPosition(sys)
  if (is.null(solver) || !(solver %in% c("euler", "RK4Fixed", "dopri", "dopriStiff"))) {
    stop("Unknown solver. Set solver to \"euler\", \"RK4Fixed\" or \"dopri\" ")
  }
  
  #DEFINE INPUTS FOR CPP FUNCTION (mostly turn matrices to 1dimensional vectors)
  ntCpp <- sys$nt
  tInitCpp <- sys$tInit
  tEndCpp <- sys$tEnd
  nVCpp <- length(sys$tderiv)
  nIVCpp <- length(sys$intermediateVar)
  y0Cpp <- sys$y0
  parmsDefaultCpp <- sys$parms; names(parmsDefaultCpp) <- NULL
  formatedExogVar <- formatExogVar(sys)
  samplesExogVarCpp <- formatedExogVar[["samplesExogVarCpp"]]
  nSamplesVarExogVarCpp <-  formatedExogVar[["nSamplesVarExogVarCpp"]]
  nVarExogVarCpp <- formatedExogVar[["nVarExogVarCpp"]]
  
  if (solver=="dopri" || solver=="dopriStiff") { #complete arguments for dopri
    if(is.null(sys$argsDopri)) sys$argsDopri <- list()
    if (is.null(sys$argsDopri$atol)) sys$argsDopri$atol <- 1e-4                          #absolute tolerance for step validation
    if (is.null(sys$argsDopri$rtol)) sys$argsDopri$rtol <- 0                             #relative tolerance for step validation
    if (is.null(sys$argsDopri$fac)) sys$argsDopri$fac <- 0.85                            #reduction factor for step length update
    if (is.null(sys$argsDopri$facMin)) sys$argsDopri$facMin <- 0.1                       #min multiplicative factor for step length update
    if (is.null(sys$argsDopri$facMax)) sys$argsDopri$facMax <- 4                         #max multiplicative factor for step length update
    if (is.null(sys$argsDopri$nStepMax)) sys$argsDopri$nStepMax <- 100                   #maximum number of intermediate steps
    if (is.null(sys$argsDopri$hInit)) sys$argsDopri$hInit <- sys$times[2] - sys$times[1] #initial step length
    if (is.null(sys$argsDopri$hMin)) sys$argsDopri$hMin <- argsDopri$hInit/100           #minimum step length
    if (is.null(sys$argsDopri$hMax)) sys$argsDopri$hMax <- 10*argsDopri$hInit            #maximum step length
    if (is.null(sys$argsDopri$nStiffMax)) sys$argsDopri$nStiffMax <- 100                   #maximum step length
    if (is.null(sys$argsDopri$nStiffSuccessiveMax)) sys$argsDopri$nStiffSuccessiveMax <- 15 #maximum step length
    
    argsDopriCpp <- sys$argsDopri
    
    if (sys$argsDopri$hMax>((sys$tEnd-sys$tInit)/(sys$nt-1))) {
      warning("hMax is higher than the step length for output, reducing hMax to sys$times[2]-sys$times[1]")
      sys$argsDopri$hMax <- (sys$tEnd-sys$tInit)/(sys$nt)
    }
  }
  
  #parms for SA
  nParmsSetSACpp <- nrow(allParmsSA)
  nParmsSACpp <- ncol(allParmsSA)
  parmsPosSACpp <-  sapply(colnames(allParmsSA), function(x) which(names(sys$parms)==x)) - 1
  allParmsSACpp <- c(t(allParmsSA))

  dim <- length(sys$tderiv); dimIv <- length(sys$intermediateVar)
  
  
  if (solver=="dopri") {
    if(!exists("SADopriForR"))
      stop("No Compiled C++ code found. Make sur you compiled the code using initSysCpp")
    
    solution <- SADopriForR(ntCpp, tInitCpp, tEndCpp, nVCpp, nIVCpp, y0Cpp, parmsDefaultCpp, samplesExogVarCpp, nSamplesVarExogVarCpp, 
                            nVarExogVarCpp, argsDopriCpp,
                            nParmsSetSACpp, nParmsSACpp, parmsPosSACpp, allParmsSACpp, fullTrajectory)      
  } else if(solver=="RK4Fixed") {
    if(!exists("SARK4FixedForR"))
      stop("No Compiled C++ code found. Make sur you compiled the code using initSysCpp")
    
    solution <- SARK4FixedForR(ntCpp, tInitCpp, tEndCpp, nVCpp, nIVCpp, y0Cpp, parmsDefaultCpp, samplesExogVarCpp, nSamplesVarExogVarCpp, 
                               nVarExogVarCpp, nParmsSetSACpp, nParmsSACpp, parmsPosSACpp, allParmsSACpp, fullTrajectory)
  } else if (solver=="euler") {
    if(!exists("SAEulerForR"))
      stop("No Compiled C++ code found. Make sur you compiled the code using initSysCpp")
    
    solution <- SAEulerForR(ntCpp, tInitCpp, tEndCpp, nVCpp, nIVCpp, y0Cpp, parmsDefaultCpp, samplesExogVarCpp, nSamplesVarExogVarCpp, 
                            nVarExogVarCpp, nParmsSetSACpp, nParmsSACpp, parmsPosSACpp, allParmsSACpp, fullTrajectory)
  } else if (solver=="dopriStiff") {
    if(!exists("SADopriStiffForR"))
      stop("No Compiled C++ code found. Make sur you compiled the code using initSysCpp")
    
    solution <- SADopriStiffForR(ntCpp, tInitCpp, tEndCpp, nVCpp, nIVCpp, y0Cpp, parmsDefaultCpp, samplesExogVarCpp, nSamplesVarExogVarCpp, 
                            nVarExogVarCpp, argsDopriCpp,
                            nParmsSetSACpp, nParmsSACpp, parmsPosSACpp, allParmsSACpp, fullTrajectory)      
    
  }
  if(fullTrajectory==TRUE) { #COMPUTE FULL TRAJECTORY
    # FORMAT SOLUTION AND RETURN IT
    solution <- array(solution, dim = c(c(dim, 2*dim, dim+dimIv, 2*dim+dimIv)[sys$returnRK4+1], 
                                        ntCpp,
                                        nParmsSetSACpp))
    solution <- aperm(solution, c(3, 2, 1))
    if(sys$returnRK4==3) {
      dimnames(solution)[[3]] <- c(names(sys$tderiv), paste0(names(sys$tderiv),"Dot"), names(sys$intermediateVar))
    } else if(sys$returnRK4==2) {
      dimnames(solution)[[3]] <- c(names(sys$tderiv), paste0(names(sys$tderiv),"Dot"))
    } else if(sys$returnRK4==1) {
      dimnames(solution)[[3]] <- c(names(sys$tderiv), names(sys$intermediateVar))
    } else {
      dimnames(solution)[[3]] <- c(names(sys$tderiv))
    }
    return(solution)
    
  } else {                    # RETURN LAST POINT ONlY
    # FORMAT SOLUTION AND RETURN IT
    solution <- matrix(solution, nParmsSetSACpp, c(dim, 2*dim, dim+dimIv, 2*dim+dimIv)[sys$returnRK4+1], byrow=TRUE)
      
    if(sys$returnRK4==3) {
      dimnames(solution)[[2]] <- c(names(sys$tderiv), paste0(names(sys$tderiv),"Dot"), names(sys$intermediateVar))
    } else if(sys$returnRK4==2) {
      dimnames(solution)[[2]] <- c(names(sys$tderiv), paste0(names(sys$tderiv),"Dot"))
    } else if(sys$returnRK4==1) {
      dimnames(solution)[[2]] <- c(names(sys$tderiv), names(sys$intermediateVar))
    } else {
      dimnames(solution)[[2]] <- c(names(sys$tderiv))
    }
    return(solution)
  }    
}


#################################
## Compute Basin of Attraction ##
#################################
cppBasinBruteForce <- function(sys, 
                               grid, varForConv=NULL,
                               yeq=NULL, tol=0.01, boundsSD=NULL, useDist=TRUE,
                               times=NULL, parms=NULL, y0=NULL, solver="dopri", argsDopri=NULL, updateY0=FALSE) {
  sys <- completeSys(c("parms", "times", "y0", "argsDopri"), sys, envir=environment())
  if (is.null(sys$samplesExogVar)) sys$samplesExogVar <- list()
  missingArguments(c("parms", "times", "y0"), sys, envir=environment())
  if (updateY0==TRUE) computeInitialPosition(sys)  
  
  if (is.null(solver) || !(solver %in% c("euler", "RK4Fixed", "dopri", "dopriStiff"))) {
    stop("Unknown solver. Set solver to \"euler\", \"RK4Fixed\" or \"dopri\" ")
  }
  
  if(is.null(yeq)) yeq <- sys$y0
  if(is.null(boundsSD)) boundsSD <- matrix()
  if(is.null(varForConv)) varForConv <- names(sys$tderiv)
  
  ## MAKE SURE elements in yeq, grid and boundsSD are in the correct order
  allVarNames <- c(names(sys$tderiv), names(sys$intermediateVar))
  yeq <- yeq[match(allVarNames, names(yeq))[which(!is.na(match(allVarNames, names(yeq))))]]
  grid <- grid[, match(allVarNames, colnames(grid))[which(!is.na(match(allVarNames, colnames(grid))))]]
  sys$samplingTimeMinDist <- sys$samplingTimeMinDist[match(allVarNames, names(sys$samplingTimeMinDist))[which(!is.na(match(allVarNames, names(sys$samplingTimeMinDist))))]]
  
  
  #DEFINE INPUTS FOR CPP FUNCTION (mostly turn matrices to 1dimensional vectors)
  ntCpp <- sys$nt
  tInitCpp <- sys$tInit
  tEndCpp <- sys$tEnd
  nVCpp <- length(sys$tderiv)
  nIVCpp <- length(sys$intermediateVar)
  y0Cpp <- sys$y0
  parmsCpp <- sys$parms; names(parmsCpp) <- NULL
  solverCpp <- 0*(solver=="dopri") + 1*(solver=="RK4Fixed") + 2*(solver=="euler") + 3*(solver=="dopriStiff")
  formatedExogVar <- formatExogVar(sys)
  samplesExogVarCpp <- formatedExogVar[["samplesExogVarCpp"]]
  nSamplesVarExogVarCpp <-  formatedExogVar[["nSamplesVarExogVarCpp"]]
  nVarExogVarCpp <- formatedExogVar[["nVarExogVarCpp"]]
  useDistCpp <- useDist
  
  if(is.null(sys$argsDopri)) sys$argsDopri <- list()
  if (is.null(sys$argsDopri$atol)) sys$argsDopri$atol <- 1e-4                          #absolute tolerance for step validation
  if (is.null(sys$argsDopri$rtol)) sys$argsDopri$rtol <- 0                             #relative tolerance for step validation
  if (is.null(sys$argsDopri$fac)) sys$argsDopri$fac <- 0.85                            #reduction factor for step length update
  if (is.null(sys$argsDopri$facMin)) sys$argsDopri$facMin <- 0.1                       #min multiplicative factor for step length update
  if (is.null(sys$argsDopri$facMax)) sys$argsDopri$facMax <- 4                         #max multiplicative factor for step length update
  if (is.null(sys$argsDopri$nStepMax)) sys$argsDopri$nStepMax <- 100                   #maximum number of intermediate steps
  if (is.null(sys$argsDopri$hInit)) sys$argsDopri$hInit <- sys$times[2] - sys$times[1] #initial step length
  if (is.null(sys$argsDopri$hMin)) sys$argsDopri$hMin <- argsDopri$hInit/100           #minimum step length
  if (is.null(sys$argsDopri$hMax)) sys$argsDopri$hMax <- 10*argsDopri$hInit            #maximum step length
  if (is.null(sys$argsDopri$nStiffMax)) sys$argsDopri$nStiffMax <- 100                   #maximum step length
  if (is.null(sys$argsDopri$nStiffSuccessiveMax)) sys$argsDopri$nStiffSuccessiveMax <- 15 #maximum step length
  
  argsDopriCpp <- sys$argsDopri
  
  
  if (solver=="dopri" || solver=="dopriStiff") { #complete arguments for dopri
    
    if (sys$argsDopri$hMax>((sys$tEnd-sys$tInit)/(sys$nt-1))) {
      warning("hMax is higher than the step length for output, reducing hMax to sys$times[2]-sys$times[1]")
      sys$argsDopri$hMax <- (sys$tEnd-sys$tInit)/(sys$nt)
    }
  }
  
  yeqCpp <- yeq;   names(yeqCpp) <- NULL
  boundsSDCpp <- c((boundsSD)); names(boundsSDCpp) <- NULL
  varForConvCpp <- sapply(varForConv, function(x) which(names(sys$tderiv)==x)-1); names(varForConvCpp) <- NULL
  gridSizeCpp <- nrow(grid)
  gridCpp <- c(t(grid))
  
  #Make sure cpp code has been compiled
  if(!exists("basinBruteForceForR"))
    stop("No Compiled C++ code found. Use the function cppMakeSys with the argument compile=TRUE to compile.")
  
  #Simulate model and reformat data as a matrix 
  solution <- basinBruteForceForR(ntCpp, tInitCpp, tEndCpp, 
                                  nVCpp, nIVCpp,
                                  yeqCpp, parmsCpp, 
                                  samplesExogVarCpp, nSamplesVarExogVarCpp, nVarExogVarCpp,
                                  useDistCpp,
                                  tol, boundsSDCpp, varForConvCpp, gridCpp, gridSizeCpp, solverCpp, argsDopriCpp)
}

#############################################################
## Estimate the model empirically using a CMA-ES algorithm ##
#############################################################
cppCMAES <- function(sys, times=NULL, y0 = NULL, samplesExogVar=NULL,
                     parmsInit=NULL, parmsToOptimize=NULL,
                     parmsLower=NULL, parmsUpper=NULL, 
                     dataMinDist, samplingTimeMinDist=NULL, pointsWeightMinDist=NULL, 
                     lambda,
                     sigma=1,
                     nIterMax=200,
                     tol=1e-9,
                     standardizeParms=TRUE, 
                     solver="dopri", 
                     nReDrawMax = 100, 
                     argsDopri=NULL, updateY0=FALSE) {
  sys <- completeSys(c("parmsInit", "parmsToOptimize", "parmsLower", "parmsUpper", "times", "y0", "samplesExogVar", "argsDopri",
                       "dataMinDist", "samplingTimeMinDist", "pointsWeightMinDist", "nReDrawMax"), sys, envir=environment())

  ## Handling initial positions as inputs to optimize (stuffed in parmsInit for now for simplicity)
  #check uniqueness of names in parmsInit
  if (length(names(sys$parmsInit))!=length(unique(names(sys$parmsInit)))) {
    errorMsg <- paste0("At least one parameter in parmsInit is defined multiple times.\\ Parameter(s): ",
                       names(sys$parmsInit)[which(sapply(names(sys$parmsInit),
                                                         function(x) length(which(names(sys$parmsInit)==x)))>1)][1]
                       )
    stop(errorMsg)
  }
  
  if (is.null(solver) || !(solver %in% c("euler", "RK4Fixed", "dopri", "dopriStiff"))) {
    stop("Unknown solver. Set solver to \"euler\", \"RK4Fixed\" or \"dopri\" ")
  }
  
  for (i in which(sys$parmsToOptimize)) {
    if (parmsInit[i]<parmsLower[i] || parmsInit[i]>parmsUpper[i]) {
      stop(paste0("Logical error in CMAES initialization: initial value of parameter \"", names(sys$parms)[i], "\" is out of the domain to explore."))
    }
  }
  
  if (is.null(sys$samplesExogVar)) sys$samplesExogVar <- list()
  
  if (is.null(sys$parmsToOptimize)) sys$parmsToOptimize <- rep(TRUE, length(sys$parmsInit))
  if (is.null(sys$samplingTimeMinDist)) sys$samplingTimeMinDist <- lapply(sys$dataMinDist, function(x) seq(sys$times[1], sys$times[length(sys$times)], length.out = length(x)))
  if (is.null(sys$pointsWeightMinDist)) sys$pointsWeightMinDist <- lapply(sys$dataMinDist, function(x) rep(1, length(x)))
  
  if(!is.list(sys$dataMinDist)) stop("Incorrect data structure: dataMinDist must be a list of (named) vectors")
  if(!is.list(sys$samplingTimeMinDist)) stop("Incorrect data structure: samplingTimeMinDist must be a list of (named) vectors")
  if(!is.list(sys$pointsWeightMinDist)) stop("Incorrect data structure: pointsWeightMinDist must be a list of (named) vectors")
  
  if (is.null(sys$parmsLower)) {
    sys$parmsLower <- rep(-Inf, length(sys$parms))
  }
  if (is.null(sys$parmsUpper)) {
    sys$parmsUpper <- rep(+Inf, length(sys$parms))
  }
  

  if(is.null(sys$argsDopri)) sys$argsDopri <- list()
  if (is.null(sys$argsDopri$atol)) sys$argsDopri$atol <- 1e-4                          #absolute tolerance for step validation
  if (is.null(sys$argsDopri$rtol)) sys$argsDopri$rtol <- 0                             #relative tolerance for step validation
  if (is.null(sys$argsDopri$fac)) sys$argsDopri$fac <- 0.85                            #reduction factor for step length update
  if (is.null(sys$argsDopri$facMin)) sys$argsDopri$facMin <- 0.1                       #min multiplicative factor for step length update
  if (is.null(sys$argsDopri$facMax)) sys$argsDopri$facMax <- 4                         #max multiplicative factor for step length update
  if (is.null(sys$argsDopri$nStepMax)) sys$argsDopri$nStepMax <- 100                   #maximum number of intermediate steps
  if (is.null(sys$argsDopri$hInit)) sys$argsDopri$hInit <- sys$times[2] - sys$times[1] #initial step length
  if (is.null(sys$argsDopri$hMin)) sys$argsDopri$hMin <- argsDopri$hInit/100           #minimum step length
  if (is.null(sys$argsDopri$hMax)) sys$argsDopri$hMax <- 10*argsDopri$hInit            #maximum step length
  if (is.null(sys$argsDopri$nStiffMax)) sys$argsDopri$nStiffMax <- 100                   #maximum step length
  if (is.null(sys$argsDopri$nStiffSuccessiveMax)) sys$argsDopri$nStiffSuccessiveMax <- 15 #maximum step length
  
    
  if (solver=="dopri" || solver=="dopriStiff") { #complete arguments for dopri    
    if (sys$argsDopri$hMax>((sys$tEnd-sys$tInit)/(sys$nt-1))) {
      warning("hMax is higher than the step length for output, reducing hMax to sys$times[2]-sys$times[1]")
      sys$argsDopri$hMax <- (sys$tEnd-sys$tInit)/(sys$nt)
    }
  }
  
  missingArguments(c("parms", "times", "y0", "parmsInit", "dataMinDist"), sys, envir=environment())
  if (updateY0==TRUE) computeInitialPosition(sys)
  
  # check all y0 are in parms, are correctly ordered, and complete/correct if not
  tempParmsInit <- c()
  tempParmsToOptimize <- c()
  tempParmsLower <- c()
  tempParmsUpper <- c()
  swapEltAndName <- function(vec, i, j) { #temporary function used only here. Performs swap of elements (and their names) in a vector. Also work on unnamed vectors.
    memName <- names(vec)[i]; memValue <- vec[i]
    vec[i] <- vec[j]; names(vec)[i] <- names(vec)[j]
    vec[j] <- memValue; names(vec)[j] <- memName
    vec
  }
  for (i in 1:length(sys$y0)) {
    if (names(sys$y0)[i]!=names(sys$parmsInit)[i]) {
      if (names(sys$y0)[i] %in% names(sys$parmsInit)) { #element is in parmsInit, but not at correct position, put it at correct position
        indexWrong <- which(names(sys$parmsInit)==names(sys$y0Text)[i])
        tempParmsInit       <- swapEltAndName(sys$parmsInit, i, indexWrong)
        tempParmsToOptimize <- swapEltAndName(sys$parmsToOptimize, i, indexWrong)
        tempParmsLower      <- swapEltAndName(sys$parmsLower, i, indexWrong)
        tempParmsUpper      <- swapEltAndName(sys$parmsUpper, i, indexWrong)
      } else { #element missing in parmsInit, add it but keep it fixed
        tempParmsInit       <- c(tempParmsInit, sys$y0[i])
        tempParmsToOptimize <- c(tempParmsToOptimize, FALSE)
        tempParmsLower      <- c(tempParmsLower, -Inf)
        tempParmsUpper      <- c(tempParmsUpper, +Inf)
      }
    } else { #element is here and at correct position
      tempParmsInit       <- c(tempParmsInit, sys$parmsInit[i])
      tempParmsToOptimize <- c(tempParmsToOptimize, sys$parmsToOptimize[i])
      tempParmsLower      <- c(tempParmsLower, sys$parmsLower[i])
      tempParmsUpper      <- c(tempParmsUpper, sys$parmsUpper[i])
    }
  }
  tempParmsInit       <- c(tempParmsInit, sys$parmsInit[!names(sys$parmsInit)%in% names(sys$y0)])
  tempParmsToOptimize <- c(tempParmsToOptimize, sys$parmsToOptimize[!names(sys$parmsInit)%in% names(sys$y0)])
  tempParmsLower      <- c(tempParmsLower, sys$parmsLower[!names(sys$parmsInit)%in% names(sys$y0)])
  tempParmsUpper      <- c(tempParmsUpper, sys$parmsUpper[!names(sys$parmsInit)%in% names(sys$y0)])
  
  sys$parmsInit       <- tempParmsInit
  sys$parmsToOptimize <- tempParmsToOptimize
  sys$parmsLower      <- tempParmsLower
  sys$parmsUpper      <- tempParmsUpper
  
  
  ## Make sure data for minDist are in the correct order
  allVarNames <- c(names(sys$tderiv), names(sys$intermediateVar))
  sys$dataMinDist <- sys$dataMinDist[match(allVarNames, names(sys$dataMinDist))[which(!is.na(match(allVarNames, names(sys$dataMinDist))))]]
  sys$samplingTimeMinDist <- sys$samplingTimeMinDist[match(allVarNames, names(sys$samplingTimeMinDist))[which(!is.na(match(allVarNames, names(sys$samplingTimeMinDist))))]]
  
  void <- lapply(sys$samplingTimeMinDist, 
                 function(i) sapply(i,
                                    function(t) {
                                      if(t<min(sys$times) || t>max(sys$times)) {
                                        stop("Invalid samplingTimeMinDist: sampling time is out of time interval")
                                      }
                                    }))
  
  #DEFINE INPUTS FOR CPP FUNCTION (mostly turn matrices to 1dimensional vectors)
  ntCpp <- sys$nt
  tInitCpp <- sys$tInit
  tEndCpp <- sys$tEnd
  nVCpp <- length(sys$tderiv)
  nIVCpp <- length(sys$intermediateVar)
  parmsCpp <- sys$parmsInit; names(parmsCpp) <- NULL
  formatedExogVar <- formatExogVar(sys)
  samplesExogVarCpp <- formatedExogVar[["samplesExogVarCpp"]]
  nSamplesVarExogVarCpp <-  formatedExogVar[["nSamplesVarExogVarCpp"]]
  nVarExogVarCpp <- formatedExogVar[["nVarExogVarCpp"]]
  nVarMinDistCpp <- length(sys$dataMinDist)
  varMinDistCpp <- which(sapply(allVarNames, function(x) x %in% names(sys$dataMinDist)))-1; names(varMinDistCpp) <- NULL
  dataMinDistCpp <- unlist(sys$dataMinDist); names(dataMinDistCpp) <- NULL
  samplingTimeMinDistCpp <- unlist(sys$samplingTimeMinDist); names(samplingTimeMinDistCpp) <- NULL
  pointsWeightMinDistCpp <- unlist(sys$pointsWeightMinDist); names(pointsWeightMinDistCpp) <- NULL
  nObsVarMinDistCpp <- sapply(sys$dataMinDist, length); names(nObsVarMinDistCpp) <- NULL
  parmsToOptimizeCpp <- sys$parmsToOptimize
  parmsLowerCpp <- sys$parmsLower[sys$parmsToOptimize]; names(parmsLowerCpp) <- NULL
  parmsUpperCpp <- sys$parmsUpper[sys$parmsToOptimize]; names(parmsUpperCpp) <- NULL
  nReDrawMaxCpp <- sys$nReDrawMax
  argsDopriCpp <- sys$argsDopri
  solverCpp <- 0*(solver=="dopri") + 1*(solver=="RK4Fixed") + 2*(solver=="euler") + 3*(solver=="dopriStiff")
  
  #Make sure cpp code has been compiled
  if(!exists("CMAESForR"))
    stop("No Compiled C++ code found. Use the function cppMakeSys with the argument compile=TRUE to compile.")
  
  #Simulate model and reformat data as a matrix 
  out <- CMAESForR(ntCpp, tInitCpp, tEndCpp, nVCpp, nIVCpp, parmsCpp,
                   samplesExogVarCpp, nSamplesVarExogVarCpp, nVarExogVarCpp, 
                   nVarMinDistCpp, varMinDistCpp, dataMinDistCpp, samplingTimeMinDistCpp, 
                   pointsWeightMinDistCpp, nObsVarMinDistCpp, parmsToOptimizeCpp, 
                   parmsLowerCpp, parmsUpperCpp, 
                   lambda, sigma, nIterMax, tol, TRUE, nReDrawMaxCpp, solverCpp, argsDopriCpp)
  tempParms <- parmsInit
  tempParms[parmsToOptimize] <- out$parms
  out$parms <- c(tempParms)
  out
}

####################################################################
## Run Hybrid Cell Mapping (HCM) with Subdivision                 ##
## to find ALL the zeroes of the system of differential equations ##
## (Uses Newton-Raphson algorithm)                                ##
####################################################################
cppSCMZeroes <- function(sys, times=NULL, parms=NULL,
                         samplesExogVar=NULL,
                         nCells, yLower, yUpper, 
                         epsilonJacobian=1e-6, tolNR=1e-6, nIterMaxNR=20,
                         postProcessing=TRUE) {
  sys <- completeSys(c("parms", "times", "nCells", "yLower", "yUpper", 
                       "epsilonJacobian", "tolNR", "nIterMaxNR", "samplesExogVar"),
                     sys,
                     envir=environment())
  if(!is.null(sys$samplesExogVar)) {
    stop("Cell-Mapping methods cannot be used on a system with exogenous variables")
  }
  if (is.null(sys$samplesExogVar)) sys$samplesExogVar <- list()
  
  missingArguments(c("parms", "nCells", "yLower", "yUpper", 
                     "epsilonJacobian", "tolNR", "nIterMaxNR"),
                   sys, envir=environment())
  
  #DEFINE INPUTS FOR CPP FUNCTION (mostly turn matrices to 1dimensional vectors)
  ntCpp <- sys$nt
  tInitCpp <- sys$tInit
  tEndCpp <- sys$tEnd
  nVCpp <- length(sys$tderiv)
  nIVCpp <-length(sys$intermediateVar)
  parmsCpp <- sys$parms; names(parmsCpp) <- NULL
  formatedExogVar <- formatExogVar(sys)
  samplesExogVarCpp <- formatedExogVar[["samplesExogVarCpp"]]
  nSamplesVarExogVarCpp <-  formatedExogVar[["nSamplesVarExogVarCpp"]]
  nVarExogVarCpp <- formatedExogVar[["nVarExogVarCpp"]]
  
  nCellsCpp <- sys$nCells; names(nCellsCpp) <- NULL
  yLowerCpp <- sys$yLower; names(yLowerCpp) <- NULL
  yUpperCpp <- sys$yUpper; names(yUpperCpp) <- NULL
  byCellsCpp <- (yUpperCpp - yLowerCpp)/nCellsCpp
  epsilonJacobianCpp <- sys$epsilonJacobian
  tolNRCpp <- sys$tolNR
  nIterMaxCpp <- sys$nIterMaxNR
  postProcessingCpp <- postProcessing
  #Make sure cpp code has been compiled
  if(!exists("SCMNR"))
    stop("No Compiled C++ code found. Use the function cppMakeSys with the argument compile=TRUE to compile.")
  
  #Simulate model and reformat data as a matrix 
  solution <- SCMNR(ntCpp, tInitCpp, tEndCpp, nVCpp, nIVCpp,
                    parmsCpp,
                    samplesExogVarCpp, nSamplesVarExogVarCpp, nVarExogVarCpp, 
                    nCellsCpp, yLowerCpp, yUpperCpp, 
                    byCellsCpp, epsilonJacobianCpp, 
                    tolNRCpp, nIterMaxCpp, postProcessingCpp)
  
}

####################################################################
## Run Hybrid Cell Mapping (HCM) with Subdivision                 ##
## to find ALL the zeroes of the system of differential equations ##
## (Uses Newton-Raphson algorithm)                                ##
####################################################################
cppFindZeroes <- function(sys, times=NULL, parms=NULL, samplesExogVar=NULL,
                          nCells, yLower, yUpper, 
                          epsilonJacobian=1e-6, tolNR=1e-6, nIterMaxNR=20,
                          nSubdivisionMax, nSubdivisionGCM=nSubdivisionMax-1,
                          gridCell, postProcessing=TRUE) {
  sys <- completeSys(c("parms", "times", "nCells", "yLower", "yUpper", 
                       "epsilonJacobian", "tolNR", "nIterMaxNR", 
                       "nSubdivisionMax", "nSubdivisionGCM", 
                       "gridCell", "samplesExogVar"),
                     sys,
                     envir=environment())
  if(!is.null(sys$samplesExogVar)) {
    stop("Cell-Mapping methods cannot be used on a system with exogenous variables")
  }
  if (is.null(sys$samplesExogVar)) sys$samplesExogVar <- list()
  
  missingArguments(c("parms", "nCells", "yLower", "yUpper", 
                     "epsilonJacobian", "tolNR", "nIterMaxNR", 
                     "nSubdivisionMax", "nSubdivisionGCM", 
                     "gridCell"), sys, envir=environment())
  
  #DEFINE INPUTS FOR CPP FUNCTION (mostly turn matrices to 1dimensional vectors)
  ntCpp <- sys$nt
  tInitCpp <- sys$tInit
  tEndCpp <- sys$tEnd
  nVCpp <- length(sys$tderiv)
  nIVCpp <- length(sys$intermediateVar)
  parmsCpp <- sys$parms; names(parmsCpp) <- NULL
  formatedExogVar <- formatExogVar(sys)
  samplesExogVarCpp <- formatedExogVar[["samplesExogVarCpp"]]
  nSamplesVarExogVarCpp <-  formatedExogVar[["nSamplesVarExogVarCpp"]]
  nVarExogVarCpp <- formatedExogVar[["nVarExogVarCpp"]]
  
  nCellsCpp <- sys$nCells; names(nCellsCpp) <- NULL
  yLowerCpp <- sys$yLower; names(yLowerCpp) <- NULL
  yUpperCpp <- sys$yUpper; names(yUpperCpp) <- NULL
  byCellsCpp <- (yUpperCpp - yLowerCpp)/nCellsCpp
  epsilonJacobianCpp <- sys$epsilonJacobian
  tolNRCpp <- sys$tolNR
  nIterMaxCpp <- sys$nIterMaxNR
  nSubdivisionMaxCpp <- sys$nSubdivisionMax
  nSubdivisionGCMCpp <- sys$nSubdivisionGCM
  gridCellCpp <- c(t(sys$gridCell))
  nEltsGridCellCpp <- nrow(sys$gridCell)
  postProcessingCpp <- postProcessing
  #Make sure cpp code has been compiled
  if(!exists("HCMSubdivisionNR"))
    stop("No Compiled C++ code found. Use the function cppMakeSys with the argument compile=TRUE to compile.")
  
  #Simulate model and reformat data as a matrix 
  solution <- HCMSubdivisionNR(ntCpp, tInitCpp, tEndCpp, nVCpp, nIVCpp, parmsCpp,
                               samplesExogVarCpp,
                               nSamplesVarExogVarCpp, nVarExogVarCpp, 
                               nCellsCpp, yLowerCpp, yUpperCpp, 
                               byCellsCpp, epsilonJacobianCpp, 
                               tolNRCpp, nIterMaxCpp, nSubdivisionMaxCpp, 
                               nSubdivisionGCMCpp, gridCellCpp, nEltsGridCellCpp, postProcessingCpp)
  
}