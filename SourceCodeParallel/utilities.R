list.of.packages <- c("Rcpp", "RcppArmadillo")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
void <- sapply(list.of.packages, function(x) library(x, character.only=T))

###########################################################################
##                           Internal function                           ##
## Load the model from an external file and create an object sys with it ##
###########################################################################
loadModel <- function(fileName, samplesExogVar=NULL) {
  options(warn=-1)
  modelFile <- file(fileName)
  modelText <- readLines(modelFile, n = -1)
  close(modelFile)
  out <- list()
  ## remove blank lines
  emptyLines <- which(modelText=="")
  if(length(emptyLines)!=0) modelText <- modelText[-emptyLines]
  ## remove simple comment lines (i.e. those that begin with a "#", but not with a "##")
  modelText <- modelText[sapply(modelText,      
                                function(line) ifelse(substr(line, 1, 1)=="#" && substr(line, 1, 2) !="##", FALSE, TRUE)
                                )]
  ## EXOGENOUS VARIABLES ##
  beginExogenousVar <- which(grepl("##exogenous variables", modelText))[1]
  if (!is.na(beginExogenousVar)) {
    endExogenousVar <- min(which(regexpr("##", modelText[-(1:beginExogenousVar)])==1)) # search next row beginning by a "##"
    
    if (length(endExogenousVar)>0 && endExogenousVar>1) {
      ExogenousVarText <- modelText[beginExogenousVar+1:(endExogenousVar-1)]
      out$samplesExogVar <- lapply(createEquations(ExogenousVarText), function(x) eval(parse(text=x)))
      out$intermediateVar <- createExogVar(names(out$samplesExogVar))
    }
  } else {
  	if(!is.null(samplesExogVar)) {
  		ExogenousVarText <- sapply(1:length(samplesExogVar), 
  															 function(i) {
  															 	paste0(names(samplesExogVar)[i], " = c(", paste0(as.character(samplesExogVar[[i]]), collapse=", "), ")")
  															 })
  		out$samplesExogVar <- lapply(createEquations(ExogenousVarText), function(x) eval(parse(text=x)))
  		out$intermediateVar <- createExogVar(names(out$samplesExogVar))
  	}
  }
  
  ## INTERMEDIATE VARIABLES ##
  beginIntermediateVar <- which(grepl("##intermediate variables", modelText))[1]
  if (!is.na(beginIntermediateVar)) {
    endIntermediateVar <- min(which(regexpr("##", modelText[-(1:beginIntermediateVar)])==1), length(modelText)-beginIntermediateVar+1) # search next row beginning by a "##"
    
    if (length(endIntermediateVar)>0 && endIntermediateVar>1) {
      intermediateVarText <- modelText[beginIntermediateVar+1:(endIntermediateVar-1)]
      out$intermediateVar <- c(out$intermediateVar, createEquations(intermediateVarText))
    }
  }
  
  ## TIMES DERIVATIVES ##
  beginTderiv <- which(grepl("##time derivatives", modelText))[1]
  # if (is.na(beginTderiv)) {
  #   stop("Missing time derivatives")
  # }
  if (!is.na(beginTderiv)) {
    endTderiv <- min(c(which(regexpr("##", modelText[-(1:beginTderiv)])==1), length(modelText)-beginTderiv+1)) # search next row beginning by a "#"
    if (length(endTderiv)!=0 && endTderiv>1) {
      tderivText <- modelText[beginTderiv+1:(endTderiv-1)]
      out$tderiv <- createEquations(tderivText)
    }
  }
  
  ## PARAMETERS ##
  beginParms <- which(grepl("##parameters", modelText))[1]
  if(!is.na(beginParms)) {
    endParms <- min(c(which(regexpr("##", modelText[-(1:beginParms)])==1), length(modelText)-beginParms+1)) # search next row beginning by a "#"
    if (length(endParms)!=0 && endParms>1) {
      parmsText <- modelText[beginParms+1:(endParms-1)]
      out$parms <- sapply(createEquations(parmsText), function(x) eval(parse(text=x)))
      sapply(1:length(out$parms), function(i) eval(call("=", names(out$parms)[i], out$parms[i]), envir=parent.env(environment())))
    }
  }
  
  ## INITIAL VALUES ##
  beginInitialValues <- which(grepl("##initial values", modelText))[1]
  if(!is.na(beginInitialValues)) {
    endInitialValues <- min(c(which(regexpr("##", modelText[-(1:beginInitialValues)])==1), length(modelText)-beginInitialValues+1)) # search next row beginning by a "#"
    if (length(endInitialValues)!=0 && endInitialValues>1) {
      y0Text <- modelText[beginInitialValues+1:(endInitialValues-1)]
      out$y0Text <- sapply(createEquations(y0Text), function(x) eval(parse(text=x)))
    }
  }
  
  ## TIME SEQUENCE ##
  beginTime <- which(modelText=="##time")[1]
  if(!is.na(beginTime)) {
    endTime <- min(c(which(regexpr("##", modelText[-(1:beginTime)])==1), length(modelText)-beginTime+1)) # search next row beginning by a "#"
    if (length(endTime)!=0 && endTime>1) {
      timeText <- modelText[beginTime+1:(endTime-1)]
      out$times <- sapply(createEquations(timeText), function(x) eval(parse(text=x)))
      out$times <- seq(out$times["begin"], out$times["end"], out$times["by"])
    }
  }
  
  ## EVENT TIME ##
  beginEventTime <- which(grepl("##events time", modelText))[1]
  if(!is.na(beginEventTime)) {
    endEventTime <- min(c(which(regexpr("##", modelText[-(1:beginEventTime)])==1), length(modelText)-beginEventTime+1)) # search next row beginning by a "#"
    if (length(endEventTime)!=0 && endEventTime>1) {
      out$eventTime <- lapply(1:(endEventTime-1), function(i) eval(parse(text=modelText[beginEventTime+1:(endEventTime-1)][i])[[1]]))
    }
  }
  
  ## EVENT VARIABLES ##
  beginEventVar <- which(grepl("##events variable", modelText))[1]
  if(!is.na(beginEventVar)) {
    endEventVar <- min(c(which(regexpr("##", modelText[-(1:beginEventVar)])==1), length(modelText)-beginEventVar+1)) # search next row beginning by a "#"
    if (length(endEventVar)!=0 && endEventVar>1) {
      out$eventVar <- lapply(1:(endEventVar-1), function(i) eval(parse(text=modelText[beginEventVar+1:(endEventVar-1)][i])[[1]]))
    }
  }
  
  ## COMPLETE TDERIV IF NEEDED ##
    for (i in 1:length(out$intermediateVar)) {
      if (length(grep("Dot", names(out$intermediateVar)[i]))>0) {
        out$tderiv <- c(out$tderiv, out$intermediateVar[i])
        names(out$tderiv)[length(names(out$tderiv))] <- substr(names(out$intermediateVar)[i], 1, nchar(names(out$intermediateVar)[i])-3)
      }
    }
  # Makes sure variable in y0Text and tderiv are in the same order
  out$y0Text <- reorderElements(toReorder = out$y0Text, orderedNames = names(out$tderiv), what = "y0Text")  
  #Remove variables from intermediate variables
  removeVariablesFromIntermediateVar <- rep(FALSE, length(out$intermediateVar))
  for (yNames in names(out$y0Text)) {
    removeVariablesFromIntermediateVar[which(names(out$intermediateVar)==paste0(yNames, "Dot"))] <- TRUE 
  }
  out$intermediateVar <- out$intermediateVar[!removeVariablesFromIntermediateVar]
  options(warn=0)
  out
}
##########################################################################
##                           Internal function                          ##
## Used by load model to turn str in equations for sys                  ##
##########################################################################
## I included (pre) formatting of data here for simplicity. I might move it somewhere else later.
## the problem is that it makes some of the formatting visible in the R object sys(in sys$intermediateVar)...
createExogVar <- function(strExogVarNames) {
  out <- c()
  if (length(strExogVarNames)>0) {
    for (i in 1:length(strExogVarNames)) {
      lineText <- paste0("thisMakeThisARightArrowgetExogVar(t, ", as.character(i-1), ");") #the first "this" is the this of C++ (this->x), not a this from formatting in R !!!
      out <- c(out, lineText)
      names(out)[length(out)] <- strExogVarNames[i]
    }
  }
  out
}
##########################################################################
##                           Internal function                          ##
## Used by load model to turn str in equations for sys                  ##
##########################################################################
createEquations <- function(strEquations) {
  out <- c()
  for (i in 1:length(strEquations)) {
    lineText <- strEquations[i]
    lineText <- sub("#.*", "", lineText) ## Remove any trailing comment
    if (!grepl("^[[:blank:]]*$", lineText)) { ## Skip empty lines or lines containing only spaces
      tempStr1Line = strsplit(lineText, "=")[[1]]
      
      #This is to manage the case when there are other = in the equation (use of logical operators)
      if(length(tempStr1Line)>2){
        tempStr <- tempStr1Line[2]
        for(iter in 3:length(tempStr1Line)){
          tempStr <- paste(tempStr,tempStr1Line[3],sep="=")
          tempStr1Line <- tempStr1Line[-3]
        }
        tempStr1Line[2] <- tempStr
      }
      #Replacing all reseverved words, for now only in-> inv
      
      tempStr1Line[1] <- gsub("[[:space:]]","",tempStr1Line[1])
      tempStr1Line[2] <- gsub("[[:space:]]","",tempStr1Line[2])
      out <- c(out, tempStr1Line[2])
      names(out)[length(out)] <- tempStr1Line[1]
    }
  }
  out
}
################################################################
##                    Internal function                       ##
## Used to initialize sys at the beginning of a function,     ##
## Takes a vector of str allNames and an environment envir    ##
## as input, and append the elements name in allNames to sys  ##
## creates the sys object from scratch if none is provoided   ##
################################################################
completeSys <- function(allNames, sys=NULL, envir=environment()) {
  if (is.null(sys))
    sys <- list()
  for (i in 1:length(allNames)) { #for all elements to complete
    if(!is.null(eval(as.name(allNames[i]), envir=envir))) { #if non default value specified
      if (identical(eval(as.name(allNames[i]), envir=envir), FALSE)) { #if non default value is FALSE
        sys[allNames[[i]]] <- list(NULL) # set element to NULL
      } else { # if non default value is anything other than FALSE
        if (allNames[i]=="parms") {
          if (is.null(sys$parms)) {
            sys$parms <- eval(eval(as.name(allNames[i]), envir=envir), envir=envir)
          } else {
            for (j in 1:length(eval(as.name(allNames[i]), envir=envir))) sys$parms[names(eval(as.name(allNames[i]), envir=envir))[j]] <- eval(as.name(allNames[i]), envir=envir)[j]
          }
        } else if (allNames[i]=="y0") {
          if (is.null(sys$y0)) {
            y0Temp <- eval(as.name(allNames[i]), envir=envir)
            #reordering:
            sys$y0 <- sapply(1:length(sys$tderiv), function(j) y0Temp[which(names(y0Temp)==names(sys$tderiv)[j])])
          } else {
            for (j in 1:length(eval(as.name(allNames[i]), envir=envir))) sys$y0[names(eval(as.name(allNames[i]), envir=envir))[j]] <- eval(as.name(allNames[i]), envir=envir)[j]
          }
        } else if (allNames[i]=="y0Text") {
          if (is.null(sys$y0Text)) {
            y0TextTemp <- eval(as.name(allNames[i]), envir=envir)
            #reordering:
            sys$y0Text <- unlist(sapply(1:length(sys$tderiv), function(j) y0TextTemp[which(names(y0TextTemp)==names(sys$tderiv)[j])]))
          } else {
            for (j in 1:length(eval(as.name(allNames[i]), envir=envir))) {
              sys$y0Text[names(eval(as.name(allNames[i]), envir=envir))[j]] <- eval(as.name(allNames[i]), envir=envir)[j]
            }
          }
        } else if (allNames[i]=="argsDopri") {
          if (is.null(sys$argsDopri)) {
            sys$argsDopri <- eval(as.name(allNames[i]), envir=envir)
          } else {
            for (j in 1:length(eval(as.name(allNames[i]), envir=envir))) sys$argsDopri[[names(eval(as.name(allNames[i]), envir=envir))[j]]] <- eval(as.name(allNames[i]), envir=envir)[[j]]
          }
        } else {
          sys[[allNames[i]]] <- eval(as.name(allNames[i]), envir=envir) #replace element in sys by non-default value
        }
      }
    }
  }
  if (any(allNames=="times")) { #specific treatment for times
    sys$tInit <- sys$times[1]
    sys$tEnd <- sys$times[length(sys$times)]
    sys$nt <- length(sys$times)
  }
  sys
}
################################################################
##                    Internal function                       ##
## Test if there is any missing input when calling a function ##
## call missingArgument at the beginning of any function      ##
## using a sys object, to check that sys contains everything  ##
## the function needs to run.                                 ##
## takes a vector of str allNames as input, which defines the ##
## needed elements                                            ##
################################################################
missingArguments <- function(allNames, sys=NULL, envir=environment()) {
  if (is.null(sys)) {
    missingArguments <-(sapply(allNames, function(x) is.null(eval(as.name(x),
                                                                  envir=envir))))
    if (any(missingArguments)) {
      stop("Missing arguments: ", paste(allNames[missingArguments], collapse=", "))
    }
  } else {
    missingArguments <-(sapply(allNames, function(x) is.null(eval(call("$", quote(sys), as.name(x)),
                                                                  envir=envir))))
    if (any(missingArguments)) {
      stop("Missing elements in sys: ", paste(paste0("sys$", allNames[missingArguments]), collapse=", "))
    }
  }
}
################################################################################
##                               Internal function                            ##
## Performs basic tests to ensure the structure of the sys object is correct  ##
##  WARNING: will only verify SYSTEM STRUCTURE (i.e. missing/extra elements)  ##
##  BUT does not verify LOGICAL CORRECTNESS (i.e. absence of  circularities)  ##
################################################################################

checkModeDefinition <- function(sys) {
  ## Test correct initial variables value (if any)
  searchMissingVariableInitialValue(sys)
  searchMissingVariableDefinition(sys)
    
  if (is.null(sys$intermediateVar)) {
    ##  Test if there is one or more variable that is defined multiple times
    searchDuplicateElements(allEltsNames=names(sys$tderiv),
                            what="variable", problem = "is defined multiple times")
  } else {
  ##  Test if there is one or more variable or intermediate variable that is defined multiple times 
  searchDuplicateElements(allEltsNames=c(names(sys$tderiv), names(sys$intermediateVar)),
                          what="variable or intermediate variable", problem="is defined multiple times")
  }
  ## Test if a parameter and a variable have the same name
  if (!is.null(sys$parms)) searchDuplicateElements(allEltsNames=c(names(sys$parms), names(sys$tderiv)),
                                                   what="element", problem="is defined as both a variable and a parameter")
  ## Test if a parameter and an intermediate variable have the same name
  if (!is.null(sys$intermediateVar) && !is.null(sys$parms)) searchDuplicateElements(allEltsNames=c(names(sys$parms), names(sys$intermediateVar)),
                                                                                    what="element", problem="is defined as both an intermediate variable and a parameter")
  
  if (!is.null(sys$y0Text)) searchDuplicateElements(allEltsNames=names(sys$y0Text),
                                                    what="initial value", problem="is defined multiple times")
  
  ## Make sure there are no undefined elements (ie variable, intermediate variable or parameter used in model equations but not defined)
  # NOT IMPLEMENTED YET (difficult to do in practice...)
}

##############################################################################
##                             Internal function                            ##
##   Search for multiple occurences of elements in a vector of elements     ##
##   and build appropriate error message                                    ## 
##############################################################################
searchDuplicateElements <- function(allEltsNames, what = "element", problem="is defined multiple times") {
  uniqueEltsNames <- unique(allEltsNames)
  nOccurUniqueEltsNames <- sapply(uniqueEltsNames, function(name) length(which(allEltsNames ==name)))
  nDuplicatedEltsNames <- length(which(nOccurUniqueEltsNames>1))
  duplicatedEltsNames <- names(nOccurUniqueEltsNames)[which(nOccurUniqueEltsNames>1)]
  if(nDuplicatedEltsNames>0) {
    errorMsg <- paste0("Incorrect model structure: at least one ", what, " ", problem, ". \n Element(s) concerned: ", collapse="")
    errorMsg <- paste0(errorMsg, paste(duplicatedEltsNames, collapse=", "))
    stop(errorMsg)
  }
}
##############################################################################
##                             Internal functions                           ##
##   Search for incorrect initial variables values definition               ##
##   i.e. either a variable defined in tderiv but not in y0Text             ##
##   or a variable defined in y0Text but not in tderiv                      ##
##############################################################################
searchMissingVariableInitialValue <- function(sys) {
  if (!is.null(sys$y0Text)) { # Make sure all variables have both a definition in tderiv and an initial position in y0Text
    missingVarsTderiv <- names(sys$tderiv)[-which(names(sys$tderiv) %in% names(sys$y0Text))]
    if (length(missingVarsTderiv)>0) { #more time derivatives than initial values
      errorMsg <- "Incorrect model structure: at least one variable is defined in time derivatives but has no initial value. \n Variable(s) concerned: "
      errorMsg <- paste0(errorMsg, paste(missingVarsTderiv, collapse =", "))
      stop(errorMsg)
    }
  }
}

searchMissingVariableDefinition <- function(sys) {
  if (!is.null(sys$y0Text)) { # Make sure all variables have both a definition in tderiv and an initial position in y0Text
    missingVarsY0Text <- names(sys$y0Text)[-which(names(sys$y0Text) %in% names(sys$tderiv))]
    if (length(missingVarsY0Text)>0) { #more initial values than time derivatives
      errorMsg <- "Incorrect model structure: at least one variable has an initial value but is not defined in time derivatives. \n Variable(s) concerned: "
      errorMsg <- paste0(errorMsg, paste(missingVarsY0Text, collapse =", "))
      stop(errorMsg)
    }
  }
}

##############################################################################
##                             Internal function                            ##
##   Reorder elements in toReorder according to the order of their names in ##
##   orderedNames                                                           ##
##############################################################################

reorderElements <- function(toReorder, orderedNames, what="") {
  # reorder elements
  if (!is.null(toReorder)) {
    #verify inputs correctness
    sapply(names(toReorder), function(name) {
      if (!(name %in% orderedNames)) {
        errorMsg <- paste0("Incorrect model structure detected when reordering elements in ",
                           what, ". \n Element: ", name, 
                           "is not in the list of ordered names.")
      }
    })
    return(sapply(1:length(orderedNames), function(i) toReorder[which(names(toReorder)==orderedNames[i])]))
  }
  NULL
}

##########################################################################
##                           Internal function                          ##
## Search circularities in the system definition                        ##
## WARNING: it is a recursive algorithm and hence might fail for        ##
## excessively large systems, returning a "C stack usage" error         ##
## If this happens, contact me and I will make a sequential version     ##
## of the algorithm                                                     ##
##########################################################################
tarjan <- function(adjacency) {
  if (length(which(diag(adjacency)!=0)!=0)) {
    errorMsg <- paste0("Logical error in model definition: variable ", colnames(adjacency)[which(diag(adjacency)!=0)][1], " is a function of itself !")
    stop(errorMsg)
  }
  n <- ncol(adjacency)
  indexVec <- rep(-1, n) #-1 for unexplored indexes
  # indexVec[which(colSums(adjacency)==0)] <- 0 #0 for points that depend on no variables (ie sinkcells)
  lowerIndexVec <- rep(-1, n)
  # lowerIndexVec[which(colSums(adjacency)==0)] <- 0 #0 for points that depend on no variables (ie sinkcells)
  onStack <- rep(FALSE, n)
  nLoopFound <- 0
  allLoops <- list()
  index <- 0
  currentPoints <- c()
  
  searchLoop <- function(j) {
    index <<- index+1
    indexVec[j] <<- index
    lowerIndexVec[j] <<- index
    currentPoints <<- c(currentPoints, j)
    onStack[j] <<- TRUE
    for (i in 1:n) {
      if (adjacency[i, j]!=0) { # For all variables i that j depends on 
        if (indexVec[i]==-1) { #if i not explored yet
          searchLoop(i)
          lowerIndexVec[j] <<- min(lowerIndexVec[j], lowerIndexVec[i])
        } else {
          if (onStack[i]==TRUE) { #if i is already on stack, i.e. we found a new circularity (i -> ... -> j -> i)
            lowerIndexVec[j] <<- min(lowerIndexVec[j], indexVec[i]) #this line is CORRECT, we take indexVec[j] and not lowerIndexVec[j] on purpose (see wiki on Tarjan alfgorithm)
          }
          # If we found a point already explored, but not on stack, then j is an edge pointing towards a circularity previously identified (or a variable properly defined), and thus can be ignored
        }
      }
    }
    if (lowerIndexVec[j]==indexVec[j]) { #Root node found (i.e. starting point from a circularity)
      pointsCircularity <- c()
      repeat {
        i <- currentPoints[length(currentPoints)]
        currentPoints <<- currentPoints[-length(currentPoints)]
        onStack[i] <<- FALSE
        pointsCircularity <- c(pointsCircularity, i)
        if (i==j) break
      }
      if(length(pointsCircularity)>1) { # ignore singletons
        nLoopFound <<- nLoopFound + 1
        allLoops[[nLoopFound]] <<- pointsCircularity
      }
    }
  }
  
  for (j in 1:n) {
    if(indexVec[j]==-1) {
      searchLoop(j)
    }
  }
  
  return(allLoops)
}
##########################################################################
##                           Internal function                          ##
## Reorder equations (both intermediate variables and time derivatives) ##
## And detect circularities (which cause an error and stop the code)    ##
##########################################################################
unloopEquations <- function(tderiv, intermediateVar = NULL) {
  allEquations <- tderiv
  names(allEquations) <- paste0(names(allEquations), "Dot")
  if (!is.null(intermediateVar))
    allEquations <- c(allEquations, intermediateVar)
  adjacency <- createAdjacency(allEquations)
  variablesRemaining <- rep(TRUE, length(allEquations))
  
  unloopedOrder <- c()
  circularities <- tarjan(adjacency)
  if (length(circularities)>0) {
    errorMsg <- "Logical error in model definition: circularities detected. \n"
    for (i in 1:length(circularities)) {
      errorMsg <- c(errorMsg, "Circularity ", i, ": ", paste(names(allEquations)[circularities[[i]]], collapse="; "), "\n")
    }
    stop(errorMsg)
  }
  while (any(variablesRemaining==TRUE)) {
    tempSum <- colSums(as.matrix(adjacency[variablesRemaining, variablesRemaining]))
    newUnloopedVariables <- (1:length(allEquations))[variablesRemaining][which(tempSum==0)]
    if (length(newUnloopedVariables)==0 || is.null(newUnloopedVariables)) {
      stop("If you see this message, it means that the Tarjan Algorithm failed to detect a circularity, which should never happen in theory. Please signal this to the library maintainer at stanislas.augier@protonmail.com. \n", 
           "Error message: bad model specification: there is a circularity with the variables: ",
           paste(names(allEquations)[variablesRemaining], collapse=" "))
    }
    unloopedOrder <- c(unloopedOrder, newUnloopedVariables)
    variablesRemaining[unloopedOrder] <- FALSE
  }
  return(allEquations[unloopedOrder])
}

#################################################################################
##                           Internal function                                 ##
## initialize vector of initial positions y0                                   ##
#################################################################################

computeInitialPosition <- function(sys) {
  y0Text <- sys$y0Text; parms <- sys$parms; tderivText <- sys$tderiv; intermediateVarText <- sys$intermediateVar
  if(!is.null(sys$samplesExogVar)) { # replace exogenous variables by their value at tInit for computation of y0
    for (i in 1:length(sys$samplesExogVar)) {
      intermediateVarText[which(names(intermediateVarText)==names(sys$samplesExogVar)[i])] <- as.character(sys$samplesExogVar[[i]][1])
    }
  }
  out <- rep(0, length(y0Text))
  names(out) <- names(y0Text)
  #INIT
  y0Call <- tderivCall <- intermediateVarCall <- c()
  for (i in 1:length(y0Text)) y0Call <- c(y0Call, parse(text=y0Text[i])[[1]])
  for (i in 1:length(tderivText)) tderivCall <- c(tderivCall, parse(text=tderivText[i])[[1]])
  if(!is.null(intermediateVarText)) for (i in 1:length(intermediateVarText)) intermediateVarCall <- c(intermediateVarCall, parse(text=intermediateVarText[i])[[1]])
  names(y0Call) <- names(y0Text)
  names(tderivCall) <- sapply(names(tderivText), function(name) paste0(name, "Dot"))
  if(!is.null(intermediateVarText)) names(intermediateVarCall) <- names(intermediateVarText)
  
  # Search for correct order to define initial position
  allEquations <- c(y0Call, tderivCall, intermediateVarCall)
  adjacency <- createAdjacency(allEquations)
  variablesRemaining <- rep(TRUE, length(allEquations))
  
  unloopedOrder <- c()
  circularities <- tarjan(adjacency)
  if (length(circularities)>0) {
    errorMsg <- "Logical error in model definition: circularities detected when defining initial variables values (y0). \n"
    for (i in 1:length(circularities)) {
      errorMsg <- c(errorMsg, "Circularity ", 
                    i, 
                    ": ",
                    paste(names(allEquations)[circularities[[i]]], collapse="; "),
                    "\n")
      errorMsg <- c(errorMsg, 
                    "Involving initial variables values: ", 
                    paste(names(allEquations)[circularities[[i]]][which(names(allEquations)[circularities[[i]]] %in% names(y0Text))], collapse="; "),
                    "\n")
    }
    for (i in 1:length(circularities)) {

    }
    stop(errorMsg)
  }
  while (any(variablesRemaining==TRUE)) {
    tempSum <- colSums(as.matrix(adjacency[variablesRemaining, variablesRemaining]))
    newUnloopedVariables <- (1:length(allEquations))[variablesRemaining][which(tempSum==0)]
    if (length(newUnloopedVariables)==0 || is.null(newUnloopedVariables)) {
      stop("If you see this message, it means that the Tarjan Algorithm failed to detect a circularity, which should never happen in theory. Please signal this to the library maintainer at stanislas.augier@protonmail.com. \n", 
           "Error message: bad model specification: there is a circularity with the variables: ",
           paste(names(allEquations)[variablesRemaining], collapse=" "))
    }
    unloopedOrder <- c(unloopedOrder, newUnloopedVariables)
    variablesRemaining[unloopedOrder] <- FALSE
  }
  currentEnv <- environment()
  useless <- sapply(1:length(parms),
         function(i) eval(parse(text=paste0(names(parms)[i], " = ", parms[i]))[[1]], envir=currentEnv))
  for (i in unloopedOrder) { #evaluate all equation in the unlooped order
    # evaluate equation in current environment
    myCall <- allEquations[i]
    names(myCall) <- NULL
    myCall <- call("=", parse(text=names(allEquations)[i])[[1]], myCall[[1]])
    eval(myCall, envir=currentEnv)
    if (names(allEquations)[i] %in% names(y0Call)) { #if equation defines an initial position, update out vector
      myCall <- allEquations[i]
      names(myCall) <- NULL
      myCall <- call("=",
                     call("[", 
                          quote(out),
                          names(allEquations)[i]),
                     myCall[[1]])
      eval(myCall, envir=currentEnv)
    }
  }
  out
}


#################################################################################
##                           Internal function                                 ##
## Create adjacency tables for the different variables                         ##
## adjacency[i, j] tells us if the variable j depends on the variable i or not ##
#################################################################################

searchVariableInEq <- function(Eq, variable) {
  if (length(Eq)>1) {
    outVec <- sapply(Eq, searchVariableInEq, variable=variable)
    return(any(outVec==TRUE))
  } else {
    return(as.character(Eq)==variable)
  }
}
createAdjacency <- function(allEquationsStr) {
  varNames <- names(allEquationsStr)
  adjacency = matrix(0, nrow = length(allEquationsStr), ncol = length(allEquationsStr), dimnames=list(varNames , varNames))
  for (j in 1:length(allEquationsStr)) {
    Eq <- parse(text=allEquationsStr[j])[[1]]
    for (i in 1:length(allEquationsStr)) { #set adjacency i, j to TRUE if variable i is in equation j (WARNING, we will rotate adjacency after) 
      adjacency[j, i] <- 1*searchVariableInEq(Eq=Eq, variable=varNames[i])
    }
  }
  t(adjacency)
}
##################################################
##            Internal function                 ##
## Create str for the cpp function of the model ##
##################################################
makeCppFunc <- function(allEquations, parms, tderiv, intermediateVar) { # Create the cpp function defining differential equations for the kernel
  # init the string with inputs
  strCppFunc <- 'void Func(const T t, const T* y, const T* parms, T* ydot, T* x) override {'
  strCppFunc <- paste0(strCppFunc, '\n')
  
  # (Pre)Format the equations
  allEquations <- sapply(allEquations,
                         function(x) {
                           funcFormatEquation(Eq=parse(text=x)[[1]], parmsNames=names(parms), 
                                              yNames=names(tderiv), ivNames = names(intermediateVar))
                         }
  )
  
  # for (i in 1:length(allEquations)) {
  #   x <- allEquations[i]
  #   funcFormatEquation(Eq=parse(text=x)[[1]], parmsNames=names(parms), 
  #                      yNames=names(tderiv), ivNames = names(intermediateVar))
  # }
  
  # Include the equations to the str
  for(i in 1:length(allEquations)) {     # add one line with "varName = varExpression"
    if(any(paste0(names(tderiv), "Dot")==names(allEquations[i]))==TRUE) { # if it is a time derivative
      indiceVar <- which(paste0(names(tderiv), "Dot")==names(allEquations[i]))
      strCppFunc <- paste(strCppFunc, 
                          paste0(paste0("ydot[", indiceVar-1, "]"),
                                 ' = ', 
                                 allEquations[i], 
                                 ';'),
                          sep="\n")
    } else { # If it is an intermediate variable
      indiceVar <- which(names(intermediateVar)==names(allEquations[i]))
      strCppFunc <- paste(strCppFunc, 
                          paste0(paste0("x[", indiceVar-1, "]"),
                                 ' = ', 
                                 allEquations[i], 
                                 ';'),
                          sep="\n")
    }
  }
  strCppFunc <- paste(strCppFunc, '}', sep="\n")
}
##########################################################################
##                           Internal function                          ##
## Create str for the events triggered by a condition on time           ##
##########################################################################
makeEventTime <- function(sys=NULL, eventTime=sys$eventTime,
                          parmsNames=names(sys$parms), yNames = names(sys$tderiv), ivNames=names(sys$intermediateVar)) {
  
  allNames <- c("parmsNames", "yNames")
  missingArguments(allNames, envir=environment())
  
  out <- 'void makeEventTime(const T t, T* parms, T* y, T* x, T h) override {\n'
  if (!is.null(eventTime)) {
    for (i in 1:length(eventTime)) {
      newEventTime <- paste0("if(abs(t - ",
                             as.character(eventTime[[i]][[1]]),
                             ")<0.5*h) { \n")
      for (j in 2:length(eventTime[[i]])) {
        newEventTime <- paste0(newEventTime,
                               deparse(funcFormatEquation(Eq=parse(text=paste0(names(eventTime[[i]])[j],
                                                                               "=", 
                                                                               eventTime[[i]][[j]]))[[1]],
                                                          parmsNames=parmsNames,
                                                          yNames=yNames, ivNames=ivNames)),
                               "; \n")
      }
      out <- paste0(out, newEventTime, "} \n")
    }
  }
  out <-  paste0(out, "} \n")
  out
}
##########################################################################
##                           Internal function                          ##
## Create str for the events triggered by a condition on variables      ##
##########################################################################
makeEventVar <- function(sys=NULL, eventVar=sys$eventVar, 
                         parmsNames=names(sys$parms), yNames = names(sys$tderiv), 
                         ivNames=names(sys$intermediateVar)) {
  allNames <- c("parmsNames", "yNames")
  missingArguments(allNames, envir=environment())
  
  out <- 'void makeEventVar(const T t,T* parms, T* y, T* x, T h) override {\n'
  if (!is.null(eventVar)) {
    for (i in 1:length(eventVar)) {
      newEventVar <- paste0("if(",
                            deparse(funcFormatEquation(Eq=parse(text=as.character(eventVar[[i]][[1]]))[[1]], 
                                                       parmsNames=parmsNames,
                                                       yNames=yNames, ivNames=ivNames)),
                            ") { \n")
      for (j in 2:length(eventVar[[i]])) {
        newEventVar <- paste0(newEventVar,
                              deparse(funcFormatEquation(Eq=parse(text=paste0(names(eventVar[[i]])[j], 
                                                                              "=", 
                                                                              eventVar[[i]][[j]]))[[1]],
                                                         parmsNames=parmsNames,
                                                         yNames=yNames, ivNames=ivNames)),
                              "; \n")
      }
      out <- paste0(out, newEventVar, "} \n")
    }
  }
  out <- paste0(out, "} \n")
  out
}
######################################################
##                  Internal function               ##
## Format the equations to translate them in C code ##
######################################################
funcFormatEquation <- function(Eq, parmsNames, yNames, ivNames) {
  
  if (length(Eq)>1) { #Eq is a call with multiple elements
    if (Eq[[1]]=="getIthExogVar") return(Eq) #Ecogenous variable, nothing to to
    for (i in 1:length(Eq)) { #apply the function to all the elements
      Eq[[i]] <- funcFormatEquation(Eq[[i]], parmsNames, yNames,ivNames)
    }
    if (Eq[[1]]==quote(ifelse)) { # Handle ifelse statement, use temporary code that is translated during the last step
      newEq <- call("somethingToRemoveLater")
      newEq[[2]] <- call("(", Eq[[2]])
      newEq[[3]] <- quote(MakeThisAndCommaAnInterrogationMark)
      newEq[[4]] <- Eq[[3]]
      newEq[[5]] <- quote(MakeThisAndCommaAColon)
      newEq[[6]] <- Eq[[4]]
      Eq <- newEq
    }
    else if (Eq[[1]]==quote(max)) {
      Eq[[1]] <- quote(makeThisAnSTDAndDoubleColonmax)
    }
    else if (Eq[[1]]==quote(min)) {
      Eq[[1]] <- quote(makeThisAnSTDAndDoubleColonmin)
    }
  } else { #Single element
    if (any(parmsNames==as.character(Eq))) { # Eq is a parameter, replace it by its position in the vector parms
      parmsIndice <- which(parmsNames==rep(as.character(Eq), length(parmsNames))) - 1 
      Eq <- paste0("parms[", as.character(parmsIndice), "]")
      Eq <- parse(text=Eq)[[1]]
    } else if (any(yNames==as.character(Eq))) { # Eq is a variable, replace it by its position in the vector y
      yIndice <- which(yNames==as.character(Eq)) - 1
      Eq <- paste0("y[", as.character(yIndice), "]")
      Eq <- parse(text=Eq)[[1]]
    } else if (any(paste0(yNames, "Dot")==as.character(Eq))) { # Eq is a time derivative, repace it by its position in the vector ydot
      yIndice <- which(paste0(yNames, "Dot")==as.character(Eq)) - 1
      Eq <- paste0("ydot[", as.character(yIndice), "]")
      Eq <- parse(text=Eq)[[1]]
    } else if (any(ivNames==as.character(Eq))) { # Eq is an intermediate variable, replace it by its position in the vector x
      ivIndice <- which(ivNames==as.character(Eq)) - 1
      Eq <- paste0("x[", as.character(ivIndice), "]")
      Eq <- parse(text=Eq)[[1]]
    } else if(is.numeric(Eq)) { # Numeric value, need to make it a float! Mark it and make it a float in the final step
      tempEq <- Eq
      if (as.integer(Eq)==Eq) { #if it is an int
        Eq <- as.name(paste0("somethingToRemoveLater", as.character(abs(Eq)), "intToFloat"))
      } else { # if it is a double
        Eq <- as.name(paste0("somethingToRemoveLater", as.character(abs(Eq)), "doubleToFloat"))
      }
      if (tempEq<0) { #Negative value, add a "-"
        Eq <- call("-", Eq)
      }
    } else if (Eq=="^" || Eq=="**") { #Exponent, need to use pow, directly replace it here
      Eq <- quote(pow)
    }
  }
  return(Eq)
}

######################################################################
##                         Internal Function                        ##
##  Turns lists/matrices into one dimensional vectors for C++ code  ##
######################################################################
formatExogVar <- function(sys) {
  if (!is.null(sys$samplesExogVar) && length(sys$samplesExogVar)>0) {
    samplesExogVarCpp  <- unlist(sys$samplesExogVar); names(samplesExogVarCpp) <- NULL
    nSamplesVarExogVarCpp <- sapply(sys$samplesExogVar, length); names(nSamplesVarExogVarCpp) <- NULL
    nVarExogVarCpp <- length(sys$samplesExogVar)
  } else {
    samplesExogVarCpp <- c(0)
    nSamplesVarExogVarCpp <- c(0)
    nVarExogVarCpp <- 0
  }  
  return(list(samplesExogVarCpp = samplesExogVarCpp, 
              nSamplesVarExogVarCpp = nSamplesVarExogVarCpp,
              nVarExogVarCpp = nVarExogVarCpp))
}


#############################################
##           Internal function             ##
## Edit the preprocRCPP.h file to define   ##
## compile-time options                    ##
#############################################
editPreprocRCPP <- function(sys) {
  allNames <- c("tderiv", "returnRK4")
  missingArguments(allNames, sys)
  verboseCMAES <- sys$verboseCMAES
  # initScipen <- getOption("scipen")
  # options(scipen=999) #avoid scientific notation, that ca cause issues with the different conversions from R to cpp
  
  
  # Init objects to pass to Cpp
  dim <- length(sys$tderiv)
  dimIv <- length(sys$intermediateVar) #will be zero if sys$intermediateVar is NULL
  nExogVar <- length(sys$samplesExogVar)
  useDistForConvergence <- ifelse(sys$BruteForceBasinConvCrit=="dist", 1, 0)
  useParallel <- ifelse(sys$useParallel==0, 0, 
                        ifelse(sys$useParallel==1, 1, 2))
  useExogVar <- ifelse(length(sys$samplesExogVar)>0, 1, 0)
  useEventTime <- ifelse(is.null(sys$eventTime), 0, 1)
  useEventVar <- ifelse(is.null(sys$eventVar), 0, 1)
  returnRK4 <- sys$returnRK4
  longIntForHCMSubdivision <- sys$longIntForHCMSubdivision
  #Load Raw Cpp Code and edit it
  cppCode <- readChar("SourceCodeParallel/cppCode/src/preprocRcpp.h", file.info("SourceCodeParallel/cppCode/src/preprocRcpp.h")$size)
  cppCode <- sub("@ADDuseParallel", useParallel, cppCode)
  cppCode <- sub("@ADDuseEventTime", useEventTime, cppCode)
  cppCode <- sub("@ADDuseEventVar", useEventVar, cppCode)
  cppCode <- sub("@ADDreturnRK4", returnRK4, cppCode)
  cppCode <- sub("@ADDTInt", longIntForHCMSubdivision, cppCode)
  cppCode <- sub("@ADDVerboseCMAES", verboseCMAES, cppCode)
  if (!is.null(sys$compileForPython==TRUE)) {
    cppCode <- sub("@ADDNtForPython", sys$nt, cppCode)
    cppCode <- sub("@ADDTInitForPython", sys$tInit, cppCode)
    cppCode <- sub("@ADDTEndForPython", sys$tEnd, cppCode)
    cppCode <- sub("@ADDNVForPython", length(sys$tderiv), cppCode)
    cppCode <- sub("@ADDNIVForPython", length(sys$intermediateVar), cppCode)
    cppCode <- sub("@ADDVarNamesForPython",
                   paste0('{',
                          paste0(sapply(names(sys$tderiv), 
                                        function(x) paste0('"', x, '"')),
                                 collapse=", "),
                          '}')
                   , cppCode)
    cppCode <- sub("@ADDIntermediateVarNamesForPython",
                   paste0('{',
                          paste0(sapply(names(sys$intermediateVar), 
                                        function(x) paste0('"', x, '"')),
                                 collapse=", "),
                          '}')
                   , cppCode)
    cppCode <- sub("@ADDParmsNamesForPython",
                   paste0('{',
                          paste0(sapply(names(sys$parms), 
                                        function(x) paste0('"', x, '"')),
                                 collapse=", "),
                          '}')
                   , cppCode)
    
    cppCode <- sub("@ADDYInitForPython",
                   paste0('{',
                          paste0(sys$y0, collapse=", "),
                          '}')
                   , cppCode)
    cppCode <- sub("@ADDParmsForPython",
                   paste0('{',
                          paste0(sys$parms, collapse=", "),
                          '}')
                   , cppCode)
    formatedExogVar <- formatExogVar(sys)
    cppCode <- sub("@ADDSamplesExogVarForPython",
                   paste0('{',
                          paste0(formatedExogVar$samplesExogVarCpp, collapse=", "),
                          '}')
                   , cppCode)
    cppCode <- sub("@ADDNSamplesVarExogVarForPython",
                   paste0('{',
                          paste0(formatedExogVar$nSamplesVarExogVarCpp, collapse=", "),
                          '}')
                   , cppCode)
    cppCode <- sub("@ADDNVarExogVarForPython",
                   paste0('{',
                          paste0(formatedExogVar$nVarExogVarCpp, collapse=", "),
                          '}')
                   , cppCode)
  } else {
    cppCode <- sub("@ADDNtForPython", "0", cppCode)
    cppCode <- sub("@ADDTInitForPython", "0", cppCode)
    cppCode <- sub("@ADDTEndForPython", "0", cppCode)
    cppCode <- sub("@ADDNVForPython", "0", cppCode)
    cppCode <- sub("@ADDNIVForPython", "0", cppCode)    
    cppCode <- sub("@ADDVarNamesForPython", "0", cppCode)
    cppCode <- sub("@ADDIntermediateVarNamesPython", "0", cppCode)
    cppCode <- sub("@ADDParmsNamesForPython", "0", cppCode)
    cppCode <- sub("@ADDYInitForPython", "0", cppCode)
    cppCode <- sub("@ADDSamplesExogVarForPython", "0", cppCode)
    cppCode <- sub("@ADDBSamplesVarExogVarForPython", "0", cppCode)
    cppCode <- sub("@ADDNVarExogVarForPython", "0", cppCode)
  }
  # options(scipen=initScipen)
  writeChar(object=cppCode, "SourceCodeParallel/cppCode/src/preprocRCPP_R.h", nchars=nchar(cppCode), eos=NULL)
}

createFuncCpp <- function(sys) {
  allNames <- c("strFunc", "strEventTime", "strEventVar")
  missingArguments(allNames, sys)
  #Load Raw Cpp Code and edit it
  cppCode <- readChar("SourceCodeParallel/cppCode/src/modelR_raw.h", file.info("SourceCodeParallel/cppCode/src/modelR_raw.h")$size)
  cppCode <- sub("@ADDFunc", sys$strFunc, cppCode)
  cppCode <- sub("@ADDEventVar", sys$strEventVar, cppCode)
  cppCode <- sub("@ADDEventTime", sys$strEventTime, cppCode)

  # Final formatting
  cppCode <- gsub("intToFloat", ".0", cppCode) #double not float
  cppCode <- gsub("doubleToFloat", "", cppCode) #double not float
  cppCode <- gsub("somethingToRemoveLater", "", cppCode)
  cppCode <- gsub("thisIsAnInt", "", cppCode)
  cppCode <- gsub(", MakeThisAndCommaAnInterrogationMark,", " ? ", cppCode)
  cppCode <- gsub(", MakeThisAndCommaAColon,", " : ", cppCode)
  cppCode <- gsub("makeThisAnSTDAndDoubleColon", "std::", cppCode)
  cppCode <- gsub("MakeThisARightArrow", "->", cppCode)
  
  
  #write the result
  writeChar(object=cppCode, "SourceCodeParallel/cppCode/src/modelR.h", nchars=nchar(cppCode), eos=NULL)
}


#################################################
## Utility function                            ##
## Change the value of the selected parameters ##
#################################################
changeParameters <- function(sys,toChange) {
  parms = sys$parms
  for(i in 1:length(toChange))
    parms[which(names(parms)==names(toChange)[i])]=as.numeric(toChange[i])
  return(parms)
}

