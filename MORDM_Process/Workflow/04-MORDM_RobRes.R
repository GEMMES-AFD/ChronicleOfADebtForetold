cat("Calling Parallelisation Libraries...\n")
source("SourceCodeParallel/SourceCode.R")
#Put the equations into the algorithm
eventTime        <- NULL  #FALSE or NULL/TRUE
eventVar         <- TRUE  #FALSE or NULL/TRUE
useParallel      <- TRUE   #FALSE or NULL/TRUE
compileForPython <- TRUE  #FALSE or NULL/TRUE
verboseCMAES     <- 0      #0 (no verbose), 1 (a bit of verbose), 2 (lots of verbose), 3 (too much verbose)
#NOTE: Setting verbose!=0 might reduce performances, use only for debug/calibration of CMAES
cat("Compiling...\n")
SOEM2 <- initSysCpp(
  fileName=MORDM_FILE,
  useParallel=useParallel,
  verboseCMAES=verboseCMAES)

cat("Setting offset for parallelisation...\n")
SOEM2$parms["tOffset"] <- 2019   # parallel mode: temp = t - 2019

if (PARALLELVANILLATEST){
  cat("Testing for consistency between vanilla and parallel...\n")
i = 1
params_all <- readRDS(robfile)[[i]]
params   <-params_all[1,]
parms_NewC2 <- recal(params, parms_NewC)
###Check for error
parms_NewC2["betaen"] == params["betaen"]

parms_NewC2['lambdatr0'] = 5                #Speed of the NDC investment path
parms_NewC2['lambdatr1']  = 6               #Initial period of the NDC investment path
parms_NewC2['lambdatr2']  = 0.011           #Target NDC investment as a share of NFC's capital stock in 2019
parms_NewC2['lambdatr0_adj'] = 5            #Speed of the NDC investment path
parms_NewC2['lambdatr1_adj']  = 12          #Initial period of the NDC investment path
parms_NewC2['lambdatr2_adj']  = 1           #Target NDC investment as a share of NFC's capital stock in 2019

for (shape in shapes){ # 1: Smooth transition, 2: Sharp Transition; 3: Delayed-Smooth; #4 Very quick
  
  if (shape == 1){
    parms_NewC['alpha_tr']=0.000709889
    parms_NewC['beta_tr']=0.09
    parms_NewC['gamma_tr']=1.6
    parms_NewC['delta_tr']=4.99377*10^-7
    event1 <- list(triggerDate=4, reducXrO=reducXro)
  } else if (shape == 2) {
    parms_NewC['alpha_tr']=0.00104006
    parms_NewC['beta_tr']=0.3
    parms_NewC['gamma_tr']=2.6
    parms_NewC['delta_tr']=4.99377*10^-7
    event1 <- list(triggerDate=4, reducXrO=reducXro)
  } else if (shape == 3){
    parms_NewC['alpha_tr']=0.00075
    parms_NewC['beta_tr']=0.18
    parms_NewC['gamma_tr']=2.11
    parms_NewC['delta_tr']=4.99377*10^-7
  } else if (shape == 4) {
    parms_NewC['alpha_tr']=0.09
    parms_NewC['beta_tr']=0.6
    parms_NewC['gamma_tr']=0.9
    parms_NewC['delta_tr']=4.99377*10^-7
  }
  
}

reduc = params["reducXrO"]
event1 <- list(triggerDate=4, reducXrO=reduc)
res <- cppRK4(SOEM, parms= parms_NewC2, times=seq(from=2019, to=2050, by=0.1), eventTime=list(event1))
res_par <- cppSA(sys=SOEM2, times=seq(from=2019, to=2050, by=0.1), allParmsSA=params_all, fullTrajectory = TRUE, solver="RK4Fixed")[1,,]

err <- max((abs(as.matrix(res["gdp"]) - as.matrix(res_par[,"gdp"])))/as.matrix(res["gdp"]))

if (err < 0.5e-1){#Relatively lagre tolerance due to make for RUNGEKUTTA numerical imprecisions
  cat("Parallel and vanilla are consistent!\n")
} else {
  stop("Parallel and vanilla are NOT consistent!")
}

}

cat("Running robustness simulations...\n")
###Getting shape 
for (shape in shapes){
  
  ##
  data <- readRDS(OPTPOL_FILE)
  
  length_robfile <- length(readRDS(ROBCAL_FILE))
  
  if (nrow(data[[1]]) != length_robfile){
    print("Optimal policy number does not match numer of robustness calibrations")
    break
  }
  
  varNames <-c('perCapita','inflation','reserves','foreignDebt','privateDebt','pubDebt','fiscalDef','hhFrag','firmsFrag','unem','CAD', "Gip", "GDP")
  
  roblist <- vector("list",length_robfile)
  
  orig.set <- mordm.get.set(data)
  #roblist <- list()
  options(mc.cores=n_cores)
  pb <- txtProgressBar(min=0, max=nrow(data[[1]]), style=3)
  for (i in 1:nrow(data[[1]])){
    setTxtProgressBar(pb, i)
    params   <-readRDS(ROBCAL_FILE)[[i]]
    SOEM2$parms["dsactive"] = 1
    SOEM2$parms["reaction"] = 1
    policies <- unique(params[,8:13])
    outcomes <- cppSA(sys=SOEM2, times=seq(from=2019, to=2050, by=0.1), allParmsSA=params, fullTrajectory = TRUE, solver="RK4Fixed")[,,varNames]
    gc()
    
    
    a <-lapply(1:dim(outcomes)[1], function(kk){
      failure = F
      slice <- outcomes[kk,,]
      slice <- slice[11:110,]
      slice[,1] <- (slice[,1] - lag(slice[,1]))/lag(slice[,1])
      
      slice[,12] <- slice[,12]/slice[,13]
      slice <- slice[,-13]
      
      colnames(slice)[1] <- "growth"
      
      
      
      
      if (max((slice[,'reserves'] <0)>0, na.rm=T)){
        failure = T
      } else {
        failure =F
      }
      
      if (failure==1){
        if (length(which(slice[,'reserves'] <0)) >0 | length(which(is.na(slice[,'reserves']))) >0){
          crisistime = min(which(slice[,'reserves'] <0), min(which(is.na(slice[,'reserves']))))
        }
      }
      
      costs <- matrix(ncol = 12, nrow = 1)
      
      for (m in 1:12){
        
        
        col <- varNames[m]
        variable <- slice[,m]
        lbound <- threshold[m, 2]
        hbound <- threshold[m, 3]
        direction <- threshold[m,4]
        center <- threshold[m, 5]
        
        ####Indicator for below or above center -- to avoid over-paying bonus if within the SOP
        up <- ifelse(variable > center, 1, 0)
        ###Indicator for above the higher bound
        above_higher <- ifelse(variable >= hbound, 1, 0)
        ###Indicator for above the lower   bound
        below_lower <- ifelse(variable <= lbound, 1, 0)
        
        if (direction == 0){
          
          
          ###Computing time-wise performance
          timeperformance <- ifelse(up > 0,
                                    ifelse(above_higher == 1, 
                                           (100*abs((variable - hbound))/(hbound-lbound)),
                                           -sqrt(100*abs(hbound - variable)/(hbound-lbound))
                                    ),
                                    ifelse(below_lower == 1,
                                           (100*abs((variable - lbound))/(hbound-lbound)),
                                           -sqrt(100*abs(variable - lbound)/(hbound-lbound))
                                    )
          )
        } else if (direction == -1){
          timeperformance = ifelse(above_higher ==0, 100*abs((variable - hbound)), -sqrt(100*abs(variable - hbound)/(hbound-lbound)))
          
        } else if (direction == 1) {
          timeperformance = ifelse(above_higher == 1, 
                                   (100*abs((variable - hbound))/(hbound-lbound)),
                                   -sqrt(100*abs((hbound - variable))/(hbound-lbound)))
        }
        
        
        timeperformance <- timeperformance[which(!is.infinite(timeperformance))]
        
        
        
        # for (k in 1:length(timeperformance)){
        #   if (failure){
        # timeperformance[k] <- min(timeperformance[k],100^2)*nrow(runfile)/crisistime
        #} else {
        # timeperformance[k] <- min(timeperformance[k],100^2)
        #}
        #}
        
        for (k in 1:length(timeperformance)){
          if (failure){
            timeperformance[k] <- 100^2*nrow(slice)/crisistime
          } else {
            timeperformance[k] <- min(timeperformance[k],100^2)
          }
        }
        
        costs[,m]  <-  sum(timeperformance, na.rm = T)/(length(which(!is.na(timeperformance)))) 
        
      }
      
      colnames(costs) <- colnames(slice)
      
      final <- cbind(t(params[i,8:13]), costs, failure)
      colnames(final)[ncol(final)] <- "Failure"
      final })
    
    cnames <- colnames(a[[1]])
    
    a <- matrix(unlist(a), ncol = ncol(a[[1]]), byrow = TRUE)
    colnames(a) <- cnames
    
    roblist[[i]] <- a
    
    rm(a)
    
    rm(outcomes)
  }
  close(pb)

  saveRDS(roblist, ROBLIST_FILE)
}
