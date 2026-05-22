library(OpenMORDM)
library(lhs)
library(rdyncall)
library(parallel)


#load(".RData")
#Disabling warnings 
options(warn = -1) 
library(readxl)
library(xtable)
library(rlist)

library("data.table")
library("ggplot2")

#Load the package
library(Rcpp)
library(rmarkdown)
library(lhs)
library(readr)
library(tidyverse)

#Load the algorithm
source("Source/SourceCode.R")
source("Source/utilities.R")
source("Extrafunctions.R")

alts <- c("Low", "High")
for (alt in alts){
###Calling parameters -- Five best
alt_bas <- read_delim(paste0("Data/X",alt,"_ord.csv"), delim = ";", 
                      escape_double = FALSE, trim_ws = TRUE) %>%
  slice(1:5)


###Generating System
xr0 = 0.025
event1 <- list(triggerDate=4, reducXrO=xr0)
SOEM <- cppMakeSys(fileName = "model_equations_MORDM.R",reportVars=3, eventTime = list(event1))


varNames<-c('perCapita','inflation','reserves','foreignDebt','privateDebt','pubDebt','fiscalDef','hhFrag','firmsFrag','unem','CAD', 'Gip', 'GDP')


###Creating SOS based on baseline performances
threshold <- as.data.frame(matrix(nrow=12, ncol = 4))
colnames(threshold) <- c("Variables", "L_Bound", "H_Bound", "direction")
threshold[,1] <- c("growth", "inflation", "reserves", "foreignDebt", "privateDebt",  "pubDebt",  "fiscalDef",    "hhFrag", "firmsFrag", "unem","CAD", "Gip")
threshold[1,2:4] <- c(0,0.2, 0)###Per capita GDP growth
threshold[2,2:4] <- c(0.01,0.04, 0)###Inflation
threshold[3,2:4] <- c(0,0.18, -1)  ###Reserves
threshold[4,2:4] <- c(0,0.4, 1) ###Foreign Debt Ratio
threshold[5,2:4] <- c(0,0.9,1)  ###Private Debt ratio
threshold[6,2:4] <- c(0,0.63, 1) ###Public debt ratio
threshold[7,2:4] <- c(0,0.035, 1)###Fiscal Deficit
threshold[8,2:4] <- c(0,0.45, 1) ###Household fragility Index
threshold[9,2:4] <- c(0,0.26, 1)  ###Firm fragility index
threshold[10,2:4] <- c(-0.01,0.12, 0)   ###Unemployment
threshold[11,2:4] <- c(0,0.05, 1)   ###Current Account Deficit
threshold[12,2:4] <- c(0,0.04, 1)   ###Government interest

threshold <- threshold %>%
  dplyr::mutate(center = ifelse(direction == 0, 0.5*(L_Bound+H_Bound), 0))



plain_met <- function(levers){
  
  costs <- matrix(ncol=11)
  
  if (length(levers)!=3){
    stop("Expected exactly 3 decision variables.")
  }
  
  parms_NewC["shrGrL"]<-levers[1]
  parms_NewC["shrGrLFx"]<-levers[2]
  parms_NewC["md_lgtr"]<- levers[3]
  parms_NewC["md_bgtr"]<-levers[3]
  
  event1 <- list(triggerDate=4, reducXrO=0.01)
  
  
  newScen <- cppRK4(SOEM, parms= parms_NewC,times=seq(from=2019, to=2050, by=0.1), eventTime=list(event1))
  
  newScen <- newScen[,varNames] %>% as.data.frame()
  
  
  
  newScen <- newScen %>% 
    dplyr::mutate(perCapita = 100*diff(perCapita)/lag(perCapita)) %>%
    filter(row_number() < 301 & row_number() > 50)
  
  
  
  if (max(newScen$reserves <0)>0){
    failure = T
  } else {
    failure =F
  }
  
  if (failure==1){
    if (length(which(newScen$reserves <0)) >0 | length(which(is.na(newScen$reserves))) >0){
      crisistime = min(which(newScen$reserves <0), min(which(is.na(newScen$reserves))))
    }
  }
  
  for (m in 1:11){
    
    
    col <- varNames[m]
    variable <- newScen[,m]
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
        timeperformance[k] <- 100^2*nrow(newScen)/crisistime
      } else {
        timeperformance[k] <- min(timeperformance[k],100^2)
      }
    }
    
    costs[,m]  <-  sum(timeperformance, na.rm = T)/(length(which(!is.na(timeperformance))))
    
  }
  constraints <- NA
  return(as.numeric(costs))
  
}


###Setting function
###Storing params
parms_NewC<-SOEM$parms
###Changing 
for (kk in 1:nrow(alt_bas)){
  for (col in colnames(alt_bas)){
    parms_NewC[col] <- alt_bas[[col]][kk]
  }
###Creating cost vector storage

costs <- matrix(ncol=nrow(threshold))


###Setting overarching NDC parameters
parms_NewC['lambdatr0'] = 5                #Speed of the NDC investment path
parms_NewC['lambdatr1']  = 6               #Initial period of the NDC investment path
parms_NewC['lambdatr2']  = 0.011           #Target NDC investment as a share of NFC's capital stock in 2019
parms_NewC['lambdatr0_adj'] = 5            #Speed of the NDC investment path
parms_NewC['lambdatr1_adj']  = 12          #Initial period of the NDC investment path
parms_NewC['lambdatr2_adj']  = 1           #Target NDC investment as a share of NFC's capital stock in 2019


####Setting NDC speed
for (shape in 2:2){ # 1: Smooth transition, 2: Sharp Transition; 3: Delayed-Smooth; #4 Very quick
  
  
  
  
  
  if (shape == 1){
    parms_NewC['alpha_tr']=0.000709889
    parms_NewC['beta_tr']=0.09
    parms_NewC['gamma_tr']=1.6
    parms_NewC['delta_tr']=4.99377*10^-7
    event1 <- list(triggerDate=4, reducXrO=0.025)
  } else if (shape == 2) {
    parms_NewC['alpha_tr']=0.00104006
    parms_NewC['beta_tr']=0.3
    parms_NewC['gamma_tr']=2.6
    parms_NewC['delta_tr']=4.99377*10^-7
    event1 <- list(triggerDate=4, reducXrO=0.025)
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
  
  
  
  #baseline_perfs <- plain_met(c(0,0,0))
  
  my_run <- function(levers){
    
    costs <- matrix(ncol=nrow(threshold))
    
    if (length(levers)!=6){
      stop("Expected exactly 6 decision variables.")
    }
    
    parms_NewC['dsactive']=1
    parms_NewC['shrGrL']= 1 
    ###Setting market reaction
    
    ###Levers
    parms_NewC['shrGrLFx'] = levers[1]
    parms_NewC['shrGrBw']= levers[2]
    parms_NewC['md_lgtr']= levers[3]
    parms_NewC['md_bgtr'] = levers[4]
    parms_NewC['mdds']= levers[5]
    parms_NewC['decds']= levers[6]
    
    event1 <- list(triggerDate=4, reducXrO=xr0)
    
    
    newScen <- cppRK4(SOEM, parms= parms_NewC,times=seq(from=2019, to=2050, by=0.1), eventTime=list(event1))
    
    newScen <- newScen[,varNames] %>% as.data.frame()
    
    
    
    newScen <- newScen %>% 
      dplyr::mutate(perCapita = 100*diff(perCapita)/perCapita) %>%
      dplyr::mutate(Gip = Gip/GDP) %>%
      dplyr::select(-GDP) %>%
      filter(row_number() < 111 & row_number() > 10)
    
    
    
    if (max(newScen$reserves <0)>0){
      failure = T
    } else {
      failure =F
    }
    
    if (failure==1){
      if (length(which(newScen$reserves <0)) >0 | length(which(is.na(newScen$reserves))) >0){
        crisistime = min(which(newScen$reserves <0), min(which(is.na(newScen$reserves))))
      }
    }
    
    for (m in 1:nrow(threshold)){
      
      col <- varNames[m]
      variable <- newScen[,m]
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
          timeperformance[k] <- 100^2*nrow(newScen)/crisistime
        } else {
          timeperformance[k] <- min(timeperformance[k],100^2)
        }
      }
      
      costs[,m]  <-  sum(timeperformance, na.rm = T)/(length(which(!is.na(timeperformance))))
      
    }
    constraints <- NA
    return(as.numeric(costs))
    
  }
  
  
  bounds <- matrix(rep(c(0, 1), 6), nrow = 2)
  ##Defining borg problem
  problem <-define.problem(my_run, nobj = 12, nvars = 6, nconstrs = 0, bounds = bounds, epsilons = rep(1e-01,12), names = c("Foreign_Loans", "Foreign _Bonds", "Greenium_Loan", "Greenium_Bonds","Interest Forgiveness", "Stock Forgiveness" ,c(varNames[-length(varNames)])))
  ###Optimiising through borg and storing
  data <- borg.optimize(problem, 1000, verbose=T)
  colnames(data[[1]]) <- c("Foreign_Loans", "Foreign _Bonds", "Greenium_Loan", "Greenium_Bonds","Interest Forgiveness", "Stock Forgiveness", varNames[-length(varNames)])
  
  filename <- paste0("Data/optpol_bastransDS",alt,"_",shape,"_robustness_",kk,".RDS")
  
  saveRDS(data, filename)
}

}

}

