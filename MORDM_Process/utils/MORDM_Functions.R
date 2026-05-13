plain_met <- function(levers, params=NULL){
  
  costs <- matrix(ncol=nrow(threshold))

  if (!(is.null(params))){
    parms_NewC <- params
  } else {
    parms_NewC <-SOEM$parms
  }
  
  if (length(levers)!=6){
    stop("Expected exactly 6 decision variables.")
  }
  
  parms_NewC['dsactive']=1
  parms_NewC['shrGrL']= 1 
  
  ###Setting reaction
  parms_NewC['reaction']=1
  parms_NewC['zetafx3']=1
  parms_NewC['rho']=0.2
  parms_NewC['gammariskFFX']=0
  parms_NewC['gammariskBFX']=0
  
  ###Levers
  parms_NewC['shrGrLFx'] = levers[1]
  parms_NewC['shrGrBw']= levers[2]
  parms_NewC['md_lgtr']= levers[3]
  parms_NewC['md_bgtr'] = levers[4]
  parms_NewC['mdds']= levers[5]
  parms_NewC['decds']= levers[6]
  
  
  
  event1 <- list(triggerDate=4, reducXrO=reducXro)
  
  
  newScen <- cppRK4(SOEM, parms= parms_NewC,times=seq(from=2019, to=2050, by=0.1), eventTime=list(event1))
  
  newScen <- newScen[,varNames] %>% as.data.frame()
  
  
  
  newScen <- newScen %>%
    mutate(perCapita = 100*diff(perCapita)/(perCapita)) %>%
    mutate(Gip = Gip/GDP) %>%
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


my_run <- function(levers){
  
  costs <- matrix(ncol=nrow(threshold))
  
  if (length(levers)!=6){
    stop("Expected exactly 6 decision variables.")
  }
  
  parms_NewC['dsactive']=1
  parms_NewC['shrGrL']= 1 
  
  ###Setting reaction
  parms_NewC['reaction']=1
  parms_NewC['zetafx3']=1
  parms_NewC['rho']=0.2
  parms_NewC['gammariskFFX']=0
  parms_NewC['gammariskBFX']=0
  
  ###Levers
  parms_NewC['shrGrLFx'] = levers[1]
  parms_NewC['shrGrBw']= levers[2]
  parms_NewC['md_lgtr']= levers[3]
  parms_NewC['md_bgtr'] = levers[4]
  parms_NewC['mdds']= levers[5]
  parms_NewC['decds']= levers[6]
  
  event1 <- list(triggerDate=4, reducXrO=reducXro)
  
  
  newScen <- cppRK4(SOEM, parms= parms_NewC,times=seq(from=2019, to=2050, by=0.1), eventTime=list(event1))
  
  newScen <- newScen[,varNames] %>% as.data.frame()
  
  
  
  newScen <- newScen %>% 
    mutate(perCapita = 100*diff(perCapita)/(perCapita)) %>%
    mutate(Gip = Gip/GDP) %>%
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

recal <- function(parms, init_par){
  for (i in 1:length(names(parms))){
    if (names(parms)[i]%in%names(init_par)){
      init_par[which(names(init_par)==names(parms)[i])] <- parms[i]
    }
  }
  return(init_par)
}
