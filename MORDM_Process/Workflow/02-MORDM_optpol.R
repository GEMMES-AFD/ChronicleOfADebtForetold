###Generating System
cat("Generating SOEM...\n")
event1 <- list(triggerDate=4, reducXrO=reducXro)
SOEM <- cppMakeSys(fileName = MORDM_FILE, reportVars=3, eventTime = list(event1))

varNames<-c('perCapita','inflation','reserves','foreignDebt','privateDebt','pubDebt','fiscalDef','hhFrag','firmsFrag','unem','CAD', 'Gip', 'GDP')

###Storing params
parms_NewC<-SOEM$parms
###Creating cost vector storage
costs <- matrix(ncol=nrow(threshold))

cat("Setting overarching NDC parameters...\n")
###Setting overarching NDC parameters
parms_NewC['lambdatr0'] = 5                #Speed of the NDC investment path
parms_NewC['lambdatr1']  = 6               #Initial period of the NDC investment path
parms_NewC['lambdatr2']  = 0.011           #Target NDC investment as a share of NFC's capital stock in 2019
parms_NewC['lambdatr0_adj'] = 5            #Speed of the NDC investment path
parms_NewC['lambdatr1_adj']  = 12          #Initial period of the NDC investment path
parms_NewC['lambdatr2_adj']  = 1           #Target NDC investment as a share of NFC's capital stock in 2019

cat("Setting NDC speed (Standard scenario)...\n")
####Setting NDC speed
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
  
  
  cat("Defining boundaries and problem...\n")
  bounds <- matrix(rep(range(0, 1), nvars), nrow = 2)
  ##Defining borg problem
  problem <-define.problem(my_run, nobj = nobj, nvars = nvars, nconstrs = nconstrs, bounds = bounds, epsilons = rep(1e-01,nobj), names = c("Foreign_Loans", "Foreign_Bonds", "Greenium_Loan", "Greenium_Bonds","Interest_Forgiveness", "Stock_Forgiveness" ,c(varNames[-length(varNames)])))
  ###Optimiising through borg and storing
  cat("Entering Borg optimisation...\n")
  data <- borg.optimize(problem, n_sim, verbose=T)
  colnames(data[[1]]) <- c("Foreign_Loans", "Foreign_Bonds", "Greenium_Loan", "Greenium_Bonds","Interest_Forgiveness", "Stock_Forgiveness", varNames[-length(varNames)])
  cat("Borg optimisation: OK!...\n")
  
  
  ####Consistency check
  cat("Running Consistency check...\n")
  test <- max(plain_met(as.numeric(data[[1]][1, 1:nvars]), parms_NewC) - data[[1]][1,col_out])
  
  if (test == 0){
    cat("No Discrepancy!\n")
    saveRDS(data, OPTPOL_FILE)
    saveRDS(parms_NewC, PARMS_BACKUP)
    cat("File export: OK!...\n")
  } else {
    cat("Discrepancy between baseline and simulations")
  }
  
  cat("Optimal Policies: OK!...\n")
  
}





