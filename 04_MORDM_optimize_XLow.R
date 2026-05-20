library(OpenMORDM)
library(lhs)
library(rdyncall)
library(parallel)

options(warn = -1) 
library(readxl)
library(xtable)
library(rlist)
library("sensobol")
library("data.table")
library("ggplot2")
library(Rcpp)
library(rmarkdown)
library(lhs)
library(readr)
library(tidyverse)

source("Source/SourceCode.R")
source("Source/utilities.R")
source("Extrafunctions.R")

###Calling parameters -- Five best
alt_bas <- read_delim("Data/XLow_ord.csv", delim = ";", 
                      escape_double = FALSE, trim_ws = TRUE) %>%
  slice(1:5)

###Generating System
xr0 = 0.025
event1 <- list(triggerDate=4, reducXrO=xr0)
SOEM <- cppMakeSys(
  fileName = "model_equations_MORDM_Unified.R",
  reportVars = 3, 
  eventTime = list(event1)
)

varNames <- c('perCapita','inflation','reserves','foreignDebt','privateDebt',
              'pubDebt','fiscalDef','hhFrag','firmsFrag','unem','CAD','Gip','GDP')

###Creating SOS based on baseline performances
threshold <- as.data.frame(matrix(nrow=12, ncol=4))
colnames(threshold) <- c("Variables", "L_Bound", "H_Bound", "direction")
threshold[,1] <- c("growth","inflation","reserves","foreignDebt","privateDebt",
                   "pubDebt","fiscalDef","hhFrag","firmsFrag","unem","CAD","Gip")
threshold[1,2:4]  <- c(0, 0.2,  0)
threshold[2,2:4]  <- c(0.01, 0.04, 0)
threshold[3,2:4]  <- c(0, 0.18, -1)
threshold[4,2:4]  <- c(0, 0.4,  1)
threshold[5,2:4]  <- c(0, 0.9,  1)
threshold[6,2:4]  <- c(0, 0.63, 1)
threshold[7,2:4]  <- c(0, 0.035,1)
threshold[8,2:4]  <- c(0, 0.45, 1)
threshold[9,2:4]  <- c(0, 0.26, 1)
threshold[10,2:4] <- c(-0.01, 0.12, 0)
threshold[11,2:4] <- c(0, 0.05, 1)
threshold[12,2:4] <- c(0, 0.04, 1)
threshold <- threshold %>%
  mutate(center = ifelse(direction == 0, 0.5*(L_Bound+H_Bound), 0))

# ── Detect cores ──────────────────────────────────────────────────────────────
n_cores <- min(nrow(alt_bas), detectCores() - 1)
cat("Running on", n_cores, "cores\n")

# ── Worker function for one kk ────────────────────────────────────────────────
run_one_kk <- function(kk) {
  
  # Each worker needs its own copies of shared objects
  # (mclapply forks the process so these are already available,
  #  but we re-source to be safe in case of copy-on-write issues)
  source("Source/SourceCode.R",         local = TRUE)
  source("Source/sourceCodeCalibration.R", local = TRUE)
  source("Source/utilities.R",          local = TRUE)
  source("Extrafunctions.R",            local = TRUE)
  
  parms_worker <- SOEM$parms  # local copy — avoids race conditions
  
  for (col in colnames(alt_bas)) {
    parms_worker[col] <- alt_bas[[col]][kk]
  }
  
  # NDC parameters
  parms_worker['lambdatr0']     <- 5
  parms_worker['lambdatr1']     <- 6
  parms_worker['lambdatr0_adj'] <- 5
  parms_worker['lambdatr1_adj'] <- 12
  parms_worker['lambdatr2_adj'] <- 1
  
  # ── Inner loop over shape (sequential within each worker) ──────────────────
  for (shape in 2:2) {
    
    if (shape == 1) {
      parms_worker['alpha_tr'] <- 0.000709889
      parms_worker['beta_tr']  <- 0.09
      parms_worker['gamma_tr'] <- 1.6
      parms_worker['delta_tr'] <- 4.99377e-7
      event_shape <- list(triggerDate=4, reducXrO=0.025)
    } else if (shape == 2) {
      parms_worker['alpha_tr'] <- 0.00104006
      parms_worker['beta_tr']  <- 0.3
      parms_worker['gamma_tr'] <- 2.6
      parms_worker['delta_tr'] <- 4.99377e-7
      event_shape <- list(triggerDate=4, reducXrO=0.025)
    } else if (shape == 3) {
      parms_worker['alpha_tr'] <- 0.00075
      parms_worker['beta_tr']  <- 0.18
      parms_worker['gamma_tr'] <- 2.11
      parms_worker['delta_tr'] <- 4.99377e-7
      event_shape <- list(triggerDate=4, reducXrO=0.025)
    } else if (shape == 4) {
      parms_worker['alpha_tr'] <- 0.09
      parms_worker['beta_tr']  <- 0.6
      parms_worker['gamma_tr'] <- 0.9
      parms_worker['delta_tr'] <- 4.99377e-7
      event_shape <- list(triggerDate=4, reducXrO=0.025)
    }
    
    # ── Objective function (closed over parms_worker) ─────────────────────────
    my_run <- function(levers) {
      
      costs <- matrix(ncol = nrow(threshold))
      
      if (length(levers) != 6) stop("Expected exactly 6 decision variables.")
      
      parms_NewC <- parms_worker  # local copy inside each evaluation
      parms_NewC['dsactive']  <- 1
      parms_NewC['shrGrL']    <- 1
      parms_NewC['shrGrLFx']  <- levers[1]
      parms_NewC['shrGrBw']   <- levers[2]
      parms_NewC['md_lgtr']   <- levers[3]
      parms_NewC['md_bgtr']   <- levers[4]
      parms_NewC['mdds']      <- levers[5]
      parms_NewC['decds']     <- levers[6]
      
      newScen <- cppRK4(SOEM, parms = parms_NewC,
                        times = seq(from=2019, to=2050, by=0.1),
                        eventTime = list(event_shape))
      
      newScen <- newScen[, varNames] %>% as.data.frame() %>%
        mutate(perCapita = 100 * diff(perCapita) / perCapita) %>%
        mutate(Gip = Gip / GDP) %>%
        dplyr::select(-GDP) %>%
        filter(row_number() < 111 & row_number() > 10)
      
      failure <- max(newScen$reserves < 0) > 0
      crisistime <- if (failure) {
        min(which(newScen$reserves < 0), 
            min(which(is.na(newScen$reserves))))
      } else NA
      
      for (m in 1:nrow(threshold)) {
        variable  <- newScen[, m]
        lbound    <- threshold[m, 2]
        hbound    <- threshold[m, 3]
        direction <- threshold[m, 4]
        center    <- threshold[m, 5]
        
        up           <- ifelse(variable > center,  1, 0)
        above_higher <- ifelse(variable >= hbound, 1, 0)
        below_lower  <- ifelse(variable <= lbound, 1, 0)
        
        timeperformance <- if (direction == 0) {
          ifelse(up > 0,
                 ifelse(above_higher == 1,
                        100*abs((variable - hbound))/(hbound-lbound),
                        -sqrt(100*abs(hbound - variable)/(hbound-lbound))),
                 ifelse(below_lower == 1,
                        100*abs((variable - lbound))/(hbound-lbound),
                        -sqrt(100*abs(variable - lbound)/(hbound-lbound))))
        } else if (direction == -1) {
          ifelse(above_higher == 0,
                 100*abs((variable - hbound)),
                 -sqrt(100*abs(variable - hbound)/(hbound-lbound)))
        } else {
          ifelse(above_higher == 1,
                 100*abs((variable - hbound))/(hbound-lbound),
                 -sqrt(100*abs((hbound - variable))/(hbound-lbound)))
        }
        
        timeperformance <- timeperformance[!is.infinite(timeperformance)]
        
        timeperformance <- sapply(timeperformance, function(tp) {
          if (failure) 100^2 * nrow(newScen) / crisistime
          else min(tp, 100^2)
        })
        
        costs[, m] <- sum(timeperformance, na.rm=TRUE) / 
          length(which(!is.na(timeperformance)))
      }
      
      return(as.numeric(costs))
    }
    
    # ── Borg optimisation ─────────────────────────────────────────────────────
    bounds  <- matrix(rep(c(0, 1), 6), nrow=2)
    problem <- define.problem(
      my_run, nobj=12, nvars=6, nconstrs=0,
      bounds   = bounds,
      epsilons = rep(1e-01, 12),
      names    = c("Foreign_Loans", "Foreign_Bonds", "Greenium_Loan",
                   "Greenium_Bonds", "Interest_Forgiveness", "Stock_Forgiveness",
                   varNames[-length(varNames)])
    )
    
    data <- borg.optimize(problem, 10000, verbose=TRUE)
    colnames(data[[1]]) <- c("Foreign_Loans", "Foreign_Bonds", "Greenium_Loan",
                             "Greenium_Bonds", "Interest_Forgiveness",
                             "Stock_Forgiveness", varNames[-length(varNames)])
    
    filename <- paste0("Data/optpol_bastransDS", shape, "_robustness_", kk, ".RDS")
    saveRDS(data, filename)
    cat("Saved:", filename, "\n")
  }
  
  return(invisible(kk))
}

# ── Launch parallel runs over kk ──────────────────────────────────────────────
results <- mclapply(
  1:nrow(alt_bas),
  run_one_kk,
  mc.cores    = n_cores,
  mc.set.seed = TRUE  
)

cat("All runs complete.\n")