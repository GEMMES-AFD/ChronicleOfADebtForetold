library(OpenMORDM)
library(lhs)
library(rdyncall)
library(parallel)
library(foreach)
library(future)
library(future.apply)
library(Rcpp)
library(rmarkdown)
library(tidyverse)
library(readxl)
library(xtable)
library(rlist)
library("sensobol")
library("data.table")
library("ggplot2")
library(viridis)
library(htmlwidgets)
library(webshot2)
library(factoextra)
library(ggforce)
library(dtwclust)
library(corrplot)
library(FactoMineR)
library(kernlab)

options(warn = -1)

#Load the algorithm
source("Source/SourceCode.R")
source("Source/sourceCodeCalibration.R")
source("Source/utilities.R")
source("Extrafunctions.R")

###Overarching parameters
reducXro = 0.025
shapes = 2


##Model Files
MORDM_FILE = "model_equations_MORDM_Unified.R"

###MORDM settings
n_sim <- 1000
nvars = 6 ##Number of policies
nobj <- 12 ##Number of constraints
nconstrs <- 0

LHS_N <- 1000
PARMS_LHC <- c("$betaen$","$v1$","$alphagw$","$alphapw$","$zetafx3$","$rho$","$reducXrO$")
POLICY_NAMES <- c("shrGrLFx", "shrGrBw", "md_lgtr","md_bgtr","decds","mdds")

##Simulations
n_cores <- 12

###Setting seed
set.seed(42)

###Calling MORDM functions
source("MORDM_Process/utils/MORDM_Functions.R")

###File names
OPTPOL_FILE <- "MORDM_Process/MORDM_Results/optpol_DS_REAC.RDS"
PARMS_BACKUP <- "MORDM_Process/MORDM_Results/parms_NewC.RDS"
ROBCAL_FILE <- "MORDM_Process/MORDM_Results/RobcalsShape_REAC.RDS"
ROBLIST_FILE <- "MORDM_Process/MORDM_Results/roblist_REAC.RDS"
REG1_FILE <- "MORDM_PROCESS/MORDM_Results/Regret1_REAC.RDS"
REG2_FILE <- "MORDM_PROCESS/MORDM_Results/Regret2_REAC.RDS"
SAT1_FILE <- "MORDM_PROCESS/MORDM_Results/Sat1_REAC.RDS"
SAT2_FILE <- "MORDM_PROCESS/MORDM_Results/Sat2_REAC.RDS"
ROBPOL_FILE <- "MORDM_PROCESS/MORDM_Results/rob_pol.RDS"

###Consistency Switch
PARALLELVANILLATEST <- T



###Regret metric switches
COMPUTE_REG1 <- F
COMPUTE_REG2 <- T
COMPUTE_SAT1 <- F
COMPUTE_SAT2 <- F

###Progress bar (set TRUE to show txtProgressBar / progressr bars in step 05)
SHOW_PROGRESS <- TRUE


