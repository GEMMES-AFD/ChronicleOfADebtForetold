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

library("sensobol")
library("data.table")
library("ggplot2")

#Load the package
library(Rcpp)
library(rmarkdown)
library(lhs)

library(tidyverse)

#Load the algorithm
source("Source/SourceCode.R")
source("Source/utilities.R")

for (shape in 2:2){
  nvars <- 3
  ##Getting shape
  filename = paste0("Data/optpol_bastransnew2",shape,".RDS")
  data <- readRDS(filename)
  
  ##Generating LHS
  N <- 1000
  params=c("$betaen$","$v1$","$alphagw$","$alphapw$","$kappa01$","$reducXrO$")
  matrices <- c("A", "B", "AB", "BA")
  first <- total <- "azzini"
  order <- "second"
  R <- 10 ^ 3
  type <- "percent"
  conf <- 0.9
  lhs_samples <- maximinLHS(N, k = 6)
  colnames(lhs_samples)<-params
  # Define the ranges
  ranges <- matrix(c(
    0.5, 3, #betaen
    0.002543549, 0.012543549, #v1
    0.015, 0.05, #alphagw
    0.015, 0.05, #alphapw
    0.1, 0.3, #kappa01
    0,0.075 #reducXrO
  ), ncol = 2, byrow = TRUE)
  
  # Scale the Latin Hypercube Samples to the desired ranges
  mat <- t(apply(lhs_samples, 1, function(x) {
    x * (ranges[, 2] - ranges[, 1]) + ranges[, 1]
  }))
  
  ###Getting policies
  mat_pol <- data[[1]][,1:3]
  rob_samples <- 1
  sd = 0
  
  ###Building robustness calibration list
  rob_list <- lapply(1:nrow(mat_pol), function (i){
    ###Getting random policies
    matall <- matrix(ncol = ncol(mat_pol), nrow = rob_samples)
    
    for (j in 1:nrow(matall)){
      matall[j,] <- rnorm(n=nvars,mean=mat_pol[i,], sd=0)
    }
    
    
    ####Merging datasets
    matall=cbind(matall,matall[,3])
    colnames(matall)<-c("shrGrL","shrGrLFx","md_lgtr","md_bgtr")
    
    matnew<-mat
    colnames(matnew)<-c("betaen","v1","alphagw","alphapw","kappa01","reducXrO")
    
    
    matfinal <- do.call(rbind, lapply(1:nrow(matnew), function(i) {
      cbind(
        matrix(rep(matnew[i, ], each = nrow(matall)), ncol = ncol(matnew)),
        matall
      )
    }))
    
    
    
    matfinal=cbind(matfinal,rep(0.00104006,nrow(matfinal)),
                   rep(0.3,nrow(matfinal)),
                   rep(2.6,nrow(matfinal)),
                   rep(4.99377*10^-7,nrow(matfinal)),
                   rep(5,nrow(matfinal)),
                   rep(6,nrow(matfinal)),
                   rep(0.011,nrow(matfinal)),
                   rep(5,nrow(matfinal)),
                   rep(12,nrow(matfinal)),
                   rep(1,nrow(matfinal)))
    colnames(matfinal)<-c("betaen","v1","alphagw","alphapw","kappa01","reducXrO","shrGrL","shrGrLFx","md_lgtr","md_bgtr","alpha_tr","beta_tr","gamma_tr","delta_tr","lambdatr0","lambdatr1","lambdatr2","lambdatr0_adj","lambdatr1_adj","lambdatr2_adj")
    matfinal
  })
  
  robfile <-paste0("Robustness/RobcalsShape_no_new2",shape,".RDS")
  
  saveRDS(rob_list, robfile)
  rm(rob_list)
  gc()
}
