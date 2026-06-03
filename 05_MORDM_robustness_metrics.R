library(OpenMORDM)
library(lhs)
library(rdyncall)
library(parallel)
library(foreach)
library(future.apply)
library(viridis)
library(htmlwidgets)
library(webshot2)
library(factoextra)
library(ggforce)
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

library(tidyverse)
library(dtwclust)
library(future)



#Load the algorithm
source("Source/SourceCode.R")
source("Source/utilities.R")
library(corrplot)
library(factoextra)
library(FactoMineR)
library(kernlab)


for (shape in 2:2){
  
  
  ###Treating Robustness
  robfile <- paste0("Robustness/roblist_no_new2",shape,".RDS")
  roblist <- readRDS(robfile)
  maxx <- lapply(roblist, function(x){
    y <- x[-which(x[,'Failure'] == T),7:18]
    maxx <- apply(y, 2, max, na.rm=T)*1.01
  }) %>% list.rbind()
  
  maxx <- apply(maxx, 2, max, na.rm=T)
  
  rob_pol <- c()
  
  roblist2 <- lapply(roblist, function(x){
    for (i in which(x[,'Failure'] == T)){
      x[i,7:18] <- maxx
    }
    x
  })
  
  
  plan(cluster)
  polfile <- paste0("Data/optpol_bastransnew2",shape,".RDS")
  set  <- mordm.get.set(readRDS(polfile))
  
  
  regret_type1 <- future.apply::future_lapply(1:length(roblist2), function(ll) {
    a<-roblist2[[ll]] %>% as.data.frame()
    max(sapply(1:ncol(a[,7:18]), function(j){
      quantile(sapply(1:nrow(a), function(k) {
        (a[k,j+3] - set[ll,j+3])/abs(set[ll,j+3])
      }), probs =0.9, na.rm=T)
    }))}
  ) %>% list.rbind()
  
  
  regret_type1 <- cbind(set, regret_type1)
  colnames(regret_type1)[ncol(regret_type1)] <- "Regret"
  
  rob_pol[1] <- which.min(regret_type1[,"Regret"])
  
  
  start_time <- Sys.time()
  regret_type2 <- future.apply::future_lapply(1:length(roblist2), function(ll){max(
    sapply(7:18, function(j) {
      quantile(sapply(1:nrow(roblist2[[ll]]), function(i){
        
        best <-min(sapply(1:length(roblist2), function(k){
          roblist2[[k]][i,j]
        }), na.rm=T)
        
        (roblist2[[ll]][i,j]- best)/abs(best)
      })
      , 0.9, na.rm=T)})
  )}
  ) %>% list.rbind() 
  end_time <- Sys.time()
  end_time - start_time
  
  
  
  regret_type2 <- cbind(set,regret_type2)
  colnames(regret_type2)[ncol(regret_type2)] <- "Regret"
  
  rob_pol[2] <-which.min(regret_type2[,"Regret"])
  
  roblist2_df <- list.rbind(roblist2)
  
  quantile_rob <- apply(roblist2_df[,5:16], 2, function(x) quantile(x, 0.15))
  
  
  
  
  satisficing_type1 <- lapply(1:length(roblist2), function(ll){
    a <- roblist2[[ll]][,7:18]
    b<- apply(a, 2, function(x, quantile_rob) x < quantile_rob, quantile_rob = quantile_rob)
    d <- apply(b, 1, sum, na.rm=T)
    dd <- length(d[which(d==11)])
    dd/length(d)
  }
  ) %>% list.rbind()
  
  satisficing_type1 <- cbind(set, satisficing_type1)
  colnames(satisficing_type1)[ncol(satisficing_type1)] <- "Satisficing"
  
  rob_pol[3] <-which.max(satisficing_type1[,"Satisficing"])
  
  calfile <- paste0("Robustness/RobcalsShape_no_new2",shape,".RDS")
  params <- readRDS(calfile)[[1]]
  
  
  
  dist <- apply(params[,1:6], 1, function(x){
    refs <- c(parms_NewC[c("betaen","v1", "alphagw", "alphapw", "kappa01")],0.025)
    for (i in 1:6){
      x[i] <- (x[i] - min(params[,i]))/(max(params[,i]) - min(params[,i]))
      refs[i] <- (refs[i] - min(params[,i]))/(max(params[,i]) - min(params[,i]))
    }
    norm(x - refs, "2")
  })
  
  satisficing_type2 <- lapply(1:length(roblist2), function(ll){
    a <- roblist2[[ll]][,7:18]
    b<- apply(a, 2, function(x, quantile_rob) x > quantile_rob, quantile_rob = quantile_rob)
    b <- apply(b, 1, sum, na.rm=T)
    f <- min(dist[b>0], na.rm=T)
    f[which(is.infinite(f))] <- NA  
    f
  }) %>% list.rbind()
  
  
  
  satisficing_type2 <- cbind(set, satisficing_type2)
  colnames(satisficing_type2)[ncol(satisficing_type2)] <- "Satisficing"
  
  rob_pol[4] <-which.max(satisficing_type2[,"Satisficing"])
  
  
  
  saveRDS(regret_type1, "Robustness/reg1_threshresnewGD2.RDS")
  saveRDS(regret_type2, "Robustness/reg2_threshresnewGD2.RDS")
  saveRDS(satisficing_type1, "Robustness/sat1_threshresnewGD2.RDS")
  saveRDS(satisficing_type2, "Robustness/sat2_threshresnewGD2.RDS")
  
  
}
