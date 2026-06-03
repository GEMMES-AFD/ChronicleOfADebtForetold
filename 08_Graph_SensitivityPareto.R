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

library("sensobol")
library("data.table")
library("ggplot2")

#Load the package
library(Rcpp)
library(rmarkdown)
library(lhs)

library(tidyverse)
library(dtwclust)
library(NbClust)
library(paletteer)


#Load the algorithm
source("Source/SourceCode.R")
source("Source/utilities.R")
library(corrplot)
library(factoextra)
library(FactoMineR)
library(kernlab)
library(ggridges)
library(eaf)
library(xtable)

labels_altsow <- c("Low"="Sloppy",
                   "High"="Stiff")

for (altSOws in c("Low", 'High')){


varNames <-c('perCapita','inflation','reserves','foreignDebt','privateDebt','pubDebt','fiscalDef','hhFrag','firmsFrag','unem','CAD', 'Gip')
shape = 2
plot_list <- list()


polfile <- paste0("MORDM_Process/MORDM_Results/optpol_DS.RDS")
set  <- mordm.get.set(readRDS(polfile)) %>% as.data.frame() %>%
  dplyr::mutate(sens = 0)
set_ref  <- mordm.get.set(readRDS(polfile)) %>% as.data.frame() 
set_ref <- set_ref %>%
  dplyr::mutate(across(1:ncol(set_ref), ~(.x - min(.x))/(max(.x)-min(.x))))


sets <- list()
sets[[1]] <- set_ref


dist_igd_out <- c()
dist_igd_pol <- c()

dist_haus_out <- c()
dist_haus_pol <- c()

for (kk in 1:5){
    ####Treating baseline
    if (altSOws == "High"){
    polfile <- paste0("Data/optpol_bastransDSHigh",shape,"_robustness_",kk,".RDS")
    } else {
      polfile <- paste0("Data/optpol_bastransDSLow",shape,"_robustness_",kk,".RDS")
    }
    set2  <- mordm.get.set(readRDS(polfile)) %>% as.data.frame()%>%
      dplyr::mutate(sens = kk)
    set <- rbind(set, set2)
    
    set2_ref <- mordm.get.set(readRDS(polfile)) %>% as.data.frame()
    set2_ref <- set2_ref  %>%
      dplyr::mutate(across(1:ncol(set2_ref), ~(.x - min(.x))/(max(.x)-min(.x))))
    
    sets[[kk+1]] <- set2_ref

if (kk >0){
  dist_igd_pol[kk] <- igd_plus(as.matrix(set_ref[,1:6]),
                         as.matrix(set2_ref[,1:6]))
  dist_haus_pol[kk] <- avg_hausdorff_dist(as.matrix(set_ref[,1:6]),
                               as.matrix(set2_ref[,1:6]))
  
  dist_igd_out[kk] <- igd_plus(as.matrix(set_ref[,7:18]),
                               as.matrix(set2_ref[,7:18]))
  dist_haus_out[kk] <- avg_hausdorff_dist(as.matrix(set_ref[,7:18]),
                                          as.matrix(set2_ref[,7:18]))
}  
}

###Saving sets for further uses
saveRDS(set, paste0("Set",altSOws, ".RDS"))

####Distance table
name_column <- c("IGD - Policies", "Hausdorff - Policy", "IGD - Outcome", "Hausdorff - Policies")
distance_matrix <- rbind(round(dist_igd_pol,2), 
                         round(dist_haus_pol,2),
                               round(dist_igd_out,2),
                                     round(dist_haus_out,2))
distance_matrix <- cbind(name_column, distance_matrix) %>% as.data.frame()
colnames(distance_matrix) <- c("Distance Metric","Alternative 1", "Alternative 2", "Alternative 3", "Alternative 4", "Alternative 5")
print(distance_matrix)

sink(paste0("Tables/LatexTable_",altSOws,".tex"))
print(xtable(distance_matrix, paste0("distance_matrix_",altSOws,".tex"), caption=paste0("Distance from reference Pareto front: ", labels_altsow[altSOws])))
sink()
saveRDS(distance_matrix, paste0("Distance_",altSOws,".RDS"))
}


###Calling sets
set_low <- readRDS("Setlow.RDS")  %>% dplyr::mutate(Params = ifelse(sens ==0, "Reference", "Sloppy"))
set_high <- readRDS("Sethigh.RDS") %>% filter(sens > 0) %>% dplyr::mutate(Params = "Stiff")

big_set <- rbind(set_low,
                 set_high)


#Version mo
sens_labels <- c(
   "Reference" = 0,
   "Alternative 1"=  1,
   "Alternative 2" = 2,
   "Alternative 3" = 3 ,
  "Alternative 4" = 4,
  "Alternative 5" = 5
)



plot_title <- paste0("Illustrative arbitrage\nacross four dimensions")

big_set$sens <- factor(big_set$sens, 
                   levels = c(0, 1, 2, 3, 4, 5), 
                   labels = c("Reference", "Alternative 1", "Alternative 2", "Alternative 3", "Alternative 4", "Alternative 5"))

plot <- ggplot(data = big_set, 
               aes(x = perCapita, 
                   y = pubDebt, 
                   fill = inflation, 
                   size = CAD, 
                   shape = Params)) +
  geom_point() +
  facet_wrap(~as.factor(sens), labeller = as_labeller(sens_labels)) +
  scale_fill_viridis("Cost on\nInflation", begin = 0.5, option = "turbo") +
  scale_size("Cost on\nCurrent\nAccount\nDeficit") +
  scale_shape_manual(values=c(24, 22, 21)) +
  xlab("Cost on Growth") +
  ylab("Cost on Public Debt") +
  ggtitle(plot_title) +
  theme(
    axis.text.x = element_text(hjust = 1, size = 15),
    axis.text.y = element_text(hjust = 1, size = 15),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 25),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 20)
  )
plot

ggsave(paste0("Figures/ComparativePFront.png"), dpi = 400, width = 7.5, height = 5)

plot <- ggplot(data = big_set, 
               aes(x = reserves , 
                   y = foreignDebt , 
                   fill = privateDebt    , 
                   size = pubDebt,
                   shape = Params )) +
  geom_point() +
  facet_wrap(~as.factor(sens), labeller = as_labeller(sens_labels)) +
  scale_fill_viridis("Cost on\nPublic\nDebt", begin = 0.5, option = "turbo") +
  scale_size("Cost on\nPublic\nDebt") +
  scale_shape_manual(values=c(24, 22, 21)) +
  xlab("Cost Reserve Ratio") +
  ylab("Cost Foreign Debt") +
  ggtitle(plot_title) +
  theme(
    axis.text.x = element_text(hjust = 1, size = 15),
    axis.text.y = element_text(hjust = 1, size = 15),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 25),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 20)
  )
plot

ggsave(paste0("Figures/ComparativePFront1.png"), dpi = 400, width = 7.5, height = 5)

plot <- ggplot(data = big_set, 
               aes(x = fiscalDef     , 
                   y = hhFrag  , 
                   fill = unem           , 
                   size = Gip,
                   shape = Params         )) +
  geom_point() +
  facet_wrap(~as.factor(sens), labeller = as_labeller(sens_labels)) +
  scale_fill_viridis("Cost on\nGov.Int.Pay", begin = 0.5, option = "turbo") +
  scale_shape_manual(values=c(24, 22, 21)) +
  scale_size("Cost on\nGov.Int.Pay") +
  xlab("Cost on Fiscal Deficit") +
  ylab("Cost on Household Deficit") +
  ggtitle(plot_title) +
  theme(
    axis.text.x = element_text(hjust = 1, size = 15),
    axis.text.y = element_text(hjust = 1, size = 15),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 25),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 20)
  )
plot

ggsave(paste0("Figures/ComparativePFront2.png"), dpi = 400, width = 7.5, height = 5)

plot_title <- paste0("Illustrative arbitrage across\nfour dimensions (Sloppy)")
plot <- ggplot(data = big_set %>% filter(Params == "Sloppy" | Params == "Reference"), 
               aes(x = fiscalDef     , 
                   y = hhFrag  , 
                   fill = unem           , 
                   size = Gip,
                   shape = Params         )) +
  geom_point() +
  facet_wrap(~as.factor(sens), labeller = as_labeller(sens_labels)) +
  scale_fill_viridis("Cost on\nGov.Int.Pay", begin = 0.5, option = "turbo") +
  scale_shape_manual(values=c(24, 22, 21)) +
  scale_size("Cost on\nGov.Int.Pay") +
  xlab("Cost on Fiscal Deficit") +
  ylab("Cost on Household Deficit") +
  ggtitle(plot_title) +
  theme(
    axis.text.x = element_text(hjust = 1, size = 15),
    axis.text.y = element_text(hjust = 1, size = 15),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 25),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 20)
  )
plot

ggsave(paste0("Figures/ComparativePFront_Sloppy.png"), dpi = 400, width = 7.5, height = 5)


plot_title <- paste0("Illustrative arbitrage across\nfour dimensions (Stiff)")
plot <- ggplot(data = big_set %>% filter(Params == "Stiff" | Params == "Reference"), 
               aes(x = fiscalDef     , 
                   y = hhFrag  , 
                   fill = unem           , 
                   size = Gip,
                   shape = Params         )) +
  geom_point() +
  facet_wrap(~as.factor(sens), labeller = as_labeller(sens_labels)) +
  scale_fill_viridis("Cost on\nGov.Int.Pay", begin = 0.5, option = "turbo") +
  scale_shape_manual(values=c(24, 22, 21)) +
  scale_size("Cost on\nGov.Int.Pay") +
  xlab("Cost on Fiscal Deficit") +
  ylab("Cost on Household Deficit") +
  ggtitle(plot_title) +
  theme(
    axis.text.x = element_text(hjust = 1, size = 15),
    axis.text.y = element_text(hjust = 1, size = 15),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 25),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 20)
  )
plot

ggsave(paste0("Figures/ComparativePFront_Stiff.png"), dpi = 400, width = 7.5, height = 5)



