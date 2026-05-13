# ── Libraries ──────────────────────────────────────────────────────────────────
library(Rcpp)
library(tidyverse)
library(rlist)
library(NbClust)
library(ggpubr)
library(cowplot)

# ── Source ─────────────────────────────────────────────────────────────────────
source("Source/SourceCode.R")
source("Source/sourceCodeCalibration.R")
source("Source/utilities.R")

# ── Setting seed ─────────────────────────────────────────────────────────────────────
set.seed(42)



# ── Generating baseline perfs ─────────────────────────────────────────────────────────
varNames <-c('perCapita','inflation','reserves','foreignDebt','privateDebt','pubDebt','fiscalDef','hhFrag','firmsFrag','unem','CAD', 'Gip')
baseline_perfs = plain_met(c(0,0,0,0,0,0)) %>% t() %>% as.data.frame() 
colnames(baseline_perfs) <- varNames
# ── Reaction switch ─────────────────────────────────────────────────────────
reaction = TRUE

# ── File name─────────────────────────────────────────────────────────
if (reaction) {
  OPTPOL_FILE <- "MORDM_Process/MORDM_Results/optpol_REAC.RDS"
} else {
  OPTPOL_FILE <- "MORDM_Process/MORDM_Results/optpol_DS.RDS"
}

# ── Generating baseline perfs ─────────────────────────────────────────────────────────
set  <- mordm.get.set(readRDS(OPTPOL_FILE))

###Correlation plot and PCA (optional)
#ggcorrplot::ggcorrplot(cor(set[,c(5:16)]), method = "circle", outline.col = "white") + ggtitle("Correlogram across costs")
#PCA <- PCA(scale(set[,c(5:16)]), graph = FALSE)
#fviz_eig(PCA, addlabels = TRUE)
#fviz_pca_var(PCA, axes = c(1, 2), col.var = "cos2",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)

###Optimal number of clusters and best partition. This function uses severla indices to determine an optimal clustering. In case of draw, use the highest number of clusters for clarity
k.opt <-NbClust(as.matrix(scale(set[,c(7:18)])), method="kmeans", min.nc = 2, max.nc = 6)$Best.partition

set <- cbind(set, k.opt)
colnames(set)[ncol(set)] <- "Cluster"
set <- set %>% as.data.frame()
#set$Cluster <- as.factor(set$Cluster)

###Backing up to avoid redoing clustering everytime
set_DS <- set

set_reord_DS <- set_DS[,c("perCapita","reserves", "inflation", "hhFrag", "firmsFrag","unem","foreignDebt", "CAD", "fiscalDef", "privateDebt","pubDebt", "Gip" ,  "Foreign_Loans", "Foreign _Bonds", "Greenium_Loan", "Greenium_Bonds", "Interest Forgiveness", "Stock Forgiveness", "Cluster")]
baseline_perfs_DS <- baseline_perfs[,c("perCapita","reserves", "inflation", "hhFrag", "firmsFrag","unem","foreignDebt", "CAD", "fiscalDef", "privateDebt","pubDebt", "Gip")]
baseline_perfs_DS <- cbind(baseline_perfs_DS, 0,0, 0,0,0,0,4)
colnames(baseline_perfs_DS) <- c("perCapita","reserves", "inflation", "hhFrag", "firmsFrag","unem","foreignDebt", "CAD", "fiscalDef", "privateDebt","pubDebt", "Gip" ,  "Foreign_Loans", "Foreign _Bonds", "Greenium_Loan", "Greenium_Bonds", "Interest Forgiveness", "Stock Forgiveness", "Cluster")
set_reord_DS1 <- rbind(set_reord_DS, baseline_perfs_DS)



scaler <- set_reord_DS1[,1:12] %>% mutate(perCapita = perCapita/10, CAD= CAD/10, fiscalDef = fiscalDef/10)
####Generating side-by-side graph
custom_breaks <- rbind(scaler, 
                       rep(min(apply(scaler[,1:12], MARGIN = 2, FUN = min)), 12)/10, 
                       rep(max(apply(scaler[,1:12], MARGIN = 2, FUN = max)), 12)/10,
                       -rep(max(apply(scaler[,1:12], MARGIN = 2, FUN = max)), 12)/10,
                       -rep(min(apply(scaler[,1:12], MARGIN = 2, FUN = min)),12 )/10)

custom_breaks <- rbind(custom_breaks,
                       rep(min(apply(custom_breaks[,1:12], MARGIN = 2, FUN = min),12)),
                       rep(max(apply(custom_breaks[,1:12], MARGIN = 2, FUN = max)),12),
                       -rep(max(apply(custom_breaks[,1:12], MARGIN = 2, FUN = max)),12),
                       -rep(min(apply(custom_breaks[,1:12], MARGIN = 2, FUN = min), 12)))




aa <- rbind(set_reord_DS1,c( rep(0, 12),rep(0, 6), 5 ))
aa <- aa %>% mutate(Cluster=as.factor(ifelse(is.na(Cluster),max(as.numeric(levels(set_reord_noDS$Cluster)))+2,Cluster)))


#version Louis
ggparcoord_ind_yaxis(data = aa %>% mutate(perCapita = perCapita/10, CAD= CAD/10, fiscalDef = fiscalDef/10), columns = 1:12, groupColumn="Cluster", linewidth = "Cluster", custom_breaks = custom_breaks, axis_normal = c("perCapita", "foreignDebt", "CAD", "fiscalDef"), axis_normal_coeff = 10, nbreaks = 5) +
  scale_x_discrete(labels = c("Per Capita Growth", "Reserves", "Inflation", "Household Fragility", "Firm Fragility", "Unemployment", "Foreign Debt", "Current Account\nDeficit", "Fiscal Deficit", "Private Debt", "Public Debt", "Government Interest\nPayment" ,"Green Debt", "Currency Financing", "Greenium")) + 
  ggtitle("") + 
  scale_alpha(guide="none") + 
  scale_color_manual("Policy Cluster", values = c("#E64B35FF", "#00A087FF","gold","black", "#3C5488FF", "red","purple", "snow4", "blue"), labels=c("Debt Reprofiling Strategies", "Greenium-Driven External Debt Strategies","Low Exposure to External Debt Instruments", "Baseline (NDC, no policy)"), breaks = 1:4) + 
  scale_linewidth_manual(guide = "none", values = c(0.2,0.2,0.2,1,1,1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = -0.1, size =8), axis.title = element_text(size=10), axis.text.y = element_blank(), plot.title = element_text(size=10), legend.text = element_text(size = 8), plot.subtitle = element_text(size=10) ,legend.title=element_text(size = 10),legend.position = "bottom") +
  xlab("") + ylab("Cost (% SOS Size)") +
  annotate("segment", x = 12.075, y = 0.55, xend = 12.075, yend = 1, arrow=arrow(length = unit(0.1, 'cm'))) +
  annotate("text", x = 12.45, y = 0.75,label="Outside SOS", angle=270, size = 3.5) +
  annotate("segment", x = 12.075, y = 0.45, xend = 12.075, yend = 0, arrow=arrow(length = unit(0.1, 'cm'))) +
  annotate("text", x = 12.45, y = 0.25,label="Inside SOS", angle=270, size = 3.5)

filename <- paste0("Figures/Greendebtperf_DS2.png", shape,".png")
ggsave(filename, dpi = 400, width = 10, height = 5)
