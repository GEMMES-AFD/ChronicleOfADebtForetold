#SOS
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
  mutate(center = ifelse(direction == 0, 0.5*(L_Bound+H_Bound), 0))

noutcomes = nrow(threshold)
col_out <- (nvars+1):(nvars+noutcomes)

if (nobj != nrow(threshold)){
  stop("Mismatch between number of objectives and SOS size")
}
