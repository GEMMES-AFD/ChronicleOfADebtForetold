#Load the package
library(Rcpp)
library(gridExtra)
library(rlist)
library(tidyverse)
library(xtable)


#Load the algorithm
source("Source/SourceCode.R")
source("Source/sourceCodeCalibration.R")
source("Source/utilities.R")

event1 <- list(triggerDate=4, reducXrO="0.025")
#,reducXbiO="0")

SOEM1 <- cppMakeSys(fileName = "model_equations_MORDM.R",reportVars=3, eventTime=list(event1))

#Baseline

parms1 = SOEM1$parms
parms1['lambdatr0'] = 5                 #Speed of the NDC investment path
parms1['lambdatr1']  = 6                #Initial period of the NDC investment path
parms1['lambdatr2']  = 0.011            #Target NDC investment as a share of NFC's capital stock in 2019
parms1['lambdatr0_adj'] = 5             #Speed of the NDC investment path
parms1['lambdatr1_adj']  = 12           #Initial period of the NDC investment path
parms1['lambdatr2_adj']  = 1            #Target NDC investment as a share of NFC's capital stock in 2019
#Louis Daumas_boxing
parms1['alpha_tr']=0.00104006
parms1['beta_tr']=0.3
parms1['gamma_tr']=2.6
parms1['delta_tr']=4.99377*10^-7
parms1['betadebtSwapFXLgFX']=4
parms1['betadebtSwapFXBgFX']=4
parms1['betadebtSwapFXLgFXtr']=4
#NDC conventional
parms1['shrGrL']= 1                   
parms1['shrGrLFx'] = 0.5

res_list <- list()

baseline<-cppRK4(SOEM1, eventTime=list(event1),parms = parms1)%>% as.data.frame() %>% mutate(Scenario = "Baseline") %>% filter(time > 2022.9 &time < 2030)
res_list[[1]]<- baseline

#1 loans with greenium

parms4=parms1
parms4['dsactive']=0
parms4['shrGrLFx']=1
parms4['md_lgtr'] = 0.8

cas1<-cppRK4(SOEM1, eventTime = list(event1),parms = parms4) %>% as.data.frame() %>% mutate(Scenario = "DR1\nGreenium on Loans") %>% filter(time > 2022.9 &time < 2030)
res_list[[2]]<- cas1

#2 Original sin

parms5=parms1
parms5['dsactive']=0
parms5['shrGrLFx'] = 0
parms5['shrGrBw']= 1
parms5['md_bgtr'] = 0

cas2<-cppRK4(SOEM1, eventTime = list(event1),parms = parms5) %>% as.data.frame() %>% mutate(Scenario = "DR2\nLCY Bonds") %>% filter(time > 2022.9 &time < 2030)
res_list[[3]] <- cas2
#3 interest rate renegociation on old debt

parms2=parms1
parms2['dsactive']=1
parms2['mdds']=0.8
parms2['decds']=0

cas3<-cppRK4(SOEM1, eventTime=list(event1),parms = parms2) %>% as.data.frame() %>% mutate(Scenario = "DR3\nInterest Rate Renegotiation")%>% filter(time > 2022.9 &time < 2030)
res_list[[4]] <- cas3
#4 Principal reduction

parms3=parms1
parms3['dsactive']=1
parms3['mdds']=0
parms3['decds']=0.8

cas4<-cppRK4(SOEM1, eventTime=list(event1),parms = parms3)  %>% as.data.frame() %>% mutate(Scenario = "DR4\nPrincipal Adjustment")%>% filter(time > 2022.9 &time < 2030)
res_list[[5]] <- cas4

#5 fiscal rule 
parms6=parms1
parms6['FRactive']=1
parms6['gammafr']=100
parms6['alphafr']=0.25
parms6['dfr']=0.6
parms6['fdr']=0.03


cas5<-cppRK4(SOEM1, eventTime = list(event1),parms = parms6) %>% as.data.frame() %>% mutate(Scenario = "FR\n Fiscal rule") %>% filter(time > 2022.9 &time < 2030)
res_list[[6]]<- cas5



diagvars <- c("GipGDP", "en", "premgd", "rsk", "FIP", "perCapita","unem","inflation","CAD","reserves",
              "hhFrag",
              #"ratFFX")
              "GDP")
diagvar_labs <- c(GipGDP = "Public Interest Payment\n(% GDP) (a)", en="Nominal Exchange \nRate (c)", premgd ="Premium on public Debt\n(Pct) (i)", rsk = "Country Risk (Pct) (k)", FIP = "Foreign International\nPosition\n(% GDP) (d)", perCapita = "GDP per capita (USD) (h)",unem = "Unemployment (l)", inflation = "Inflation (g)", CAD = "Current account deficit (b)", reserves = "Foreign reserves\n(% GDP)(j)", hhFrag = "Households fragility (f)", #ratFFX = "Credit rationing")
                   GDP ="GDP (e)" )


res_df <- res_list %>% list.rbind() %>% mutate(Scenario = as.factor(Scenario)) %>% mutate(ehat = 100*endot/en, 
                                                                                          GipGDP = 100*Gip/GDP, 
                                                                                          Debtcost = 100*(GipNoDS/Gip)/(iktr*pktr),
                                                                                          totresgdp = (Rfx*en + Dfxw*en)/GDP,
                                                                                          curdemgdp =  Dfxw*en, 
                                                                                          Gipextgdp = 100*GipExt/GDP, 
                                                                                          GipextnoDSGDP = 100*GipExtNoDS/GDP,
                                                                                          DebtCostExt = 100*(GipExtNoDS-GipExt)/(iktr*pktr),
                                                                                          Iktr=iktr*pktr,
                                                                                          debtSwapgdp = 100*debtSwapFXLgFX/GDP,
                                                                                          foreignDebtGDP = 100*(ibgfx*Bgfx*en + ilgfx*Lgfx*en + Lgfxtr*ilgfxtr*en)/GDP,
                                                                                          Foreignbondsgdp = 100*(Bgfx*en)/GDP,
                                                                                          ForeignLoansGDP = 100*(Lgfx*en)/GDP,
                                                                                          GovernmentDeficitGDP = 100*(TFNG)/GDP,
                                                                                          TotalLiabilities  = FIP*GDP/en,
                                                                                          CADGDP = 100*(X-IM)/GDP,
                                                                                          en = en*4300/1.271753,
                                                                                          rsk=100*rsk,
                                                                                          FIP = 100*FIP, 
                                                                                          premgd = 100*premgd

) %>% dplyr::select(Scenario,time, all_of(diagvars)) %>% pivot_longer(3:(length(diagvars)+2), names_to = "Variable", values_to = "Value") %>% mutate(Variable = as.factor(Variable)) %>% mutate(Variable = relevel(Variable, ref = "GipGDP"))

res_df$Scenario <- factor(res_df$Scenario, levels = c("Baseline", "DR1\nGreenium on Loans", "DR2\nLCY Bonds", "DR3\nInterest Rate Renegotiation", "DR4\nPrincipal Adjustment","FR\n Fiscal rule"))

# Plot
ggplot(res_df, aes(x = time, y = Value, 
                   group = Scenario, 
                   color = Scenario, 
                   linetype = Scenario)) + 
  geom_line(linewidth = 0.7) + 
  facet_wrap(~as.factor(Variable), 
             scales = "free_y", 
             labeller = as_labeller(diagvar_labs), 
             ncol = 3) +   # <- controls number of columns (facets per row)
  scale_colour_manual(values = c("black", "#2E5F8DFF", "#F8A049FF", "#D82632FF", "#A50021FF","darkgrey")) + 
  theme_classic() + 
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank(), 
        strip.text = element_text(colour = 'black', size = 11), 
        strip.background = element_blank(), 
        axis.text = element_text(size = 10), 
        legend.text = element_text(size = 10), 
        legend.title = element_blank(),
        legend.box = "vertical",
        legend.spacing.y = unit(1,"cm"),
        legend.key.height =unit(1,"cm")) + 
  
  ylab("") + 
  xlab("Year")


ggsave("Images/clean_plot_res18.png", width = 10, height = 5, dpi = 400)




# Indicators

res_list_named <- list(
  "Baseline" = baseline,
  "DR1\nGreenium on Loans" = cas1,
  "DR2\nLCY Bonds" = cas2,
  "DR3\nInterest Rate Renegotiation" = cas3,
  "DR4\nPrincipal Adjustment" = cas4,
  "FR\n Fiscal rule" = cas5
)

# Filter time and annualize by summing (quarterly -> annual)
res_list_ann <- lapply(res_list_named, function(x){
  x %>% 
    filter(time > 2022.9 & time < 2030) %>%
    mutate(year = trunc(time)) %>% 
    group_by(year) %>%
    select(where(is.numeric)) %>%
    summarise_all(sum)
})

compute_indicators <- function(df) {
  avg_annGIP     <- mean(df$Gip, na.rm = TRUE)
  total_GIP      <- sum(df$Gip, na.rm = TRUE)
  total_GIP_GDP  <- sum(df$Gip/df$GDP, na.rm = TRUE)
  tot_green_investment <- sum(df$Iktr, na.rm = TRUE)
  avg_trade_deficit    <- mean(100*(df$X - df$IM)/df$GDP, na.rm = TRUE)
  currency_variation   <- mean(as.numeric(df$endot)/as.numeric(df$en))
  avg_Interest_fxdebt  <- mean(df$ibgfx + df$ilgfx + df$ilgfxtr + df$ibgtr, na.rm = TRUE)
  
  return(data.frame(
    Total_Green_Investment = tot_green_investment,
    Average_Debt_Cost      = avg_annGIP,
    Total_Debt_Cost        = total_GIP,
    Total_Debt_Cost_GDP    = total_GIP_GDP,
    Average_Trade_Deficit  = avg_trade_deficit,
    Currency_Variation     = currency_variation,
    Average_interest_fxdebt = avg_Interest_fxdebt
  ))
}

results <- lapply(res_list_ann, compute_indicators)
results_df <- do.call(rbind, results)
rownames(results_df) <- c("Baseline", "DR1 Greenium on Loans", "DR2 LCY Bonds", 
                          "DR3 Interest Rate Renegotiation", "DR4 Principal Adjustment", 
                          "FR Fiscal Rule")

baseline_gip <- results_df["Baseline", "Average_Debt_Cost"]
results_df$Avg_Add_fiscal_space_baseline <- 1 - (results_df$Average_Debt_Cost / baseline_gip)
results_df["Baseline", "Avg_Add_fiscal_space_baseline"] <- NA

baseline_total_gip <- results_df["Baseline", "Total_Debt_Cost"]
results_df$Contribution_Green_Investment <- (baseline_total_gip - results_df$Total_Debt_Cost) / results_df$Total_Green_Investment
results_df["Baseline", "Contribution_Green_Investment"] <- NA

print(results_df)

results_df_t <- as.data.frame(t(results_df))
rownames(results_df_t) <- c(
  "Total Green Invest.", 
  "Av. Debt Cost", 
  "Tot. Debt Cost", 
  "Tot. Debt Cost / GDP", 
  "Av. Trade Def. /GDP", 
  "Currency Variation", 
  "Av. Interest Fx Debt",
  "Avg Add fiscal space", 
  "Contrib Green Invest"
)

# Latex
# latex_table <- xtable(results_df_t, 
#                       digits = 2,
#                       caption = "Comparison of DR Mechanisms",
#                       label = "tab:results")
# 
# print(latex_table, 
#       include.rownames = TRUE,
#       booktabs = TRUE,
#       sanitize.text.function = identity)
