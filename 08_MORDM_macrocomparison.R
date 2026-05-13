# 08_MORDM_macrocomparison.R
#
# Generates the macrocomparison figure:
#   Top    — lever intensity boxplots by cluster (Wilcoxon significance)
#   Bottom — macro trajectory ribbons (mean ± range) per cluster vs. baseline
#
# Inputs (produced by 04_MORDM_optimize_XLow.R / 05_MORDM_optimize_XHigh.R):
#   Data/optpol_bastransDS2_robustness_{1..5}.RDS
#   Data/optpol_bastransDShigh2_robustness_{1..5}.RDS
#
# Output:
#   Images/macrocomparison.png

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

# ── Model system ───────────────────────────────────────────────────────────────
event1 <- list(triggerDate = 4, reducXrO = "0.025")
SOEM   <- cppMakeSys(fileName   = "model_equations_MORDM.R",
                     reportVars = 3,
                     eventTime  = list(event1))

# ── Baseline parameters (shape = 2, following 01_run_scenarios.R) ─────────────
parms_base <- SOEM$parms
parms_base['lambdatr0']          <- 5
parms_base['lambdatr1']          <- 6
parms_base['lambdatr2']          <- 0.011
parms_base['lambdatr0_adj']      <- 5
parms_base['lambdatr1_adj']      <- 12
parms_base['lambdatr2_adj']      <- 1
parms_base['alpha_tr']           <- 0.00104006
parms_base['beta_tr']            <- 0.3
parms_base['gamma_tr']           <- 2.6
parms_base['delta_tr']           <- 4.99377e-7
parms_base['betadebtSwapFXLgFX'] <- 4
parms_base['betadebtSwapFXBgFX'] <- 4
parms_base['betadebtSwapFXLgFXtr']<- 4
parms_base['shrGrL']             <- 1
parms_base['shrGrLFx']           <- 0.5

# ── 1. Load and consolidate all Borg Pareto sets ───────────────────────────────
rds_files <- c(
  paste0("Data/optpol_bastransDS2_robustness_",     1:5, ".RDS"),
  paste0("Data/optpol_bastransDShigh2_robustness_", 1:5, ".RDS")
)
set_all <- do.call(rbind,
                   lapply(rds_files, function(f) mordm.get.set(readRDS(f))))
set_all <- as.data.frame(set_all)

# Column layout (from 04_MORDM_optimize_XLow.R, colnames set at line 207):
#   1  Foreign_Loans         → shrGrLFx
#   2  Foreign_Bonds         → shrGrBw
#   3  Greenium_Loan         → md_lgtr
#   4  Greenium_Bonds        → md_bgtr
#   5  Interest_Forgiveness  → mdds
#   6  Stock_Forgiveness     → decds
#   7–18 performance indicators (perCapita … Gip)
lever_cols <- 1:6
perf_cols  <- 7:18

# ── 2. Cluster on scaled performance indicators ────────────────────────────────
k.opt <- NbClust(
  as.matrix(scale(set_all[, perf_cols])),
  method = "kmeans", min.nc = 2, max.nc = 6
)$Best.partition

set_all$Cluster <- k.opt
clusters        <- sort(unique(k.opt))
n_clusters      <- length(clusters)

# ── 3. Baseline trajectory (following generateComparisonGraphs.R) ──────────────
baseline  <- cppRK4(SOEM,
                    times     = seq(from = 2019, to = 2030, by = 0.1),
                    eventTime = list(event1),
                    parms     = parms_base) %>%
             as.data.frame()
timeframe <- which(baseline$time >= 2024 & baseline$time <= 2030)

# ── 4. Run all MORDM policies and collect trajectories ─────────────────────────
# Variable list for the bottom panel (12 variables, matching paper figure)
var_run <- c('perCapita', 'inflation', 'reserves', 'foreignDebt',
             'pubDebt',   'fiscalDef', 'hhFrag',   'firmsFrag',
             'unem',      'CAD',       'Gip',       'en')

allResults <- list()
for (nm in var_run)
  allResults[[nm]] <- data.frame(time = baseline$time[timeframe])

for (i in 1:nrow(set_all)) {
  parms_i                      <- parms_base
  parms_i['dsactive']          <- 1
  parms_i['shrGrL']            <- 1
  parms_i['shrGrLFx']          <- set_all[i, "Foreign_Loans"]
  parms_i['shrGrBw']           <- set_all[i, "Foreign_Bonds"]
  parms_i['md_lgtr']           <- set_all[i, "Greenium_Loan"]
  parms_i['md_bgtr']           <- set_all[i, "Greenium_Bonds"]
  parms_i['mdds']              <- set_all[i, "Interest_Forgiveness"]
  parms_i['decds']             <- set_all[i, "Stock_Forgiveness"]

  tmp <- cppRK4(SOEM,
                parms     = parms_i,
                times     = seq(from = 2019, to = 2030, by = 0.1),
                eventTime = list(event1)) %>%
         as.data.frame()

  for (nm in var_run)
    allResults[[nm]][[as.character(i)]] <- tmp[[nm]][timeframe]
}

# ── 5. Build ribbon data per cluster (min / max / mean) ───────────────────────
var_labels <- c(
  perCapita   = "GDP per capita",
  inflation   = "Inflation rate",
  reserves    = "Foreign reserves",
  foreignDebt = "Foreign debt",
  pubDebt     = "Public debt",
  fiscalDef   = "Fiscal deficit",
  hhFrag      = "Household fragility",
  firmsFrag   = "Firm fragility",
  unem        = "Unemployment",
  CAD         = "Current account deficit",
  Gip         = "Government interest payment",
  en          = "Nominal exchange rate"
)

# Colour palette: first two entries match the paper (teal + orange);
# extends automatically for more clusters.
cluster_palette <- c("#2D7B73", "#E8883A", "#5E5EA5", "#C55A4F",
                     "#6DAD6D", "#8B6BAE")
named_palette   <- setNames(cluster_palette[seq_along(clusters)],
                            paste0("Cluster ", clusters))

ribbon_list <- list()
for (cl in clusters) {
  cl_ids <- 1L + which(set_all$Cluster == cl)  # +1: col 1 is time
  for (nm in var_run) {
    mat <- allResults[[nm]][, cl_ids, drop = FALSE]
    ribbon_list[[paste0(nm, "_", cl)]] <- data.frame(
      time     = allResults[[nm]]$time,
      ymin     = apply(mat, 1, min),
      ymax     = apply(mat, 1, max),
      ymean    = apply(mat, 1, mean),
      Cluster  = paste0("Cluster ", cl),
      Variable = var_labels[nm]
    )
  }
}

ribbon_df <- bind_rows(ribbon_list) %>%
  mutate(
    Cluster  = factor(Cluster,  levels = paste0("Cluster ", clusters)),
    Variable = factor(Variable, levels = var_labels)
  )

# Baseline in long format for bottom panel
baseline_bottom <- bind_rows(lapply(var_run, function(nm) {
  data.frame(
    time     = baseline$time[timeframe],
    value    = baseline[[nm]][timeframe],
    Variable = var_labels[nm]
  )
})) %>%
  mutate(Variable = factor(Variable, levels = var_labels))

# ── 6. TOP PANEL: lever intensity boxplots by cluster ─────────────────────────
# Display name for each lever column (matching levers_combined figure in the paper)
lever_display <- c(
  Foreign_Loans        = "Foreign Loans",
  Foreign_Bonds        = "LCY Bonds",
  Greenium_Loan        = "Greenium Loans",
  Greenium_Bonds       = "Greenium Bonds",
  Interest_Forgiveness = "Ir Renegotiation",
  Stock_Forgiveness    = "Principal Adjustment"
)
# Facet order matching levers_combined (row 1: LCY Bonds, Foreign Loans, Greenium Bonds;
#                                        row 2: Greenium Loans, Ir Renegotiation, Principal Adjustment)
lever_order <- c("LCY Bonds",    "Foreign Loans",    "Greenium Bonds",
                 "Greenium Loans", "Ir Renegotiation", "Principal Adjustment")

lever_long <- set_all %>%
  select(all_of(names(lever_display)), Cluster) %>%
  mutate(Cluster = factor(paste0("Cluster ", Cluster),
                          levels = paste0("Cluster ", clusters))) %>%
  pivot_longer(-Cluster, names_to = "Lever", values_to = "Intensity") %>%
  mutate(Lever = factor(lever_display[Lever], levels = lever_order))

# All pairwise comparisons (works for 2 or more clusters)
comparisons <- combn(paste0("Cluster ", clusters), 2, simplify = FALSE)

top_panel <- ggboxplot(
  lever_long,
  x            = "Cluster",
  y            = "Intensity",
  fill         = "Cluster",
  facet.by     = "Lever",
  ncol         = 3,
  outlier.size = 0.3,
  palette      = named_palette
) +
  stat_compare_means(
    comparisons = comparisons,
    method      = "wilcox.test",
    label       = "p.signif"
  ) +
  scale_y_continuous(limits = c(0, 1.25), breaks = seq(0, 1, 0.25)) +
  ylab("Intensity") +
  xlab("") +
  theme(
    legend.position = "none",
    strip.text      = element_text(size = 9),
    axis.text       = element_text(size = 8)
  )

# ── 7. BOTTOM PANEL: macro trajectory ribbons ─────────────────────────────────
bottom_panel <- ggplot() +
  geom_ribbon(
    data  = ribbon_df,
    aes(x = time, ymin = ymin, ymax = ymax, fill = Cluster),
    alpha = 0.35
  ) +
  geom_line(
    data = ribbon_df,
    aes(x = time, y = ymean, color = Cluster),
    linewidth = 0.6
  ) +
  geom_line(
    data      = baseline_bottom,
    aes(x = time, y = value),
    color     = "black",
    linewidth = 0.6
  ) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 3) +
  scale_fill_manual(values  = named_palette) +
  scale_color_manual(values = named_palette) +
  scale_x_continuous(breaks = seq(2024, 2030, 2)) +
  scale_y_continuous(labels = scales::label_number(accuracy = NULL,
                                                    big.mark = ",")) +
  theme_minimal() +
  theme(
    strip.text       = element_text(size = 9),
    axis.text        = element_text(size = 7),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 9)
  ) +
  xlab("") +
  ylab("") +
  labs(fill = "Cluster", color = "Cluster")

# ── 8. Combine and save ────────────────────────────────────────────────────────
combined <- plot_grid(top_panel, bottom_panel,
                      nrow = 2, rel_heights = c(2, 5))

ggsave("Images/macrocomparison.png", combined,
       width = 18, height = 20, dpi = 300)
