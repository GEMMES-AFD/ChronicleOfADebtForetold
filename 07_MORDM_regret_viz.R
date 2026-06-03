# 09_MORDM_regret_viz.R
#
# Generates the Regret II scatter plot (Figure: robustness analysis).
#   x-axis: Average Greenium (%) = mean(Greenium_Loan, Greenium_Bonds)
#   y-axis: Average Reprofiling (%) = mean(Interest_Forgiveness, Stock_Forgiveness)
#   size:   Regret category (larger = lower regret = more robust)
#   color:  FX Loans vs. Peso Bonds (%) = Foreign_Loans − Foreign_Bonds
#
# Input:  Robustness/reg2_threshresnewGD2.RDS  (from 07_MORDM_robustness_metrics.R)
# Output: Images/regret2.png

# ── Libraries ──────────────────────────────────────────────────────────────────
library(tidyverse)
REG2_FILE <- "MORDM_Process/MORDM_Results/Regret2_REAC.RDS" 
# ── 1. Load regret data ────────────────────────────────────────────────────────
# reg2 = cbind(set, Regret): columns 1–6 are levers, 7–18 are performance,
# last column is Regret (90th-percentile relative degradation, already in % units)
reg2 <- readRDS(REG2_FILE) %>% as.data.frame()
colnames(reg2)[2] <- "Foreign_Bonds"
colnames(reg2)[5] <- "Interest_Forgiveness"
colnames(reg2)[6] <- "Stock_Forgiveness"
# ── 2. Categorise regret into size bins ────────────────────────────────────────
max_regret <- ceiling(max(reg2$Regret, na.rm = TRUE))

reg2_df <- reg2 %>%
  mutate(
    Regret_cat = cut(
      Regret,
      breaks        = c(0, 2, 5, 8, max_regret),
      include.lowest = TRUE,
      right         = FALSE,
      labels        = c("0–2", "2–5", "5–8",
                        paste0("8–", max_regret))
    )
  ) %>%
  arrange(desc(Regret_cat))   # plot low-regret (large) points last so they appear on top

# ── 3. Build plot ──────────────────────────────────────────────────────────────
# Column names follow 04_MORDM_optimize_XLow.R conventions (underscores, no spaces):
#   Foreign_Loans, Foreign_Bonds, Greenium_Loan, Greenium_Bonds,
#   Interest_Forgiveness, Stock_Forgiveness

p_regret2 <- ggplot(
  data = reg2_df,
  aes(
    x     = 100 * (Greenium_Loan   + Greenium_Bonds)       / 2,
    y     = 100 * (Stock_Forgiveness + Interest_Forgiveness) / 2,
    color = 100 * (Foreign_Loans   - Foreign_Bonds)         / 2,
    size  = Regret_cat
  )
) +
  geom_point() +
  scale_size_manual(
    name   = "Regret",
    values = setNames(
      c(6, 3, 2, 1.5),
      levels(reg2_df$Regret_cat)
    )
  ) +
  scale_color_gradient(
    name  = "FX Loans vs. Peso Bonds(%)",
    low   = "#F8A049FF",   # orange = FX-loan-heavy
    high  = "#2E5F8DFF",   # blue   = peso-bond-heavy
    guide = guide_colorbar(ticks = FALSE, label = FALSE)
  ) +
  scale_fill_discrete(guide = "none") +
  xlab("Average Greenium (%)") +
  ylab("Average Reprofiling (%)") +
  ggtitle("Regret II") +
  theme(
    plot.title       = element_text(size = 10),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90"),
    axis.title       = element_text(size = 10),
    legend.title     = element_text(size = 10, margin = margin(b = 10)),
    axis.text        = element_text(size = 10),
    legend.text      = element_text(size = 10)
  )
p_regret2

# ── 4. Save ────────────────────────────────────────────────────────────────────
ggsave("Figures/regret2.png", p_regret2, dpi = 400, width = 7.5, height = 5)
