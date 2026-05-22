# 02_sensitivity_baseline.R
#
# Three-step sensitivity and baseline replication analysis.
#
# Step 1 — Identify stiff and sloppy parameters
#   K independent LHS samples (n = 5 000 each, ±25% around baseline).
#   Scalar output: full time-series Euclidean distance to baseline
#                  (summed over all time steps × 11 variables).
#   OT index (gsaot::ot_indices_1d) per sample → rank stability across K.
#
# Step 2 — Alternative calibrations via focused LHS (n = 20 000, ±25%)
#   Stiff params varied, sloppy fixed  → sort by distance → XHigh_ord.csv
#   Sloppy params varied, stiff fixed  → sort by distance → XLow_ord.csv
#   These files feed scripts 04 / 05 (Borg optimisation).
#
# Step 3 — Baseline replication figures (appendix)
#   Top 5 alternatives from each set plotted against baseline.
#
# Inputs:  sourced from 01_run_scenarios.R
# Outputs: Data/XHigh_ord.csv, Data/XLow_ord.csv, Figures/*.png

library(lhs)
library(dplyr)
library(gsaot)

source("01_run_scenarios.R")

# ── Shared setup ──────────────────────────────────────────────────────────────
output_names <- c(
  "perCapita", "inflation", "reserves", "foreignDebt", "pubDebt",
  "fiscalDef", "hhFrag", "firmsFrag", "unem", "CAD", "GipGDP"
)

param_names <- names(parms1)
baseline    <- as.numeric(parms1)
names(baseline) <- param_names
cat("Parameter vector length:", length(baseline), "\n")

# LHS bounds: ±25% symmetric
lower_raw <- baseline * 0.75
upper_raw <- baseline * 1.25
lower <- pmin(lower_raw, upper_raw)
upper <- pmax(lower_raw, upper_raw)
eps0  <- 1e-6
lower[baseline == 0] <- -eps0
upper[baseline == 0] <-  eps0

# Parameters held constant in all LHS (switches, levers, scenario-specific)
constantParms <- c(
  "dsactive", "FRactive", "reaction",
  "alpha_tr", "beta_tr", "gamma_tr", "delta_tr",
  "lambdatr0", "lambdatr1", "lambdatr2",
  "lambdatr0_adj", "lambdatr1_adj", "lambdatr2_adj",
  "lambdaicf",
  "shrGrTax", "shrGrIC", "shrDon",
  "shrGrBw", "shrGrL", "shrGrLFx",
  "md_lgtr", "md_bgtr",
  "betadebtSwapFXLgFX", "betadebtSwapFXLgFXtr", "betadebtSwapFXBgFX",
  "mdds", "decds",
  "taumtr", "tauothktr", "betasigmamktr",
  "zetafx2", "K_0",
  "rho", "zetafx3", "gammariskFFX", "gammariskBFX",
  "triggerFX", "speedFX", "shareFX",
  "speedTaum", "triggerTaum", "maxTaum",
  "speedWg",   "triggerWG",  "maxWg",
  "speedCg",   "triggerCg",  "maxCg",
  "triggerSTg","speedSTg",   "maxSTg",
  "tauCBAM0", "tauCBAM1",
  "sigmaxnpNew", "sigmaxnSpeed", "sigmaxnInit",
  "sigmapcNew", "sigmapicNew", "sigmapkNew",
  "sigmamSpeed", "sigmamInit",
  "sigmamktr0", "sigmaaktr", "epsilon2ktr",
  "atr0", "atr1", "scenInv"
)
constantParms <- constantParms[constantParms %in% param_names]
cat("Constant:", length(constantParms),
    "| Varied:", length(param_names) - length(constantParms), "\n")

# Baseline trajectory (reference for all distance calculations)
bas_df <- cppRK4(SOEM1, parms = parms1, eventTime = list(event1)) %>%
  as.data.frame() %>%
  filter(time > 2022.9 & time < 2030) %>%
  mutate(CAD = (X - IM) / GDP, GipGDP = Gip / GDP)

# Output statistics names: per-variable mean + dist_baseline
stat_names <- c(paste0(output_names, "_mean"), "dist_baseline")

# Model runner: returns per-variable means + full time-series distance
run_stats <- function(par_vec_named) {
  parms_i                       <- baseline
  parms_i[names(par_vec_named)] <- par_vec_named
  out <- setNames(rep(NA_real_, length(stat_names)), stat_names)
  sim <- tryCatch(
    cppRK4(SOEM1, parms = parms_i, eventTime = list(event1)) %>%
      as.data.frame() %>%
      filter(time > 2022.9 & time < 2030) %>%
      mutate(CAD = (X - IM) / GDP, GipGDP = Gip / GDP),
    error = function(e) NULL
  )
  if (is.null(sim) || nrow(sim) != nrow(bas_df)) return(out)
  for (v in output_names)
    out[paste0(v, "_mean")] <- mean(sim[[v]], na.rm = TRUE)
  eps <- 1e-5; dist_sq <- 0
  for (v in output_names)
    dist_sq <- dist_sq +
      sum(((sim[[v]] - bas_df[[v]]) / (abs(bas_df[[v]]) + eps))^2, na.rm = TRUE)
  out["dist_baseline"] <- sqrt(dist_sq)
  out
}

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Multi-criteria OT sensitivity (K independent LHS samples)
#   Criteria: dist_baseline + univariate OT per output + multivariate OT
#   Selection: parameters consistently influential across ≥ 2 criteria
# ═══════════════════════════════════════════════════════════════════════════════

K      <- 5
n_runs <- 5000

# criteria: dist_baseline + all outputs (univariate) + "multi"
criteria    <- c("dist_baseline", output_names, "multi")
res_all     <- vector("list", K)   # res_all[[k]][[criterion]] = named OT index vector
valid_pcts  <- numeric(K)

for (k in seq_len(K)) {
  cat(sprintf("\n── LHS sample %d / %d ──\n", k, K))
  set.seed(k * 42)

  U <- randomLHS(n_runs, length(baseline))
  X <- sweep(U, 2, (upper - lower), `*`)
  X <- sweep(X, 2, lower,           `+`)
  colnames(X) <- param_names
  for (cp in constantParms) X[, cp] <- baseline[cp]

  Y <- matrix(NA_real_, nrow = n_runs, ncol = length(stat_names))
  colnames(Y) <- stat_names
  for (i in seq_len(n_runs)) {
    par_i        <- X[i, ]
    names(par_i) <- param_names
    Y[i, ]       <- run_stats(par_i)
  }

  valid <- is.finite(Y[, "dist_baseline"])
  valid_pcts[k] <- 100 * mean(valid)
  cat(sprintf("  Valid: %d / %d (%.1f%%)\n", sum(valid), n_runs, valid_pcts[k]))

  x_df <- as.data.frame(X[valid, , drop = FALSE])
  colnames(x_df) <- make.unique(colnames(x_df))
  x_df <- x_df[, apply(x_df, 2, var) > 0, drop = FALSE]
  Y_clean <- Y[valid, , drop = FALSE]

  res_k <- list()

  # OT on dist_baseline
  res_dist         <- ot_indices_1d(x = x_df, y = Y_clean[, "dist_baseline"], M = 20)
  res_k[["dist_baseline"]] <- res_dist$indices
  plot(res_dist, ranking = 20, main = paste0("OT dist_baseline — s", k))

  # OT univariate per output
  for (v in output_names) {
    y_v       <- Y_clean[, paste0(v, "_mean")]
    if (var(y_v, na.rm = TRUE) == 0) next
    res_v     <- ot_indices_1d(x = x_df, y = y_v, M = 20)
    res_k[[v]] <- res_v$indices
  }

  # OT multivariate (all outputs jointly)
  y_mat       <- Y_clean[, paste0(output_names, "_mean"), drop = FALSE]
  res_multi   <- ot_indices(x = x_df, y = y_mat, M = 20)
  res_k[["multi"]] <- res_multi$indices
  plot(res_multi, ranking = 20, main = paste0("OT multivariate — s", k))

  res_all[[k]] <- res_k
}

# ── Rank stability and multi-criteria aggregation ─────────────────────────────
# For each criterion, compute rank of each parameter across K samples.
# A parameter is "stiff" if it ranks in the top 30 across ≥ 2 criteria.

common_params <- Reduce(intersect,
  lapply(res_all, function(rk) Reduce(intersect, lapply(rk, names))))

# For each criterion: mean rank and sd across K samples
crit_stability <- lapply(criteria, function(cr) {
  rmat <- do.call(cbind, lapply(res_all, function(rk) {
    idx <- rk[[cr]]
    if (is.null(idx)) return(rep(NA_real_, length(common_params)))
    rank(-idx[common_params], ties.method = "average")
  }))
  rownames(rmat) <- common_params
  data.frame(
    param     = common_params,
    criterion = cr,
    rank_mean = rowMeans(rmat, na.rm = TRUE),
    rank_sd   = apply(rmat, 1, sd, na.rm = TRUE)
  )
})
crit_df <- bind_rows(crit_stability)

# Count in how many criteria each parameter ranks in the top 30
top30_counts <- crit_df %>%
  group_by(param) %>%
  summarise(
    n_top30    = sum(rank_mean <= 30, na.rm = TRUE),
    rank_dist  = rank_mean[criterion == "dist_baseline"],
    rank_multi = rank_mean[criterion == "multi"],
    rank_sd_dist = rank_sd[criterion == "dist_baseline"],
    .groups = "drop"
  ) %>%
  arrange(desc(n_top30), rank_dist)

cat("\n══ Multi-criteria summary (top 50) ══\n")
cat("n_top30 = number of criteria (out of", length(criteria),
    ") for which mean rank ≤ 30\n\n")
print(as.data.frame(head(top30_counts, 50)), row.names = FALSE)

# ── SET AFTER REVIEWING ───────────────────────────────────────────────────────
# Stiff candidates: n_top30 ≥ 2  (consistent across multiple criteria)
# Sloppy candidates: n_top30 = 0 AND rank_dist high
# Adjust manually if domain knowledge justifies it.

top_30_stiff  <- top30_counts$param[top30_counts$n_top30 >= 2][1:30]
top_30_sloppy <- top30_counts %>%
  filter(n_top30 == 0) %>%
  arrange(desc(rank_dist)) %>%
  pull(param) %>%
  head(30)

cat("\n── Stiff (n_top30 ≥ 2, top 30) ──\n");  dput(top_30_stiff)
cat("\n── Sloppy (n_top30 = 0, bottom 30 on dist) ──\n"); dput(top_30_sloppy)

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2 — Focused LHS (n = 20 000) on stiff and sloppy separately
#           → alternative calibrations sorted by distance to baseline
# ═══════════════════════════════════════════════════════════════════════════════

run_focused_lhs <- function(params_sel, n = 20000, seed = 1) {
  set.seed(seed)
  baseline_sel <- baseline[params_sel]
  lo <- pmin(baseline_sel * 0.75, baseline_sel * 1.25)
  hi <- pmax(baseline_sel * 0.75, baseline_sel * 1.25)
  lo[baseline_sel == 0] <- -eps0
  hi[baseline_sel == 0] <-  eps0

  U <- randomLHS(n, length(params_sel))
  X <- sweep(U, 2, (hi - lo), `*`)
  X <- sweep(X, 2, lo,        `+`)
  colnames(X) <- params_sel

  y <- numeric(n)
  for (i in seq_len(n)) {
    par_i        <- X[i, ]
    names(par_i) <- params_sel
    y[i]         <- run_stats(par_i)["dist_baseline"]
  }

  list(X = X, y = y)
}

# Stiff params → XHigh
cat("\n── Step 2a: LHS on stiff parameters (n = 20 000) ──\n")
res_high <- run_focused_lhs(top_30_stiff, n = 20000, seed = 101)
valid_high <- is.finite(res_high$y)
cat(sprintf("Valid: %d / 20000 (%.1f%%)\n", sum(valid_high), 100*mean(valid_high)))
X_high_ord <- res_high$X[order(res_high$y), ]   # sorted ascending (closest first)
write.table(X_high_ord, file = "Data/XHigh_ord.csv", sep = ";", dec = ",")
cat("Saved → Data/XHigh_ord.csv\n")

# Sloppy params → XLow
cat("\n── Step 2b: LHS on sloppy parameters (n = 20 000) ──\n")
res_low <- run_focused_lhs(top_30_sloppy, n = 20000, seed = 202)
valid_low <- is.finite(res_low$y)
cat(sprintf("Valid: %d / 20000 (%.1f%%)\n", sum(valid_low), 100*mean(valid_low)))
X_low_ord <- res_low$X[order(res_low$y), ]
write.table(X_low_ord, file = "Data/XLow_ord.csv", sep = ";", dec = ",")
cat("Saved → Data/XLow_ord.csv\n")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3 — Baseline replication figures (top 5 alternatives vs baseline)
# ═══════════════════════════════════════════════════════════════════════════════

make_replication_df <- function(X_ord, label_prefix) {
  res_list    <- list()
  res_list[[1]] <- cppRK4(SOEM1, parms = parms1, eventTime = list(event1)) %>%
    as.data.frame() %>%
    mutate(Scenario = "Baseline") %>%
    filter(time > 2022.9 & time < 2030)

  for (i in 1:5) {
    parms_i                       <- parms1
    parms_i[names(X_ord[i, ])]   <- X_ord[i, ]
    sim <- cppRK4(SOEM1, parms = parms_i, eventTime = list(event1)) %>%
      as.data.frame() %>%
      mutate(Scenario = paste0(label_prefix, i)) %>%
      filter(time > 2022.9 & time < 2030)
    res_list[[i + 1]] <- sim
  }

  res_list %>%
    list.rbind() %>%
    mutate(
      Scenario = as.factor(Scenario),
      GipGDP   = 100 * Gip / GDP,
      en       = en * 4300 / 1.271753,
      rsk      = 100 * rsk,
      FIP      = 100 * FIP,
      premgd   = 100 * premgd
    ) %>%
    dplyr::select(Scenario, time, all_of(diagvars)) %>%
    pivot_longer(3:(length(diagvars) + 2),
                 names_to = "Variable", values_to = "Value") %>%
    mutate(Variable = factor(Variable, levels = diagvars)) %>%
    mutate(Scenario = factor(Scenario,
                             levels = c("Baseline",
                                        paste0(label_prefix, 1:5))))
}

plot_replication <- function(res_df, title) {
  ggplot(res_df,
         aes(x = time, y = Value, group = Scenario,
             color = Scenario, linetype = Scenario)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ Variable, scales = "free_y",
               labeller = as_labeller(diagvar_labs), ncol = 3) +
    scale_y_continuous(labels = scales::label_number(accuracy = NULL,
                                                     big.mark = ",")) +
    theme_classic() +
    theme(
      strip.text        = element_text(colour = "black", size = 10),
      strip.background  = element_blank(),
      axis.text         = element_text(size = 8),
      legend.text       = element_text(size = 9),
      legend.title      = element_blank(),
      legend.box        = "vertical",
      legend.spacing.y  = unit(0.5, "cm"),
      legend.key.height = unit(0.8, "cm")
    ) +
    ylab("") + xlab("Year") + ggtitle(title)
}

cat("\n── Step 3a: Replication figure — stiff alternatives ──\n")
df_high <- make_replication_df(X_high_ord, "Stiff_")
print(plot_replication(df_high, "Baseline replication — stiff parameter alternatives"))
ggsave("Figures/clean_plot_alternativeHigh.png", width = 10, height = 5, dpi = 400)

cat("── Step 3b: Replication figure — sloppy alternatives ──\n")
df_low <- make_replication_df(X_low_ord, "Sloppy_")
print(plot_replication(df_low, "Baseline replication — sloppy parameter alternatives"))
ggsave("Figures/clean_plot_alternativeLow.png", width = 10, height = 5, dpi = 400)

cat("\nDone. XHigh_ord.csv and XLow_ord.csv are ready for scripts 04 and 05.\n")
