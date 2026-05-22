# Chronicle of a Debt Foretold

Replication code for the paper **"Chronicle of a Debt Foretold"**.

This repository contains the model, analysis pipeline, and robustness analysis for GEMMES Colombia,
a Stock-Flow Consistent (SFC) macro-financial model of the Colombian economy, used to assess green
debt restructuring levers and their macroeconomic implications. Attribution 4.0 International CC BY 4.0.

---

## Model

The model is a multi-sector SFC model calibrated on Colombian national accounts and
financial data. It tracks real-financial interactions across households, non-financial
corporations, the government, commercial banks, and the rest of the world, with an
explicit green investment block representing NDC-aligned transition expenditures.

Decision variables (debt restructuring levers):

| Parameter | Lever |
|---|---|
| `shrGrLFx` | Share of new green financing via foreign loans |
| `shrGrBw` | Share via local-currency bonds |
| `md_lgtr` / `md_bgtr` | Greenium on loans / bonds |
| `mdds` | Interest rate reduction on existing debt (debt swap) |
| `decds` | Principal reduction on existing debt (debt swap) |

Performance objectives (12): GDP per capita, inflation, foreign reserves, foreign debt,
public debt, fiscal deficit, household fragility, firm fragility, unemployment, current
account deficit, government interest payments, nominal exchange rate.

---

## Repository structure

```
.
├── model_equations_MORDM.R               # Unified SFC model equations — baseline, DR scenarios,
│                                         # Fiscal Rule and MORDM optimization (switches inside)
│
├── 01_run_scenarios.R                    # Run baseline + DR1–DR4 + fiscal rule; produce Figure 1
├── 02_sensitivity_baseline.R             # LHS + OT sensitivity (K=5 runs, ±25%)
│                                         #   Step 1 — stiff/sloppy parameter identification
│                                         #   Step 2 — focused LHS → Data/XHigh_ord.csv, XLow_ord.csv
│                                         #   Step 3 — appendix replication figures
│
├── 04_MORDM_optimize_XLow.R              # Borg MOEA on sloppy alternative calibrations (XLow)
│                                         #   → sensitivity of Pareto front to sloppy parameters
├── 05_MORDM_optimize_XHigh.R             # Borg MOEA on stiff alternative calibrations (XHigh)
│                                         #   → sensitivity of Pareto front to stiff parameters
├── 06_MORDM_robustness_calibration.R     # LHS sampling for deep-uncertainty robustness evaluation
├── 07_MORDM_robustness_metrics.R         # Regret + satisficing robustness metrics; cluster analysis
│
├── 08_MORDM_PLOT_macrocomparison.R       # Macrocomparison figure (lever boxplots + ribbon trajectories)
├── 09_MORDM_regret_viz.R                 # Regret II scatter plot
│
├── Extrafunctions.R                      # ggplot parallel coordinates and visualization helpers
│
├── Source/
│   ├── SourceCode.R                      # Core engine: cppMakeSys(), cppRK4(), cppCompileRK4()
│   ├── utilities.R                       # mordm.get.set() and MORDM support functions
│   ├── CppCodeRK4.cpp                    # C++ RK4 solver (generated at runtime; committed as reference)
│   ├── RawCppCodeRK4.cpp                 # C++ RK4 solver template
│   └── RawCppCodeMinDist.cpp             # C++ minimum-distance solver template
│
├── Data/
│   ├── XLow_ord.csv                      # Alternative calibrations (sloppy params varied), sorted by
│   │                                     # distance to baseline; input for script 04
│   ├── XHigh_ord.csv                     # Alternative calibrations (stiff params varied); input for 05
│   ├── parms.csv                         # Baseline parameter values
│   ├── parms1.csv                        # Calibration parameter set
│   └── res_indices.csv                   # Variable index reference
│
├── Figures/                              # Generated figures (committed)
│   ├── Figure1.png                       # Scenario comparison (output of 01)
│   └── clean_plot_res18.png              # Baseline time series (output of 01)
│
├── MORDM_Process/MORDM_Results/          # Authoritative pre-computed Borg optimization results
│   ├── optpol_DS.RDS                     #   Main Pareto-optimal policy set (debt swap strategies)
│   ├── optpol_REAC.RDS                   #   Reactive policy set (produced by MORDM_Process workflow)
│   └── Regret2_REAC.RDS                 #   Regret II robustness metrics for the reactive policies
│                                         # Scripts 08 and 09 read from this directory.
│                                         # The workflow scripts in MORDM_Process/Workflow/ that
│                                         # produced these results are archived there for reference.
│
└── Robustness/                           # Generated RDS outputs from scripts 06–07 (gitignored)
```

> **Two distinct MORDM tracks.** Scripts 04–07 test the *sensitivity of the Pareto front* to
> parametric uncertainty (sloppy/stiff parameter subspaces from script 02). They are not the
> source of the main paper figures. The authoritative Pareto-optimal strategies and robustness
> results used in all main figures are pre-computed in `MORDM_Process/MORDM_Results/` and were
> produced by the workflow archived in `MORDM_Process/Workflow/`. Scripts 08 and 09 read
> exclusively from `MORDM_Process/MORDM_Results/`.

---

## Pipeline

### Track A — Model scenarios and parametric sensitivity (fully replicable)

```
01_run_scenarios.R
    └─ sources: model_equations_MORDM.R
    └─ outputs: Figures/Figure1.png
                Figures/clean_plot_res18.png

02_sensitivity_baseline.R
    └─ sources: 01_run_scenarios.R
    └─ step 1 outputs: console + OT index plots (K=5 samples)
    └─ step 2 outputs: Data/XHigh_ord.csv   (stiff alternative calibrations)
                       Data/XLow_ord.csv    (sloppy alternative calibrations)
    └─ step 3 outputs: Figures/clean_plot_alternativeHigh.png
                       Figures/clean_plot_alternativeLow.png
```

### Track B — Pareto-front sensitivity to parametric uncertainty (requires Borg license)

These scripts test whether the Pareto front structure is robust across the sloppy and stiff
parameter subspaces identified in script 02. They are not the source of the main paper figures.

```
04_MORDM_optimize_XLow.R              ┐  independent,
05_MORDM_optimize_XHigh.R             ┘  can run in parallel
    └─ sources: model_equations_MORDM.R
    └─ inputs:  Data/XLow_ord.csv (script 04) / Data/XHigh_ord.csv (script 05)
    └─ outputs: Data/optpol_bastransDSLow_{shape}_robustness_{kk}.RDS
                Data/optpol_bastransDSHigh_{shape}_robustness_{kk}.RDS
                (gitignored; large intermediate files)

06_MORDM_robustness_calibration.R
    └─ inputs:  Data/ consolidated RDS from scripts 04/05
    └─ outputs: Robustness/RobcalsShape_no_new2{shape}.RDS

07_MORDM_robustness_metrics.R
    └─ inputs:  Robustness/RobcalsShape_no_new2{shape}.RDS + Data/ consolidated RDS
    └─ outputs: Robustness/reg1_threshresnewGD2.RDS  ...  sat2_threshresnewGD2.RDS
```

### Track C — Main MORDM figures (pre-computed; requires Borg license to re-run)

The authoritative optimization results are archived in `MORDM_Process/MORDM_Results/`.
The workflow that produced them is archived in `MORDM_Process/Workflow/`. Scripts 08 and 09
read from `MORDM_Process/MORDM_Results/` and do not depend on scripts 04–07.

```
08_MORDM_PLOT_macrocomparison.R
    └─ inputs:  MORDM_Process/MORDM_Results/optpol_DS.RDS
                MORDM_Process/MORDM_Results/optpol_REAC.RDS
                MORDM_Process/MORDM_Results/Regret2_REAC.RDS
    └─ outputs: Figures/macrocomparison.png
                Figures/Reg2_1v4{shape}.png

09_MORDM_regret_viz.R
    └─ inputs:  MORDM_Process/MORDM_Results/Regret2_REAC.RDS
    └─ outputs: Figures/regret2.png
```

---

## Dependencies

### R packages

```r
install.packages(c(
  "Rcpp", "tidyverse", "data.table", "ggplot2", "gridExtra",
  "rlist", "xtable", "readxl", "readr",
  "lhs", "gsaot",
  "parallel", "foreach", "future", "future.apply",
  "NbClust", "factoextra", "dtwclust",
  "ggpubr", "cowplot", "viridis", "ggforce",
  "htmlwidgets", "webshot2"
))
```

### OpenMORDM and Borg MOEA

> **Replication note:** Scripts 04–09 require the **Borg MOEA** library, which is distributed
> under a license that prohibits redistribution. The compiled Borg library is therefore not
> included in this repository. A free license can be obtained directly from the authors at
> http://borgmoea.org. Without it, Tracks B and C cannot be re-run from scratch; however,
> the pre-computed results in `MORDM_Process/MORDM_Results/` allow scripts 08 and 09 (the
> figure-generation scripts) to run without a Borg license.

Scripts 04–07 additionally require **OpenMORDM**, an R package wrapping the Borg MOEA and
providing MORDM utilities. It is not on CRAN; installation instructions are available at:

- Borg MOEA: http://borgmoea.org
- OpenMORDM: https://github.com/dhadka/OpenMORDM

`rdyncall` is also required and available on CRAN.

### C++ compilation

The SFC model is compiled at runtime by `Source/SourceCode.R` via `cppMakeSys()`,
which generates and compiles C++ code from the R model equation files. No manual
compilation step is required beyond having a working C++ toolchain (Rtools on Windows).

---

## Scenarios

| Label | Description |
|---|---|
| Baseline | No debt restructuring; NDC investment financed at market rates |
| DR1 — Greenium on Loans | New green loans at preferential rate (greenium) |
| DR2 — LCY Bonds | Green financing via local-currency bonds (original sin reduction) |
| DR3 — Interest Rate Renegotiation | Reduction of interest rate on existing foreign debt |
| DR4 — Principal Adjustment | Partial principal cancellation on existing foreign debt |
| FR — Fiscal Rule | Fiscal consolidation rule |

---

## Authors

Morgane Gonon, Antoine Godin, Louis Daumas
