# Chronicle of a Debt Foretold

Replication code for the paper **"Chronicle of a Debt Foretold"**.

This repository contains the model, analysis pipeline, and robustness analysis for GEMMES Colombia,
Stock-Flow Consistent (SFC) macro-financial model of the Colombian economy, used to assess green 
debt levers and their macroeconomic implications.

---

## Model

The model is a multi-sector SFC model calibrated on Colombian national accounts and
financial data. It tracks real-financial interactions across households, non-financial
corporations, the government, commercial banks, and the rest of the world, with an
explicit green investment block representing NDC-aligned transition expenditures.

Decision variables (debt restructuring levers):
- **shrGrLFx** — share of new green financing via foreign loans
- **shrGrBw** — share via local-currency bonds
- **md_lgtr / md_bgtr** — greenium on loans / bonds
- **mdds** — interest rate reduction on existing debt (debt swap)
- **decds** — principal reduction on existing debt (debt swap)

Performance objectives (12): GDP per capita growth, inflation, foreign reserves, foreign
debt, private debt, public debt, fiscal deficit, household fragility, firm fragility,
unemployment, current account, government interest payments.

---

## Repository structure

```
.
├── model_equations_DebtSwap.R        # SFC model equations — baseline + DR scenarios
├── model_equations_MORDM.R           # SFC model equations — variant for MORDM optimization
│
├── 01_run_scenarios.R                # Run baseline + DR1–DR4 + fiscal rule; produce figures
├── 02_sensitivity_baseline.R         # LHS + OT sensitivity analysis on baseline parameters
├── 03_sensitivity_scenarios.R        # LHS + OT sensitivity analysis across DR scenarios
│
├── 04_MORDM_optimize_XLow.R          # Borg MOEA optimization (parallel) — low baseline set
├── 05_MORDM_optimize_XHigh.R         # Borg MOEA optimization — high baseline set
├── 06_MORDM_robustness_calibration.R # LHS sampling for deep-uncertainty robustness evaluation
├── 07_MORDM_robustness_metrics.R     # Regret + satisficing robustness metrics; cluster analysis
│
├── Extrafunctions.R                  # ggplot parallel coordinates and visualization helpers
│
├── Source/
│   ├── SourceCode.R                  # Core engine: cppMakeSys(), cppRK4(), cppCompileRK4()
│   ├── sourceCodeCalibration.R       # Calibrated parameter vector (parms_NewC)
│   ├── utilities.R                   # mordm.get.set() and MORDM support functions
│   ├── CppCodeRK4.cpp                # C++ RK4 solver (generated; do not edit)
│   ├── RawCppCodeRK4.cpp             # C++ RK4 solver template
│   ├── CppCodeMinDist.cpp            # C++ minimum-distance solver (generated)
│   └── RawCppCodeMinDist.cpp         # C++ minimum-distance solver template
│
├── Data/
│   ├── XLow_ord.csv                  # 5 alternative baseline parameter sets (low)
│   ├── XHigh_ord.csv                 # 5 alternative baseline parameter sets (high)
│   ├── parms.csv                     # Baseline parameter values
│   ├── parms1.csv                    # Calibration parameter set
│   └── res_indices.csv               # Variable index reference
│
├── Robustness/                       # Generated outputs (RDS) — see pipeline below
└── Images/                           # Generated figures (PNG)
```

---

## Pipeline

Scripts are numbered in execution order. Steps 4–5 are computationally intensive and
can be run in parallel (they operate on independent baseline sets).

```
01_run_scenarios.R
    └─ sources: model_equations_DebtSwap.R
    └─ outputs: Images/clean_plot_res*.png

02_sensitivity_baseline.R
    └─ sources: 01_run_scenarios.R
    └─ outputs: console / plots

03_sensitivity_scenarios.R
    └─ sources: 01_run_scenarios.R
    └─ outputs: console / plots

04_MORDM_optimize_XLow.R              ┐  independent,
05_MORDM_optimize_XHigh.R             ┘  can run in parallel
    └─ sources: model_equations_MORDM.R
    └─ inputs:  Data/XLow_ord.csv / Data/XHigh_ord.csv
    └─ outputs: Data/optpol_bastransDS2_robustness_{1..5}.RDS
                Data/optpol_bastransDShigh2_robustness_{1..5}.RDS

06_MORDM_robustness_calibration.R
    └─ inputs:  Data/optpol_bastransnew2*.RDS  (consolidated from step 4/5)
    └─ outputs: Robustness/RobcalsShape_no_new2*.RDS

07_MORDM_robustness_metrics.R
    └─ inputs:  Robustness/roblist_no_new2*.RDS
                Robustness/RobcalsShape_no_new2*.RDS
                Data/optpol_bastransnew2*.RDS
    └─ outputs: Robustness/reg1_threshresnewGD2.RDS
                Robustness/reg2_threshresnewGD2.RDS
                Robustness/sat1_threshresnewGD2.RDS
                Robustness/sat2_threshresnewGD2.RDS
```

---

## Dependencies

### R packages

```r
install.packages(c(
  "Rcpp", "tidyverse", "data.table", "ggplot2", "gridExtra",
  "rlist", "xtable", "readxl", "readr", "rmarkdown",
  "lhs", "sensobol", "gsaot",
  "parallel", "foreach", "future", "future.apply",
  "dtwclust", "factoextra", "FactoMineR", "kernlab",
  "corrplot", "viridis", "htmlwidgets", "webshot2", "ggforce"
))
```

### OpenMORDM and Borg MOEA

Scripts 04–07 require **OpenMORDM** and the **Borg MOEA** C library. OpenMORDM is
an R package wrapping the Borg multi-objective evolutionary algorithm and providing
MORDM utilities. It is not on CRAN; installation instructions are available at:

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

