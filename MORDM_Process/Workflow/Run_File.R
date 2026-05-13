cat("Step 0: Calling Globals...\n")
source("MORDM_Process/Workflow/00-Globals.R")

cat("Step 1: Creating SOS...\n")
source("MORDM_Process/Workflow/01-SOS.R")

cat("Step 2: Getting optimal policies...\n")
source("MORDM_Process/Workflow/02-MORDM_optpol.R")

cat("Step 3: Generating robustness calibration...\n")
source("MORDM_Process/Workflow/03-MORDM_Robcal.R")

cat("Step 4: Running robustness Runs...\n")
source("MORDM_Process/Workflow/04-MORDM_RobRes.R")

cat("Step 5: Computing Regret Metrics\n")
source("MORDM_Process/Workflow/05-MORDM_RobMetrics.R")