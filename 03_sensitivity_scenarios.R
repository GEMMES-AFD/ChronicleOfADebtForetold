library(lhs)
library(dplyr)
library(gsaot)

source("01_run_scenarios.R")

n_runs <- 5000

output_names <- c(
  "perCapita", "inflation", "reserves", "foreignDebt", "privateDebt", "pubDebt",
  "fiscalDef", "hhFrag", "firmsFrag", "unem", "CAD", "GipGDP"
)

output_stats <- as.vector(t(outer(output_names, c("min", "mean", "max"), paste, sep = "_")))


scenario_list <- list(
  baseline = list(parms = parms1, event = event1),
  DR1 = list(parms = parms4, event = event1),
  DR2 = list(parms = parms5, event = event1),
  DR3 = list(parms = parms2, event = event1),
  DR4 = list(parms = parms3, event = event1)
)

all_results <- list()

for (sc_name in names(scenario_list)) {
  
  sc        <- scenario_list[[sc_name]]
  baseline  <- as.numeric(sc$parms)
  param_names <- names(sc$parms)
  sc_event  <- sc$event
  
  # LHS sampling (same logic as before)
  U <- randomLHS(n_runs, length(baseline))
  lower_raw <- baseline * 0.9
  upper_raw <- baseline * 1.1
  lower <- pmin(lower_raw, upper_raw)
  upper <- pmax(lower_raw, upper_raw)
  eps0 <- 1e-6
  idx0 <- baseline == 0
  lower[idx0] <- -eps0
  upper[idx0] <-  eps0
  X <- sweep(U, 2, (upper - lower), `*`)
  X <- sweep(X, 2, lower, `+`)
  colnames(X) <- param_names
  
  # Model runner using the scenario-specific event
  run_model <- function(par_vec_named) {
    sim <- cppRK4(
      SOEM1,
      eventTime = list(sc_event),   # <-- scenario-specific event
      parms = par_vec_named
    ) %>%
      as.data.frame() %>%
      filter(time > 2022.9 & time < 2030) %>%
      mutate(
        CAD    = (X - IM) / GDP,
        GipGDP = Gip / GDP
      )
    
    out <- numeric(0)
    for (v in output_names) {
      out[paste0(v, "_min")]  <- min(sim[[v]],  na.rm = TRUE)
      out[paste0(v, "_mean")] <- mean(sim[[v]], na.rm = TRUE)
      out[paste0(v, "_max")]  <- max(sim[[v]],  na.rm = TRUE)
    }
    out
  }
  
  # Run simulations
  Y <- matrix(NA_real_, nrow = n_runs, ncol = length(output_stats))
  colnames(Y) <- output_stats
  for (i in seq_len(n_runs)) {
    par_i <- X[i, ]
    names(par_i) <- param_names
    out_i <- tryCatch(
      run_model(par_i),
      error = function(e) setNames(rep(NA_real_, length(output_stats)), output_stats)
    )
    Y[i, ] <- out_i[output_stats]
  }
  
  # Store
  all_results[[sc_name]] <- list(X = X, Y = Y)
  
  cat("Done:", sc_name, "\n")
}
#---------------------------
# 5) Summary stats per output
#---------------------------

# After the loop, pick your scenario first:
sc_name <- "DR3"
X <- all_results[[sc_name]]$X
Y <- all_results[[sc_name]]$Y

# THEN run summary_out on that Y:
summary_out <- t(apply(Y, 2, function(v) {
  c(
    mean   = mean(v, na.rm = TRUE),
    median = median(v, na.rm = TRUE),
    cv     = sd(v, na.rm = TRUE) / mean(v, na.rm = TRUE)
  )
}))
summary_out <- as.data.frame(summary_out)
summary_out$output <- rownames(summary_out)
rownames(summary_out) <- NULL


#Clean explosive runs 
valid_runs <- complete.cases(Y) & apply(Y, 1, function(r) all(is.finite(r)))

X_clean <- X[valid_runs, , drop = FALSE]
Y_clean <- Y[valid_runs, , drop = FALSE]



# Package gsoat -
x_df <- as.data.frame(X_clean)
colnames(x_df) <- make.unique(colnames(x_df))  # renames duplicate to betadfxb.1
x_df <- x_df %>%
  select(-c(lambdaicf, beta_tr, dsactive, FRactive))

# 8a) One output
toCheck<-paste(c( "inflation","fiscalDef", "firmsFrag", "unem", "GipGDP","hhFrag"),"_mean",sep="")
# j'ai enlevé hhFrag
allRes<-matrix(NA,nrow=6,ncol=243)

for(i in 1:length(toCheck)){
  name=toCheck[i]
  y_vec <- Y_clean[, name]
  
  M <- 20
  
  res <- ot_indices_1d(
    x = x_df,
    y = y_vec,
    M = M
  )
  
  summary(res, ranking = 15)
  allRes[i,]=res$indices
  plot(res, ranking = 15,main=name)
}





# 8b) Multivariate output

y_mat <- Y_clean[, paste(c( "inflation","fiscalDef", "firmsFrag", "unem", "GipGDP", "hhFrag"),"_mean",sep="")]
#j'ai aussi enlevé hhfrag

res_multi <- ot_indices(
  x = x_df,
  y = y_mat,
  M = M
)

summary(res_multi, ranking = 15)
plot(res_multi, ranking = 15)
