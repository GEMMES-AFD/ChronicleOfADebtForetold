library(lhs)
library(dplyr)
library(gsaot)
library(sensobol)
library(dplyr)

source("01_run_scenarios.R")

n_runs <- 5000

output_names <- c(
  "perCapita", "inflation", "reserves", "foreignDebt", "privateDebt", "pubDebt",
  "fiscalDef", "hhFrag", "firmsFrag", "unem", "CAD", "GipGDP"
)

stopifnot(is.numeric(parms1), length(parms1) == 246, !is.null(names(parms1)))

param_names <- names(parms1)
baseline <- as.numeric(parms1)
names(baseline)<-param_names

# 2) LHS Sampling ----

U <- randomLHS(n_runs, length(baseline)) 
lower_raw <- baseline * 0.75
upper_raw <- baseline * 1.5

lower <- pmin(lower_raw, upper_raw)
upper <- pmax(lower_raw, upper_raw)

# handle baseline == 0 (otherwise no variation)
eps0 <- 1e-6
idx0 <- baseline == 0
lower[idx0] <- -eps0
upper[idx0] <-  eps0


# Build X
X <- sweep(U, 2, (upper - lower), `*`)
X <- sweep(X, 2, lower, `+`)
colnames(X) <- param_names


constantParms<-c("lambdaicf","beta_tr","dsactive","FRactive","gamma_tr","alpha_tr","delta_tr","lambdatr0_adj","lambdatr1_adj","lambdatr2_adj","shrGrTax","shrGrIC","shrDon","shrGrBw","shrGrL","shrGrLFx","betadebtSwapFXLgFX","betadebtSwapFXLgFXtr","betadebtSwapFXBgFX","md_lgtr","md_bgtr","K_0","lambdatr0","lambdatr1","lambdatr2","taumtr","tauothktr","betasigmamktr","zetafx2")
for(constPar in constantParms){
  X[,constPar]=baseline[constPar]
}

# 3) Model runner ----

run_model <- function(par_vec_named,bas = NULL) {
  
  sim <- cppRK4(
    SOEM1,
    eventTime = list(event1),
    parms = par_vec_named
  ) %>%
    as.data.frame() %>%
    filter(time > 2022.9 & time < 2030)%>%
    mutate(
      CAD    = (X - IM) / GDP,
      GipGDP = Gip / GDP
    )
  
  out <- numeric(0)
  
  for (v in output_names) {
    out[paste0(v, "_min")]  <- min(sim[[v]], na.rm = TRUE)
    out[paste0(v, "_mean")] <- mean(sim[[v]], na.rm = TRUE)
    out[paste0(v, "_max")]  <- max(sim[[v]], na.rm = TRUE)
  }
  if (!is.null(bas)) {
    eps <- 1e-5
    dist_sq <- 0
    for (v in output_names) {
      sim_v <- out[paste0(v, "_mean")]
      bas_v <- bas[paste0(v, "_mean")]
      dist_sq <- dist_sq + ((sim_v - bas_v) / (abs(bas_v) + eps))^2
    }
    out["dist_baseline"] <- sqrt(dist_sq)  # Euclidean in relative space
  }
  
  return(out)
}

# 4) Run simulations ----

output_stats <- c(
  as.vector(t(outer(output_names, c("min", "mean", "max"), paste, sep = "_"))),
  "dist_baseline"
)

Y <- matrix(NA_real_, nrow = n_runs, ncol = length(output_stats))
colnames(Y) <- output_stats

# Compute baseline output once
bas_out <- run_model(setNames(baseline, param_names))

for (i in seq_len(n_runs)) {
  par_i <- X[i, ]
  names(par_i) <- param_names
  
  out_i <- tryCatch(
    run_model(par_i, bas = bas_out),
    error = function(e) setNames(rep(NA_real_, length(output_stats)), output_stats)
  )
  
  Y[i, ] <- out_i[output_stats]
}

# 5) Summary stats per output ----

summary_out <- t(apply(Y, 2, function(v) {
  c(
    mean   = mean(v, na.rm = TRUE),
    median = median(v, na.rm = TRUE),
    cv    = sd(v, na.rm = TRUE)/mean(v, na.rm = TRUE)
  )
}))

summary_out <- as.data.frame(summary_out)
summary_out$output <- rownames(summary_out)
rownames(summary_out) <- NULL

# Results:
dim(X)
dim(Y)  
summary_out
summary_out[order(summary_out$cv),]

#clean 

valid_runs <- apply(Y, 1, function(r) all(is.finite(r) & abs(r)<1e3 ))
X_clean <- X[valid_runs, , drop = FALSE]
Y_clean <- Y[valid_runs, , drop = FALSE]
Y <- cbind(Y, valid = as.integer(valid_runs))

# Package gsoat ----

x_df_clean <- as.data.frame(X_clean)

# 8a) One output
toCheck<-paste(c( "inflation","fiscalDef", "firmsFrag", "unem", "GipGDP","hhFrag","perCapita","CAD", "foreignDebt"),"_mean",sep="")
allRes<-matrix(NA,nrow=length(toCheck),ncol=length(x_df_clean))
rownames(allRes) <- toCheck

for(i in 1:length(toCheck)){
  name=toCheck[i]
  y_vec <- Y_clean[, name]
  
  M <- 20
  
  res <- ot_indices_1d(
    x = x_df_clean,
    y = y_vec,
    M = M
  )
  
  summary(res, ranking = 15)
  allRes[i,]=res$indices
  plot(res, ranking = 15,main=name)
}

#paramètres responsables de crash

M <- 20

x_df_all <- as.data.frame(X)
colnames(x_df_all) <- make.unique(colnames(x_df_all))  # renames duplicate to betadfxb.1
  res <- ot_indices_1d(
    x = x_df_all,
    y = Y[,"valid"],
    M = M
  )
  
  summary(res, ranking = 15)
  plot(res, ranking = 15,main=name)


# 8b) Multivariate output

y_mat <- Y_clean[, toCheck]

res_multi <- ot_indices(
  x = x_df_clean,
  y = y_mat,
  M = M,
)

summary(res_multi, ranking = 15)
plot(res_multi, ranking = 15)


# Sensitivity analysis on distance baseline 

y_dist <- Y_clean[, "dist_baseline"]

res_dist <- ot_indices_1d(x = x_df_clean, y = y_dist, M = 20)
summary(res_dist, ranking = 15)
plot(res_dist, ranking = 15, main = "Parameters driving deviation from baseline")


# Filter the runs that are closer to the baseline

threshold <- quantile(y_dist, 0.10)  # 10% les plus proches
similar_runs <- y_dist <= threshold

X_similar <- X_clean[similar_runs, ]
Y_similar <- Y_clean[similar_runs, ]

cat("Runs- Short distance to the baseline", sum(similar_runs), "\n")


# Compare output distribution 

par(mfrow = c(4, 3))
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "",
     main = "Output distributions: All runs vs Similar runs")  # ← global title

for (v in c("inflation", "fiscalDef", "firmsFrag", "unem", "GipGDP", "hhFrag","privateDebt","CAD")) {
  col <- paste0(v, "_mean")
  
  plot(density(Y_clean[, col],   na.rm = TRUE), col = "grey70", lwd = 2,
       main = v, xlab = "")
  lines(density(Y_similar[, col], na.rm = TRUE), col = "steelblue", lwd = 2)
  abline(v = bas_out[col], col = "red", lty = 2, lwd = 2)
  legend("topright", legend = c("All runs", "Similar runs", "Baseline"),
         col = c("grey70", "steelblue", "red"), lwd = 2, lty = c(1,1,2), cex = 0.7)
}

#Values of parameters in similar runs 

top_params <- names(sort(res_dist$indices, decreasing = TRUE))[1:10]

baseline_named <- setNames(baseline, param_names)

par(mfrow = c(2, 5))
for (p in top_params) {
  plot(density(X_clean[, p]),    col = "grey70",    lwd = 2, main = p, xlab = "")
  lines(density(X_similar[, p]), col = "steelblue", lwd = 2)
  abline(v = baseline_named[p],  col = "red", lty = 2, lwd = 2)
  legend("topright", legend = c("All runs", "Similar runs", "Baseline"),
         col = c("grey70", "steelblue", "red"), lwd = 2, lty = c(1, 1, 2), cex = 0.6)
}


#topparams: 

top_params_gsaot <- names(sort(res_dist$indices, decreasing = TRUE))[1:30]
newConstantParms <- param_names[!param_names %in% top_params_gsaot]
  
  
top_params_sensobol <- c(
  # Très robustes (2+ sources)
  "mu0", "sigmapic", "sigmapc", "omegaf3", "betapremf", "zeta0",
  "mpcUB", "betaiwst", "rho0", "ipsilon0w",
  # Multivariate fort
  "rhofx1", "sigmaRem", "sigmaG1", "fi5", "LBFFX",
  "rho2", "betasigmamk", "varsigmafdi2", "tauf",
  # Crash critiques
  "mpcLB", "phisc", "v2", "alphapw", "chi0",
  "betariskFFX", "betaen", "sigmaRfxb",
  # Univariate robustes restants
  "sigmaxnp", "etadbfx", "kappa03"
  )

#pour remplacer dans Sensobol si on veut éviter que ce soit trop lourd 


# Sensobol

param_names_unique <- make.unique(param_names)

N <-5000

# sobol_matrices expects a character vector
mat <- sobol_matrices(N = N, params = top_params_sensobol)
# mat is [0,1], rescale manually
for (j in seq_len(length(top_params_sensobol))) {
  mat[, j] <- lower[top_params_sensobol[j]] + mat[, j] * (upper[top_params_sensobol[j]] - lower[top_params_sensobol[j]])
}

# Force constant parameters to baseline
#for (constPar in constantParms) {
#  mat[, constPar] <- baseline[constPar]
#}

# Run the model
n_total <- nrow(mat)

Y_sobol <- matrix(NA_real_, nrow = n_total, ncol = length(output_stats))
colnames(Y_sobol) <- output_stats

for (i in seq_len(n_total)) {
  par_i <- baseline                          # vecteur complet baseline (247 params)
  par_i[top_params_sensobol] <- mat[i, ]    # override uniquement les 30 sélectionnés
  names(par_i) <- param_names
  Y_sobol[i, ] <- tryCatch(
    run_model(par_i, bas = bas_out),
    error = function(e) setNames(rep(NA_real_, length(output_stats)), output_stats)
  )
}

# Compute Sobol indices per output
results_list <- list()

for (name in toCheck) {
  y_vec <- Y_sobol[, name]
  
  # sobol_indices ne tolère pas les NA
  if (sum(is.finite(y_vec)) < n_total * 0.9) {
    cat("Too many NAs for", name, "— skipping\n")
    next
  }
  
  # Remplace les rares NA par la médiane (ou stoppe si trop)
  y_vec[!is.finite(y_vec)] <- median(y_vec, na.rm = TRUE)
  
  ind <- tryCatch(
    sobol_indices(
      Y      = y_vec,
      N      = N,
      params = top_params_sensobol,
      order  = "first",
      boot   = TRUE,      # intervalles de confiance
      R      = 100        # bootstrap replicates
    ),
    error = function(e) { cat("Error for", name, ":", e$message, "\n"); NULL }
  )
  
  results_list[[name]] <- ind
  
  if (!is.null(ind)) {
    p <- plot(ind, ranking = 15) + ggtitle(name)
    print(p)
  }
}
