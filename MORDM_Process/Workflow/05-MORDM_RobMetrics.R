if (SHOW_PROGRESS) {
  library(progressr)
  handlers("txtprogressbar")
}

for (shape in shapes){

  cat("Calling robustness file...\n")
  roblist <- readRDS(ROBLIST_FILE)
  maxx <- lapply(roblist, function(x){
    y <- x[-which(x[,'Failure'] == T),7:18]
    apply(y, 2, max, na.rm=T)*1.01
  }) %>% list.rbind()
  maxx <- apply(maxx, 2, max, na.rm=T)

  rob_pol <- c()

  cat("Attributing Failure Scores...\n")
  roblist2 <- lapply(roblist, function(x){
    for (i in which(x[,'Failure'] == T)) x[i,7:18] <- maxx
    x
  })

  cat("Setting up clustering...\n")
  plan(cluster)
  set <- mordm.get.set(readRDS(OPTPOL_FILE))

  if (COMPUTE_REG1){
    cat("Computing Regret I...\n")

    .reg1_fn <- function(ll) {
      a <- roblist2[[ll]] %>% as.data.frame()
      max(sapply(1:ncol(a[,7:18]), function(j){
        quantile(sapply(1:nrow(a), function(k){
          (a[k,j+6] - set[ll,j+6])/abs(set[ll,j+6])
        }), probs=0.9, na.rm=T)
      }))
    }

    if (SHOW_PROGRESS) {
      with_progress({
        p <- progressor(steps=length(roblist2))
        regret_type1 <- future.apply::future_lapply(1:length(roblist2), function(ll){
          on.exit(p())
          .reg1_fn(ll)
        }) %>% list.rbind()
      })
    } else {
      regret_type1 <- future.apply::future_lapply(1:length(roblist2), .reg1_fn) %>% list.rbind()
    }

    regret_type1 <- cbind(set, regret_type1)
    colnames(regret_type1)[ncol(regret_type1)] <- "Regret"
    rob_pol[1] <- which.min(regret_type1[,"Regret"])
    saveRDS(regret_type1, REG1_FILE)
  }

  if (COMPUTE_REG2){
    cat("Computing Regret II...\n")
    start_time <- Sys.time()

    .reg2_fn <- function(ll) {
      max(sapply(7:18, function(j){
        quantile(sapply(1:nrow(roblist2[[ll]]), function(i){
          best <- min(sapply(1:length(roblist2), function(k) roblist2[[k]][i,j]), na.rm=T)
          (roblist2[[ll]][i,j] - best)/abs(best)
        }), 0.9, na.rm=T)
      }))
    }

    if (SHOW_PROGRESS) {
      with_progress({
        p <- progressor(steps=length(roblist2))
        regret_type2 <- future.apply::future_lapply(1:length(roblist2), function(ll){
          on.exit(p())
          .reg2_fn(ll)
        }) %>% list.rbind()
      })
    } else {
      regret_type2 <- future.apply::future_lapply(1:length(roblist2), .reg2_fn) %>% list.rbind()
    }

    cat("Regret II elapsed:", format(Sys.time() - start_time), "\n")

    regret_type2 <- cbind(set, regret_type2)
    colnames(regret_type2)[ncol(regret_type2)] <- "Regret"
    rob_pol[2] <- which.min(regret_type2[,"Regret"])
    saveRDS(regret_type2, REG2_FILE)
  }
  
  if (COMPUTE_REG2SD){
    cat("Computing Regret II...\n")
    start_time <- Sys.time()
    
    .reg2_fn <- function(ll) {
      max(sapply(7:18, function(j){
        quantile(sapply(1:nrow(roblist2[[ll]]), function(i){
          best <- min(sapply(1:length(roblist2), function(k) roblist2[[k]][i,j]), na.rm=T)
          sd_SOW <- sd(sapply(1:length(roblist2), function(k){
            roblist2[[k]][i,j]
          }), na.rm = T)
          (roblist2[[ll]][i,j] - best)/sd_SOW
        }), 0.9, na.rm=T)
      }))
    }
    
    if (SHOW_PROGRESS) {
      with_progress({
        p <- progressor(steps=length(roblist2))
        regret_type2 <- future.apply::future_lapply(1:length(roblist2), function(ll){
          on.exit(p())
          .reg2_fn(ll)
        }) %>% list.rbind()
      })
    } else {
      regret_type2 <- future.apply::future_lapply(1:length(roblist2), .reg2_fn) %>% list.rbind()
    }
    
    cat("Regret II elapsed:", format(Sys.time() - start_time), "\n")
    
    regret_type2 <- cbind(set, regret_type2)
    colnames(regret_type2)[ncol(regret_type2)] <- "Regret"
    rob_pol[2] <- which.min(regret_type2[,"Regret"])
    saveRDS(regret_type2, REG2SD_FILE)
  }

  if (COMPUTE_SAT1 | COMPUTE_SAT2){
    roblist2_df <- list.rbind(roblist2)
    quantile_rob <- apply(roblist2_df[,5:16], 2, function(x) quantile(x, 0.15))
  }

  if (COMPUTE_SAT1){
    cat("Computing Satisfying I...\n")
    if (SHOW_PROGRESS) pb <- txtProgressBar(min=0, max=length(roblist2), style=3)
    satisficing_type1 <- lapply(1:length(roblist2), function(ll){
      if (SHOW_PROGRESS) setTxtProgressBar(pb, ll)
      a <- roblist2[[ll]][,7:18]
      b <- apply(a, 2, function(x, quantile_rob) x < quantile_rob, quantile_rob=quantile_rob)
      d <- apply(b, 1, sum, na.rm=T)
      length(d[which(d==11)])/length(d)
    }) %>% list.rbind()
    if (SHOW_PROGRESS) close(pb)

    satisficing_type1 <- cbind(set, satisficing_type1)
    colnames(satisficing_type1)[ncol(satisficing_type1)] <- "Satisficing"
    rob_pol[3] <- which.max(satisficing_type1[,"Satisficing"])
    saveRDS(satisficing_type1, SAT1_FILE)
  }

  if (COMPUTE_SAT2){
    cat("Computing Satisfying II...\n")
    params <- readRDS(ROBCAL_FILE)[[1]]
    dist <- apply(params[,1:6], 1, function(x){
      refs <- c(parms_NewC[c("betaen","v1","alphagw","alphapw","zetafx3","rho")], 0.025)
      for (i in 1:6){
        x[i]    <- (x[i]    - min(params[,i]))/(max(params[,i]) - min(params[,i]))
        refs[i] <- (refs[i] - min(params[,i]))/(max(params[,i]) - min(params[,i]))
      }
      norm(x - refs, "2")
    })

    if (SHOW_PROGRESS) pb <- txtProgressBar(min=0, max=length(roblist2), style=3)
    satisficing_type2 <- lapply(1:length(roblist2), function(ll){
      if (SHOW_PROGRESS) setTxtProgressBar(pb, ll)
      a <- roblist2[[ll]][,7:18]
      b <- apply(a, 2, function(x, quantile_rob) x > quantile_rob, quantile_rob=quantile_rob)
      b <- apply(b, 1, sum, na.rm=T)
      f <- min(dist[b>0], na.rm=T)
      f[which(is.infinite(f))] <- NA
      f
    }) %>% list.rbind()
    if (SHOW_PROGRESS) close(pb)

    satisficing_type2 <- cbind(set, satisficing_type2)
    colnames(satisficing_type2)[ncol(satisficing_type2)] <- "Satisficing"
    rob_pol[4] <- which.max(satisficing_type2[,"Satisficing"])
    saveRDS(satisficing_type2, SAT2_FILE)
  }

  saveRDS(rob_pol, ROBPOL_FILE)
}
