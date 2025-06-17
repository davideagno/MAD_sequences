
rm(list = ls())

library(mvtnorm)
library(truncnorm)
library(fastmatrix)
library(parallel)
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(pROC)

Rcpp::sourceCpp("Functions_Class.cpp")
source("Functions_Class_lambda.R")

# Initial settings
beta_true <- c(-3,   2, -4, 3, -3,   1, -3,   2, -3,  3, -2)
prob_true <- c(0.45, 0.65, 0.7, 0.4,  0.4, 0.6, 0.7, 0.3, 0.55, 0.55)
cov_names <- c("X.1","X.2", "X.3", "X.4", "X.5", "X.6", "X.7", "X.8", "X.9", "X.10")
N <- 150
nSim <- 50
cores_mad <- 20

ll <- list(); for (i in 1:10) ll[[i]] <- 0:1
covs <- as.matrix(expand.grid(ll)); colnames(covs) <- NULL
p_true <- rep(NA, NROW(covs))
for (j in 1:NROW(covs)) {
  lp <- beta_true * c(1, covs[j,])
  p_true[j] <- plogis(sum(lp[1:5]) + lp[6]*lp[7] + sqrt(abs(sum(lp[8:9]))) + sum(lp[10:11])^2)
}

# Data simulation
y_sim <- X_sim <- list(length = nSim)
for (s in 1:nSim) {
  set.seed(110 + s)
  X <- matrix(nrow = N, ncol = length(prob_true))
  for (i in 1:NCOL(X)) {
    X[,i] <- rbinom(N, 1, prob_true[i])
  }
  bx <- model.matrix(~ X) * matrix(beta_true, nrow = N, ncol = NCOL(X)+1, byrow = TRUE)
  mean_y <- plogis(bx[,1] + rowSums(bx[,2:5]) + bx[,6]*bx[,7] + sqrt(abs(rowSums(bx[,8:9]))) + rowSums(bx[,10:11])^2)
  y_sim[[s]] <- rbinom(N, 1, prob = mean_y)
  X_sim[[s]] <- X
}

# Out-of-sample data simulation
N_new <- 10000
y_new_sim <- X_new_sim <- list(length = nSim)
for (s in 1:nSim) {
  set.seed(nSim + s)
  X_new <- matrix(nrow = N_new, ncol = NCOL(X))
  for (i in 1:NCOL(X_new)) {
    X_new[,i] <- rbinom(N_new, 1, prob_true[i])
  }
  colnames(X_new) <- NULL
  bx_new <- model.matrix(~ X_new) * matrix(beta_true, nrow = N_new, ncol = length(beta_true), byrow = TRUE)
  mean_y_new <- plogis(bx_new[,1] + rowSums(bx_new[,2:5]) + bx_new[,6]*bx_new[,7] + sqrt(abs(rowSums(bx_new[,8:9]))) + rowSums(bx_new[,10:11])^2)
  X_new_df <- as.data.frame(X_new)
  colnames(X_new_df) <- cov_names
  y_new_sim[[s]] <- rbinom(N_new, 1, prob = mean_y_new)
  X_new_sim[[s]] <- X_new
}


# DP
B <- 200
delta_opt <- auc_mad <- rep(NA, nSim)
pred_mad_sim <- uq_mad_pc_ps <- list(length = nSim)
for (s in 1:nSim) {
  print(s)
  X <- X_sim[[s]]
  y <- y_sim[[s]]
  p_0 <- array(1/(2^(NCOL(X)+1)), dim = rep(2, NCOL(X)+1))
  alpha <- 1
  p_n <- (alpha*p_0 + N*rel_freq_class(y, X)) / (alpha + N)
  
  # Prediction
  X_new <- X_new_sim[[s]]
  y_new <- y_new_sim[[s]]
  pred <- predict_mad_class(p_n, X_new)
  auc_mad[s] <- auc(roc(y_new, pred))
  pred_mad_sim[[s]] <- pred
  print("predictive done")
  
  # UQ
  ps <- pred_res_class(p_n, N, alpha, 1, c(1e-10, 1e-10), 1000, B, cores_mad, 1)
  print("resampling done")
  pc_ps <-  matrix(nrow = NROW(covs), ncol = B)
  for (b in 1:B) {
    for (j in 1:NROW(pc_ps)) {
      i_1 <- sum(c(1, covs[j,]) * 2^c(0:10)) + 1
      pc_ps[j, b] <- ps[i_1, b] / (ps[i_1 - 1, b] + ps[i_1, b])
    }
  }
  uq_mad_pc_ps[[s]] <- pc_ps
  
  name_file <- paste0("SIM_class_", N, "_DP.RData")
  save(auc_mad, pred_mad_sim, uq_mad_pc_ps, delta_opt,
       file = file(name_file, blocking = TRUE))
}
