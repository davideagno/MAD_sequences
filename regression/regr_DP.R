
rm(list = ls())
setwd("/Users/davideagnoletto/Dropbox/sim_martingales/code/regression")

library(mvtnorm)
library(truncnorm)
library(fastmatrix)
library(coda)

Rcpp::sourceCpp("Functions_Regr.cpp")
source("Functions_Regr_ada.R")

# Initial settings
beta_true <- c(1,   -0.5, 1.5, 1, 0.5, -0.5,    -0.7, 0.5, 0.7, -0.3, -0.3)
prob_true <- c(0.5, 0.5, 0.5, 0.5, 0.5,    0.5, 0.5, 0.5, 0.5, 0.5)
cov_names <- c("X.1","X.2", "X.3", "X.4", "X.5", "X.6", "X.7", "X.8", "X.9", "X.10")
N <- 80
nSim <- 50

ll <- list(); for (i in 1:10) ll[[i]] <- 0:1; 
covs <- as.matrix(expand.grid(ll)); colnames(covs) <- NULL
mu_true <- rep(NA, NROW(covs))
for (j in 1:NROW(covs)) {
  lp <- beta_true * c(1, covs[j,])
  mu_true[j] <- exp(lp[1] + sqrt(abs(sum(lp[2:6]))) + sum(lp[7:11])^2)
}

# Data simulation
y_sim <- X_sim <- list(length = nSim)
for (s in 1:nSim) {
  set.seed(s)
  X <- matrix(nrow = N, ncol = length(prob_true))
  for (i in 1:NCOL(X)) {
    X[,i] <- rbinom(N, 1, prob_true[i])
  }
  X_mm <- model.matrix(~ X)
  bx <- X_mm * matrix(beta_true, nrow = NROW(X_mm), ncol = NCOL(X_mm), byrow = TRUE)
  mean_y <- exp(bx[,1] + sqrt(abs(rowSums(bx[,2:6]))) + rowSums(bx[,7:11])^2)
  y_sim[[s]] <- rpois(N, lambda = mean_y)
  X_sim[[s]] <- X
}


# Out-of-sample data simulation
N_new <- 10000
y_new_sim <- X_new_sim <- list(length = nSim)
for (s in 1:nSim) {
  set.seed(123 + s)
  X_new <- matrix(NA, nrow = N_new, ncol = NCOL(X))
  for (i in 1:NCOL(X_new)) {
    X_new[,i] <- rbinom(N_new, 1, prob_true[i])
  }
  colnames(X_new) <- NULL
  bx_new <- model.matrix(~ X_new) * matrix(beta_true, nrow = N_new, ncol = length(beta_true), byrow = TRUE)
  mean_y_new <- exp(bx_new[,1] + sqrt(abs(rowSums(bx_new[,2:6]))) + rowSums(bx_new[,7:11])^2)
  y_new_sim[[s]] <- rpois(N_new, lambda = mean_y_new)
  X_new_sim[[s]] <- X_new
}


# DP
mse_dp <- mae_dp <- rep(NA, nSim)
p_n_dp <- pred_dp_sim <- cev_ps_dp_sim <- list()
for (s in 1:nSim) {
  print(s)
  X <- X_sim[[s]]
  y <- y_sim[[s]]
  
  supp <- 0:100
  p_0 <- array(1/(length(supp) * 2^NCOL(X)), dim = c(length(supp), rep(2, NCOL(X))))
  alpha <- 1
  p_n_dp <- (alpha*p_0 + length(y)*rel_freq_regr(y, X, supp)) / (alpha + length(y))
  
  # Prediction
  X_new <- X_new_sim[[s]]
  y_new <- y_new_sim[[s]]
  pred_dp <- predict_mad_regr(p_n_dp, X_new)
  mse_dp[s] <- mean((y_new - pred_dp)^2)
  mae_dp[s] <- mean(abs(y_new - pred_dp))
  pred_dp_sim[[s]] <- pred_dp
  
  # UQ
  ps_dp <- pred_res_regr(p_n_dp, length(y), alpha, 1, 500, 1e-10, 1e-10, supp, 1000, 200, cores = 10, seed = 1, dp = TRUE)
  cev_ps <- matrix(nrow = NROW(covs), ncol = NCOL(ps_dp))
  for (i in 1:NROW(covs)) {
    cev_ps[i,] <- cond_exp_val_ps(ps_dp, covs[i,], supp)
  }
  cev_ps_dp_sim[[s]] <- cev_ps
  
  name_file <- paste0("SIM_regr_", N, "_DP.RData")
  save(N, y_sim, X_sim, y_new_sim, X_new_sim, p_n_dp, mse_dp, mae_dp, cev_ps_dp_sim, pred_dp_sim, 
       file = file(name_file, blocking = TRUE))
}