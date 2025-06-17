
rm(list = ls())
setwd("/Users/davideagnoletto/Dropbox/sim_martingales/code/regression")

library(mvtnorm)
library(truncnorm)
library(fastmatrix)
library(parallel)
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)

Rcpp::sourceCpp("Functions_Regr.cpp")
source("Functions_Regr_1.R")

# Initial settings
beta_true <- c(1,   -0.5, 1.5, 1, 0.5, -0.5,    -0.7, 0.5, 0.7, -0.3, -0.3)
prob_true <- c(0.5, 0.5, 0.5, 0.5, 0.5,    0.5, 0.5, 0.5, 0.5, 0.5)
cov_names <- c("X.1","X.2", "X.3", "X.4", "X.5", "X.6", "X.7", "X.8", "X.9", "X.10")
N <- 80
cores_mad <- 20
nSim <- 50

ll <- list(); for (i in 1:10) ll[[i]] <- 0:1; 
covs <- as.matrix(expand.grid(ll)); colnames(covs) <- NULL
mu_true <- rep(NA, NROW(covs))
for (j in 1:NROW(covs)) {
  lp <- beta_true * c(1, covs[j,])
  mu_true[j] <- exp(lp[1] + sqrt(abs(sum(lp[2:6]))) + sum(lp[7:11])^2)
}


# Data simulation
print("Data simulation")
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
print("Out-of-sample data simulation")
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


# MAD
print("MAD sequence")
B <- 200
sigma_opt <- mse_mad <- mae_mad <- rep(NA, nSim)
p_n_mad <- pred_mad_sim <- cev_ps_mad_sim <- list()
for (s in 17:nSim) {
  print(s)
  X <- X_sim[[s]]
  y <- y_sim[[s]]
  
  supp <- 0:100
  p_0 <- array(1/(length(supp) * 2^NCOL(X)), dim = c(length(supp), rep(2, NCOL(X))))
  alpha <- 1
  
  opt <- nlminb(start = 3.5, lower = 1e-10, upper = Inf,
                objective = function(h) -preq_cond_regr(y, X, alpha, p_0, c(h, 0.25), supp, 10, 10),
                control = list(rel.tol = 1e-2))
  print(opt$message)
  sigma_opt[s] <- opt$par
  p_n <- pred_regr(y, X, alpha, p_0, c(sigma_opt[s], 0.25), supp, 10, cores_mad)
  p_n_mad[[s]] <- p_n
  print("MAD - predictive done")
  
  # Prediction
  X_new <- X_new_sim[[s]]
  y_new <- y_new_sim[[s]]
  pred <- predict_mad_regr(p_n, X_new)
  mse_mad[s] <- mean((y_new - pred)^2)
  mae_mad[s] <- mean(abs(y_new - pred))
  pred_mad_sim[[s]] <- pred
  
  # UQ 
  ps <- pred_res_regr(p_n, N, alpha, sigma_opt[s], 0.25, supp, 1000, B, cores_mad, 1)
  print("MAD - resampling done")
  cev_ps <- matrix(nrow = NROW(covs), ncol = NCOL(ps))
  for (i in 1:NROW(covs)) {
    cev_ps[i,] <- cond_exp_val_ps(ps, covs[i,], supp)
  }
  cev_ps_mad_sim[[s]] <- cev_ps
  
  name_file <- paste0("SIM_regr_", N, "_1_cond.RData")
  save(N, y_sim, X_sim, y_new_sim, X_new_sim, p_n_mad, sigma_opt,
       mse_mad, mae_mad, cev_ps_mad_sim, pred_mad_sim, 
       file = file(name_file, blocking = TRUE))
}