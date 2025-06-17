
rm(list = ls())

library(mvtnorm)
library(truncnorm)
library(fastmatrix)
library(ggplot2)
library(coda)

library(predtools)
library(magrittr)
library(dplyr)

library(parallel)
library(rstan)
library(rstanarm)
library(BART)
library(caret)
library(surfin)

Rcpp::sourceCpp("Functions_Regr.cpp")
source("Functions_Regr_ada.R")

# Initial settings
beta_true <- c(1,   -0.5, 1.5, 1, 0.5, -0.5,    -0.7, 0.5, 0.7, -0.3, -0.3)
prob_true <- c(0.5, 0.5, 0.5, 0.5, 0.5,    0.5, 0.5, 0.5, 0.5, 0.5)
cov_names <- c("X.1","X.2", "X.3", "X.4", "X.5", "X.6", "X.7", "X.8", "X.9", "X.10")
N <- 40
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


# Random Forest
mse_rf <- mae_rf <- rep(NA, nSim)
pred_rf_sim <- uq_rf <- list()
for (s in 1:nSim) {
  print(s)
  X <- X_sim[[s]]
  y <- y_sim[[s]]
  X_rf <- apply(X, 2, as.factor)
  colnames(X_rf) <- cov_names

  ctrl_mtry <- trainControl(method = "repeatedcv", number = 10, repeats = 5, search = "grid")
  rf_tune <- train(x = X_rf, y = y, method = "rf", tuneLength = 10, trControl = ctrl_mtry)
  mtry_opt <- as.numeric(rf_tune$bestTune)
  rf <- forest(X, y, var.type = "ustat", mtry = mtry_opt, ntree = 500, B = 25)

  # Prediction
  X_new <- X_new_sim[[s]]
  y_new <- y_new_sim[[s]]
  X_new_rf <- apply(X_new, 2, as.factor)
  colnames(X_new_rf) <- cov_names
  pred_rf <- predict(rf, X_new)
  mse_rf[s] <- mean((y_new - pred_rf)^2)
  mae_rf[s] <- mean(abs(y_new - pred_rf))
  pred_rf_sim[[s]] <- pred_rf

  # UQ
  pred_rf_uq <- predict(rf, covs, individualTrees = TRUE)
  ustat <- forest.varU(pred_rf_uq$predictedAll, rf)
  cv_rf <- rep(0, NROW(covs))
  for (i in 1:length(cv_rf)) {
    ci <- ustat[i,1] + c(-1,1)*qnorm(0.975)*ustat[i,2]
    if (mu_true[i] < ci[2] & mu_true[i] > ci[1]) cv_rf[i] <- 1
  }
  uq_rf[[s]] <- cv_rf
}
# save.image("~/Dropbox/sim_martingales/sim/regr/sim_regr_40_rf.RData")


# BART
mse_bart <- mae_bart <- rep(NA, nSim)
pred_bart_sim <- uq_bart <- list()
for (s in 1:nSim) {
  print(s)
  X <- X_sim[[s]]
  y <- y_sim[[s]]
  # colnames(X) <- cov_names
  data <- data.frame(X = X, y = y)

  bart <- wbart(x.train = data[,-11], y.train = data$y, nskip = 1000, ndpost = 1000)

  # Prediction
  X_new <- X_new_sim[[s]]
  y_new <- y_new_sim[[s]]
  X_new_df <- as.data.frame(X_new)
  colnames(X_new_df) <- cov_names
  pred <- predict(bart, X_new_df, mc.cores = 1)
  pred_bart <- colMeans(pred)
  mse_bart[s] <- mean((y_new - pred_bart)^2)
  mae_bart[s] <- mean(abs(y_new - pred_bart))
  pred_bart_sim[[s]] <- pred_bart

  # UQ
  pred_bart_uq <- predict(bart, covs, mc.cores = 1)
  cv_bart <- rep(0, NROW(covs))
  for (i in 1:length(cv_bart)) {
    ci <- robustBLME::hpd(pred_bart_uq[,i])
    if (mu_true[i] <= ci[2] & mu_true[i] >= ci[1]) cv_bart[i] <- 1
  }
  uq_bart[[s]] <- cv_bart
}
# save.image("~/Dropbox/sim_martingales/sim/regr/sim_regr_40_rf_bart.RData")


# Poisson regression
mse_glm <- mae_glm <- rep(NA, nSim)
pred_glm_sim <- uq_glm <- list()
for (s in 1:nSim) {
  print(s)
  X <- X_sim[[s]]
  y <- y_sim[[s]]
  data <- data.frame(X = X, y = y)

  mod <- stan_glm(y ~ ., data = data, family = poisson(link = "log"),
                  warmup = 1000, iter = 2000, chains = 4, cores = 4, seed = 1)

  # Prediction
  X_new <- X_new_sim[[s]]
  y_new <- y_new_sim[[s]]
  X_new_df <- as.data.frame(X_new)
  colnames(X_new_df) <- cov_names
  post_pred_glm <- posterior_predict(mod, newdata = X_new_df)
  pred_glm <- colMeans(post_pred_glm)
  mse_glm[s] <- mean((y_new - pred_glm)^2)
  mae_glm[s] <- mean(abs(y_new - pred_glm))
  pred_glm_sim[[s]] <- pred_glm

  # UQ
  covs_n <- data.frame(covs); colnames(covs_n) <- cov_names
  covs_post <- posterior_predict(mod, newdata = covs_n)
  cv_glm <- rep(0, NROW(covs))
  for (i in 1:length(cv_glm)) {
    ci <- robustBLME::hpd(covs_post[,i])
    if (mu_true[i] <= ci[2] & mu_true[i] >= ci[1]) cv_glm[i] <- 1
  }
  uq_glm[[s]] <- cv_glm
}
# save.image("~/Dropbox/sim_martingales/sim/regr/sim_regr_40_rf_bart_glm.RData")


