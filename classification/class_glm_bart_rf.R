
rm(list = ls())
setwd("/Users/davideagnoletto/Dropbox/sim_martingales/code")

library(mvtnorm)
library(truncnorm)
library(fastmatrix)
library(ggplot2)
library(coda)

library(predtools)
library(magrittr)
library(dplyr)
library(plotROC)
library(pROC)

library(parallel)
library(rstan)
library(rstanarm)
library(BART)
library(caret)
library(surfin)


# Initial settings
beta_true <- c(-3,   2, -4, 3, -3,   1, -3,   2, -3,  3, -2)
prob_true <- c(0.45, 0.65, 0.7, 0.4,  0.4, 0.6, 0.7, 0.3, 0.55, 0.55)
cov_names <- c("X.1","X.2", "X.3", "X.4", "X.5", "X.6", "X.7", "X.8", "X.9", "X.10")
N <- 300
nSim <- 50
pred_mad_sim <- pred_dp_sim <- pred_glm_sim <- pred_bart_sim <- pred_rf_sim <- list(length = nSim)
auc_mad <- auc_dp <- auc_glm <- auc_bart <- auc_rf <- rep(NA, nSim)
uq_mad <- uq_dp <- uq_glm <- uq_bart <- uq_rf <- list(length = nSim)

ll <- list(); for (i in 1:10) ll[[i]] <- 0:1; 
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



# Logistic regression
for (s in 1:nSim) {
  print(s)
  X <- X_sim[[s]]
  y <- y_sim[[s]]
  data <- data.frame(X = X, y = y)
  
  mod <- stan_glm(y ~ ., data = data, family = binomial(link = "logit"),
                  warmup = 1000, iter = 2000, chains = 4, cores = 4, seed = 1)
  
  # Prediction
  X_new <- X_new_sim[[s]]
  y_new <- y_new_sim[[s]]
  X_new_df <- as.data.frame(X_new)
  colnames(X_new_df) <- cov_names
  post_pred_glm <- posterior_predict(mod, newdata = X_new_df)
  pred_glm <- colMeans(post_pred_glm)
  auc_glm[s] <- auc(roc(y_new, pred_glm))
  pred_glm_sim[[s]] <- pred_glm
  
  # UQ
  covs_n <- data.frame(covs); colnames(covs_n) <- cov_names
  lp_post <- posterior_linpred(mod, newdata = covs_n)
  cv_glm <- rep(0, NROW(covs))
  for (i in 1:length(cv_glm)) {
    ci <- robustBLME::hpd(plogis(lp_post[,i]))
    if (p_true[i] <= ci[2] & p_true[i] >= ci[1]) cv_glm[i] <- 1
  }
  uq_glm[[s]] <- cv_glm
}
# save.image("~/Dropbox/sim_martingales/file_buoni/SIM_Class_150_mad_glm.RData")


# BART
for (s in 1:nSim) {
  print(s)
  X <- X_sim[[s]]
  y <- y_sim[[s]]
  data <- data.frame(X = X, y = y)
  
  bart <- lbart(x.train = data[,-11], y.train = data$y, nskip = 1000, ndpost = 1000)
  
  # Prediction
  X_new <- X_new_sim[[s]]
  y_new <- y_new_sim[[s]]
  X_new_df <- as.data.frame(X_new)
  colnames(X_new_df) <- cov_names
  pred <- predict(bart, X_new_df, mc.cores = 4)
  pred_bart <- pred$prob.test.mean
  auc_bart[s] <- auc(roc(y_new, pred_bart))
  pred_bart_sim[[s]] <- pred_bart
  
  # UQ
  pred_bart_uq <- predict(bart, covs)
  cv_bart <- rep(0, NROW(covs))
  for (i in 1:length(cv_bart)) {
    ci <- robustBLME::hpd(pred_bart_uq$prob.test[,i])
    if (p_true[i] <= ci[2] & p_true[i] >= ci[1]) cv_bart[i] <- 1
  }
  uq_bart[[s]] <- cv_bart
}
# save.image("~/Dropbox/sim_martingales/file_buoni/SIM_Class_150_mad_glm_bart.RData")


# Random Forest
for (s in 1:nSim) {
  print(s)
  X <- X_sim[[s]]
  y <- y_sim[[s]]
  X_rf <- apply(X, 2, as.factor)
  colnames(X_rf) <- cov_names
  
  ctrl_mtry <- trainControl(method = "repeatedcv", number = 10, repeats = 5, search = "grid")
  rf_tune <- train(x = X_rf, y = as.factor(y), method = "rf", tuneLength = 10, trControl = ctrl_mtry)
  mtry_best <- as.numeric(rf_tune$bestTune)
  rf <- forest(X, y, var.type = "ustat", mtry = mtry_best, ntree = 500, B = 25)
  
  # Prediction
  X_new <- X_new_sim[[s]]
  y_new <- y_new_sim[[s]]
  pred_rf <- predict(rf, X_new)
  auc_rf[s] <- auc(roc(y_new, pred_rf))
  pred_rf_sim[[s]] <- pred_rf
  
  # UQ
  pred_rf_uq <- predict(rf, covs, individualTrees = TRUE)
  ustat <- forest.varU(pred_rf_uq$predictedAll, rf)
  cv_rf <- rep(0, NROW(covs)) 
  for (i in 1:length(cv_rf)) {
    ci <- ustat[i,1] + c(-1,1)*qnorm(0.975)*ustat[i,2]
    if (p_true[i] < ci[2] & p_true[i] > ci[1]) cv_rf[i] <- 1
  }
  uq_rf[[s]] <- cv_rf
}
# save.image("~/Dropbox/sim_martingales/file_buoni/SIM_Class_150_mad_glm_bart_rf.RData")

