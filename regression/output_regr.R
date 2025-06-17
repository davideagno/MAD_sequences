
rm(list = ls())
setwd("/Users/davideagnoletto/Dropbox/sim_martingales/code/regression")

library(ggplot2)

load("sim_regr_40_rf_bart_glm.RData")
mse_glm_40 <- mse_glm; mse_bart_40 <- mse_bart; mse_rf_40 <- mse_rf
cv_rf_40 <- rowMeans(matrix(unlist(uq_rf), ncol = nSim, nrow = NROW(covs)))
cv_bart_40 <- rowMeans(matrix(unlist(uq_bart), ncol = nSim, nrow = NROW(covs)))
cv_glm_40 <- rowMeans(matrix(unlist(uq_glm), ncol = nSim, nrow = NROW(covs)))

load("SIM_regr_40_a34_cond.RData")
mse_mad_40_a34 <- mse_mad
uq_mad_40_a34 <- matrix(0, nrow = NROW(covs), ncol = length(p_n_mad))
for (s in 1:length(p_n_mad)) {
  cev_mad_ps <- cev_ps_mad_sim[[s]]
  for (i in 1:NROW(cev_mad_ps)) {
    ci <- robustBLME::hpd(cev_mad_ps[i,])
    if (mu_true[i] <= ci[2] & mu_true[i] >= ci[1]) uq_mad_40_a34[i,s] <- 1
  }
}
cv_mad_40_a34 <- rowMeans(matrix(unlist(uq_mad_40_a34), ncol = length(p_n_mad), nrow = NROW(covs)))

load("SIM_regr_40_23_cond.RData")
mse_mad_40_23 <- mse_mad
uq_mad_40_23 <- matrix(0, nrow = NROW(covs), ncol = length(p_n_mad))
for (s in 1:length(p_n_mad)) {
  cev_mad_ps <- cev_ps_mad_sim[[s]]
  for (i in 1:NROW(cev_mad_ps)) {
    ci <- robustBLME::hpd(cev_mad_ps[i,])
    if (mu_true[i] <= ci[2] & mu_true[i] >= ci[1]) uq_mad_40_23[i,s] <- 1
  }
}
cv_mad_40_23 <- rowMeans(matrix(unlist(uq_mad_40_23), ncol = length(p_n_mad), nrow = NROW(covs)))

load("SIM_regr_40_1_cond.RData")
mse_mad_40_1 <- mse_mad
uq_mad_40_1 <- matrix(0, nrow = NROW(covs), ncol = length(p_n_mad))
for (s in 1:length(p_n_mad)) {
  cev_mad_ps <- cev_ps_mad_sim[[s]]
  for (i in 1:NROW(cev_mad_ps)) {
    ci <- robustBLME::hpd(cev_mad_ps[i,])
    if (mu_true[i] <= ci[2] & mu_true[i] >= ci[1]) uq_mad_40_1[i,s] <- 1
  }
}
cv_mad_40_1 <- rowMeans(matrix(unlist(uq_mad_40_1), ncol = length(p_n_mad), nrow = NROW(covs)))

load("SIM_regr_40_dpm_cond.RData")
mse_mad_40_dpm <- mse_mad
uq_mad_40_dpm <- matrix(0, nrow = NROW(covs), ncol = length(p_n_mad))
for (s in 1:length(p_n_mad)) {
  cev_mad_ps <- cev_ps_mad_sim[[s]]
  for (i in 1:NROW(cev_mad_ps)) {
    ci <- robustBLME::hpd(cev_mad_ps[i,])
    if (mu_true[i] <= ci[2] & mu_true[i] >= ci[1]) uq_mad_40_dpm[i,s] <- 1
  }
}
cv_mad_40_dpm <- rowMeans(matrix(unlist(uq_mad_40_dpm), ncol = length(p_n_mad), nrow = NROW(covs)))

load("SIM_regr_40_DP.RData")
mse_dp_40 <- mse_dp
uq_dp_40 <- matrix(0, nrow = NROW(covs), ncol = length(p_n_dp))
for (s in 1:length(p_n_mad)) {
  cev_dp_ps <- cev_ps_dp_sim[[s]]
  for (i in 1:NROW(cev_dp_ps)) {
    ci <- robustBLME::hpd(cev_dp_ps[i,])
    if (mu_true[i] <= ci[2] & mu_true[i] >= ci[1]) uq_dp_40[i,s] <- 1
  }
}
cv_dp_40 <- rowMeans(matrix(unlist(uq_dp_40), ncol = length(p_n_dp), nrow = NROW(covs)))

mean(mse_glm_40); sd(mse_glm_40)
mean(mse_bart_40); sd(mse_bart_40)
mean(mse_rf_40); sd(mse_rf_40)
mean(mse_dp_40); sd(mse_dp_40)
mean(mse_mad_40_1); sd(mse_mad_40_1)
mean(mse_mad_40_23); sd(mse_mad_40_23)
mean(mse_mad_40_dpm); sd(mse_mad_40_dpm)
mean(mse_mad_40_a34); sd(mse_mad_40_a34)
boxplot(cv_glm_40, cv_bart_40, cv_rf_40, cv_dp_40, cv_mad_40_1, cv_mad_40_23, cv_mad_40_dpm, cv_mad_40_a34); abline(h = 0.95, col = 2)





load("sim_regr_80_rf_bart_glm.RData")
mse_glm_80 <- mse_glm; mse_bart_80 <- mse_bart; mse_rf_80 <- mse_rf
cv_rf_80 <- rowMeans(matrix(unlist(uq_rf), ncol = nSim, nrow = NROW(covs)))
cv_bart_80 <- rowMeans(matrix(unlist(uq_bart), ncol = nSim, nrow = NROW(covs)))
cv_glm_80 <- rowMeans(matrix(unlist(uq_glm), ncol = nSim, nrow = NROW(covs)))

load("SIM_regr_80_a34_cond.RData")
mse_mad_80_a34 <- mse_mad
uq_mad_80_a34 <- matrix(0, nrow = NROW(covs), ncol = length(p_n_mad))
for (s in 1:length(p_n_mad)) {
  cev_mad_ps <- cev_ps_mad_sim[[s]]
  for (i in 1:NROW(cev_mad_ps)) {
    ci <- robustBLME::hpd(cev_mad_ps[i,])
    if (mu_true[i] <= ci[2] & mu_true[i] >= ci[1]) uq_mad_80_a34[i,s] <- 1
  }
}
cv_mad_80_a34 <- rowMeans(matrix(unlist(uq_mad_80_a34), ncol = length(p_n_mad), nrow = NROW(covs)))

load("SIM_regr_80_23_cond.RData")
mse_mad_80_23 <- mse_mad
uq_mad_80_23 <- matrix(0, nrow = NROW(covs), ncol = length(p_n_mad))
for (s in 1:length(p_n_mad)) {
  cev_mad_ps <- cev_ps_mad_sim[[s]]
  for (i in 1:NROW(cev_mad_ps)) {
    ci <- robustBLME::hpd(cev_mad_ps[i,])
    if (mu_true[i] <= ci[2] & mu_true[i] >= ci[1]) uq_mad_80_23[i,s] <- 1
  }
}
cv_mad_80_23 <- rowMeans(matrix(unlist(uq_mad_80_23), ncol = length(p_n_mad), nrow = NROW(covs)))

load("SIM_regr_80_1_cond.RData")
mse_mad_80_1 <- mse_mad
uq_mad_80_1 <- matrix(0, nrow = NROW(covs), ncol = length(p_n_mad))
for (s in 1:length(p_n_mad)) {
  cev_mad_ps <- cev_ps_mad_sim[[s]]
  for (i in 1:NROW(cev_mad_ps)) {
    ci <- robustBLME::hpd(cev_mad_ps[i,])
    if (mu_true[i] <= ci[2] & mu_true[i] >= ci[1]) uq_mad_80_1[i,s] <- 1
  }
}
cv_mad_80_1 <- rowMeans(matrix(unlist(uq_mad_80_1), ncol = length(p_n_mad), nrow = NROW(covs)))

load("SIM_regr_80_dpm_cond.RData")
mse_mad_80_dpm <- mse_mad
uq_mad_80_dpm <- matrix(0, nrow = NROW(covs), ncol = length(p_n_mad))
for (s in 1:length(p_n_mad)) {
  cev_mad_ps <- cev_ps_mad_sim[[s]]
  for (i in 1:NROW(cev_mad_ps)) {
    ci <- robustBLME::hpd(cev_mad_ps[i,])
    if (mu_true[i] <= ci[2] & mu_true[i] >= ci[1]) uq_mad_80_dpm[i,s] <- 1
  }
}
cv_mad_80_dpm <- rowMeans(matrix(unlist(uq_mad_80_dpm), ncol = length(p_n_mad), nrow = NROW(covs)))

load("SIM_regr_80_DP.RData")
mse_dp_80 <- mse_dp
uq_dp_80 <- matrix(0, nrow = NROW(covs), ncol = length(p_n_dp))
for (s in 1:length(p_n_mad)) {
  cev_dp_ps <- cev_ps_dp_sim[[s]]
  for (i in 1:NROW(cev_dp_ps)) {
    ci <- robustBLME::hpd(cev_dp_ps[i,])
    if (mu_true[i] <= ci[2] & mu_true[i] >= ci[1]) uq_dp_80[i,s] <- 1
  }
}
cv_dp_80 <- rowMeans(matrix(unlist(uq_dp_80), ncol = length(p_n_dp), nrow = NROW(covs)))

mean(mse_glm_80); sd(mse_glm_80)
mean(mse_bart_80); sd(mse_bart_80)
mean(mse_rf_80); sd(mse_rf_80)
mean(mse_dp_80); sd(mse_dp_80)
mean(mse_mad_80_1); sd(mse_mad_80_1)
mean(mse_mad_80_23); sd(mse_mad_80_23)
mean(mse_mad_80_dpm); sd(mse_mad_80_dpm)
mean(mse_mad_80_a34); sd(mse_mad_80_a34)
boxplot(cv_glm_80, cv_bart_80, cv_rf_80, cv_dp_80, cv_mad_80_1, cv_mad_80_23, cv_mad_80_dpm, cv_mad_80_a34); abline(h = 0.95, col = 2)



cv_df <- data.frame(cv = c(cv_glm_40, cv_bart_40, cv_rf_40, cv_dp_40, cv_mad_40_1, cv_mad_40_23, cv_mad_40_FHW, cv_mad_40_a34,
                           cv_glm_80, cv_bart_80, cv_rf_80, cv_dp_80, cv_mad_80_1, cv_mad_80_23, cv_mad_80_FHW, cv_mad_80_a34),
                    model = rep(c("1", "2", "3", "4", "5", "6", "7", "8"), 2, each = 1024),
                    n = rep(c("n = 40", "n = 80"), each = 8*1024))
# save(cv_df, file = file("coverage_regr.RData", blocking = TRUE))
