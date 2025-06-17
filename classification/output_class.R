
rm(list = ls())

library(mvtnorm)
library(truncnorm)
library(fastmatrix)
library(ggplot2)
library(coda)

load("SIM_Class_150_mad_glm_bart_rf.RData")
auc_glm_150 <- auc_glm; auc_bart_150 <- auc_bart; auc_rf_150 <- auc_rf
cv_rf_150 <- rowMeans(matrix(unlist(uq_rf), ncol = nSim, nrow = NROW(covs)))
cv_bart_150 <- rowMeans(matrix(unlist(uq_bart), ncol = nSim, nrow = NROW(covs)))
cv_glm_150 <- rowMeans(matrix(unlist(uq_glm), ncol = nSim, nrow = NROW(covs)))

load("SIM_Class_150_mad_a34.RData")
auc_mad_a34_150 <- auc_mad
uq_mad_a34_150 <- matrix(0, nrow = NROW(covs), ncol = nSim)
for (s in 1:nSim) {
  cev_mad_ps <- uq_mad_pc_ps[[s]]
  for (i in 1:NROW(cev_mad_ps)) {
    ci <- robustBLME::hpd(cev_mad_ps[i,])
    if (p_true[i] <= ci[2] & p_true[i] >= ci[1]) uq_mad_a34_150[i,s] <- 1
  }
}
cv_mad_a34_150 <- rowMeans(matrix(unlist(uq_mad_a34_150), ncol = nSim, nrow = NROW(covs)))

load("SIM_Class_150_mad_1.RData")
auc_mad_1_150 <- auc_mad
uq_mad_1_150 <- matrix(0, nrow = NROW(covs), ncol = nSim)
for (s in 1:nSim) {
  cev_mad_ps <- uq_mad_pc_ps[[s]]
  for (i in 1:NROW(cev_mad_ps)) {
    ci <- robustBLME::hpd(cev_mad_ps[i,])
    if (p_true[i] <= ci[2] & p_true[i] >= ci[1]) uq_mad_1_150[i,s] <- 1
  }
}
cv_mad_1_150 <- rowMeans(matrix(unlist(uq_mad_1_150), ncol = nSim, nrow = NROW(covs)))

load("SIM_Class_150_mad_23.RData")
auc_mad_23_150 <- auc_mad
uq_mad_23_150 <- matrix(0, nrow = NROW(covs), ncol = nSim)
for (s in 1:nSim) {
  cev_mad_ps <- uq_mad_pc_ps[[s]]
  for (i in 1:NROW(cev_mad_ps)) {
    ci <- robustBLME::hpd(cev_mad_ps[i,])
    if (p_true[i] <= ci[2] & p_true[i] >= ci[1]) uq_mad_23_150[i,s] <- 1
  }
}
cv_mad_23_150 <- rowMeans(matrix(unlist(uq_mad_23_150), ncol = nSim, nrow = NROW(covs)))

load("SIM_Class_150_mad_dpm.RData")
auc_mad_dpm_150 <- auc_mad
uq_mad_dpm_150 <- matrix(0, nrow = NROW(covs), ncol = nSim)
for (s in 1:nSim) {
  cev_mad_ps <- uq_mad_pc_ps[[s]]
  for (i in 1:NROW(cev_mad_ps)) {
    ci <- robustBLME::hpd(cev_mad_ps[i,])
    if (p_true[i] <= ci[2] & p_true[i] >= ci[1]) uq_mad_dpm_150[i,s] <- 1
  }
}
cv_mad_dpm_150 <- rowMeans(matrix(unlist(uq_mad_dpm_150), ncol = nSim, nrow = NROW(covs)))

load("SIM_Class_150_DP.RData")
auc_dp_150 <- auc_mad
uq_dp_150 <- matrix(0, nrow = NROW(covs), ncol = nSim)
for (s in 1:nSim) {
  cev_mad_ps <- uq_mad_pc_ps[[s]]
  for (i in 1:NROW(cev_mad_ps)) {
    ci <- robustBLME::hpd(cev_mad_ps[i,])
    if (p_true[i] <= ci[2] & p_true[i] >= ci[1]) uq_dp_150[i,s] <- 1
  }
}
cv_dp_150 <- rowMeans(matrix(unlist(uq_dp_150), ncol = nSim, nrow = NROW(covs)))

mean(auc_glm_150); sd(auc_glm_150)
mean(auc_bart_150); sd(auc_bart_150)
mean(auc_rf_150); sd(auc_rf_150)
mean(auc_dp_150); sd(auc_dp_150)
mean(auc_mad_1_150); sd(auc_mad_1_150)
mean(auc_mad_23_150); sd(auc_mad_23_150)
mean(auc_mad_dpm_150); sd(auc_mad_dpm_150)
mean(auc_mad_a34_150); sd(auc_mad_a34_150)
boxplot(cv_glm_150, cv_bart_150, cv_rf_150, cv_dp_150, cv_mad_1_150, cv_mad_23_150, cv_mad_dpm_150, cv_mad_a34_150, ylim = c(0,1)); abline(h = 0.95, col = 2)






load("SIM_Class_300_mad_glm_bart_rf.RData")
auc_glm_300 <- auc_glm; auc_bart_300 <- auc_bart; auc_rf_300 <- auc_rf
cv_rf_300 <- rowMeans(matrix(unlist(uq_rf), ncol = nSim, nrow = NROW(covs)))
cv_bart_300 <- rowMeans(matrix(unlist(uq_bart), ncol = nSim, nrow = NROW(covs)))
cv_glm_300 <- rowMeans(matrix(unlist(uq_glm), ncol = nSim, nrow = NROW(covs)))

load("SIM_Class_300_mad_a34.RData")
auc_mad_a34_300 <- auc_mad
uq_mad_a34_300 <- matrix(0, nrow = NROW(covs), ncol = nSim)
for (s in 1:nSim) {
  cev_mad_ps <- uq_mad_pc_ps[[s]]
  for (i in 1:NROW(cev_mad_ps)) {
    ci <- robustBLME::hpd(cev_mad_ps[i,])
    if (p_true[i] <= ci[2] & p_true[i] >= ci[1]) uq_mad_a34_300[i,s] <- 1
  }
}
cv_mad_a34_300 <- rowMeans(matrix(unlist(uq_mad_a34_300), ncol = nSim, nrow = NROW(covs)))

load("SIM_Class_300_mad_1.RData")
auc_mad_1_300 <- auc_mad
uq_mad_1_300 <- matrix(0, nrow = NROW(covs), ncol = nSim)
for (s in 1:nSim) {
  cev_mad_ps <- uq_mad_pc_ps[[s]]
  for (i in 1:NROW(cev_mad_ps)) {
    ci <- robustBLME::hpd(cev_mad_ps[i,])
    if (p_true[i] <= ci[2] & p_true[i] >= ci[1]) uq_mad_1_300[i,s] <- 1
  }
}
cv_mad_1_300 <- rowMeans(matrix(unlist(uq_mad_1_300), ncol = nSim, nrow = NROW(covs)))

load("SIM_Class_300_mad_23.RData")
auc_mad_23_300 <- auc_mad
uq_mad_23_300 <- matrix(0, nrow = NROW(covs), ncol = nSim)
for (s in 1:nSim) {
  cev_mad_ps <- uq_mad_pc_ps[[s]]
  for (i in 1:NROW(cev_mad_ps)) {
    ci <- robustBLME::hpd(cev_mad_ps[i,])
    if (p_true[i] <= ci[2] & p_true[i] >= ci[1]) uq_mad_23_300[i,s] <- 1
  }
}
cv_mad_23_300 <- rowMeans(matrix(unlist(uq_mad_23_300), ncol = nSim, nrow = NROW(covs)))

load("SIM_Class_300_mad_dpm.RData")
auc_mad_dpm_300 <- auc_mad
uq_mad_dpm_300 <- matrix(0, nrow = NROW(covs), ncol = nSim)
for (s in 1:nSim) {
  cev_mad_ps <- uq_mad_pc_ps[[s]]
  for (i in 1:NROW(cev_mad_ps)) {
    ci <- robustBLME::hpd(cev_mad_ps[i,])
    if (p_true[i] <= ci[2] & p_true[i] >= ci[1]) uq_mad_dpm_300[i,s] <- 1
  }
}
cv_mad_dpm_300 <- rowMeans(matrix(unlist(uq_mad_dpm_300), ncol = nSim, nrow = NROW(covs)))

load("SIM_Class_300_DP.RData")
auc_dp_300 <- auc_mad
uq_dp_300 <- matrix(0, nrow = NROW(covs), ncol = nSim)
for (s in 1:nSim) {
  cev_mad_ps <- uq_mad_pc_ps[[s]]
  for (i in 1:NROW(cev_mad_ps)) {
    ci <- robustBLME::hpd(cev_mad_ps[i,])
    if (p_true[i] <= ci[2] & p_true[i] >= ci[1]) uq_dp_300[i,s] <- 1
  }
}
cv_dp_300 <- rowMeans(matrix(unlist(uq_dp_300), ncol = nSim, nrow = NROW(covs)))

mean(auc_glm_300); sd(auc_glm_300)
mean(auc_bart_300); sd(auc_bart_300)
mean(auc_rf_300); sd(auc_rf_300)
mean(auc_dp_300); sd(auc_dp_300)
mean(auc_mad_1_300); sd(auc_mad_1_300)
mean(auc_mad_23_300); sd(auc_mad_23_300)
mean(auc_mad_dpm_300); sd(auc_mad_dpm_300)
mean(auc_mad_a34_300); sd(auc_mad_a34_300)
boxplot(cv_glm_300, cv_bart_300, cv_rf_300, cv_dp_300, cv_mad_1_300, cv_mad_23_300, cv_mad_dpm_300, cv_mad_a34_300, ylim = c(0,1)); abline(h = 0.95, col = 2)




cv_df <- data.frame(cv = c(cv_glm_150, cv_bart_150, cv_rf_150, cv_dp_150, cv_mad_1_150, cv_mad_23_150, cv_mad_FHW_150, cv_mad_a34_150,
                           cv_glm_300, cv_bart_300, cv_rf_300, cv_dp_300, cv_mad_1_300, cv_mad_23_300, cv_mad_FHW_300, cv_mad_a34_300),
                    model = rep(c("1", "2", "3", "4", "5", "6", "7", "8"), 2, each = 1024),
                    n = rep(c("n = 150", "n = 300"), each = 8*1024))
# save(cv_df, file = file("coverage_class.RData", blocking = TRUE))
