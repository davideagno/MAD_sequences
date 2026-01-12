
rm(list = ls())
setwd("/hpc/group/dunsonlab/da252/review")

library(Rcpp)
library(RcppArmadillo)
library(SparseArray)

Rcpp::sourceCpp("appl_function_2.cpp")
source("Functions_appl.R")
D <- 4
supp <- 0:50
supp_temp <- 0:2
supp_hab <- 1:5
p_0 <- SVT_SparseArray(dim = c(rep(length(supp), D), length(supp_temp), length(supp_hab)))
p_0_az <- 1 / (length(supp)^D * length(supp_temp) * length(supp_hab))
alpha <- 1
lambda <- 3/4
N_star <- 500

# sigma_opt <- c(3.7991043, 1.5306404, 1.0818673, 1.2445321)
# delta_opt <- c(0.4439974, 0.4421497)
# load("corv_opt_a.RData")
load("corv_pp.RData")

# Predictive resampling
ps_pr <- pred_res_4_FB(pp, N = 1000, B = 50, cores = 10, seed = 1)
save(ps_pr, file = file("ps_prova.RData", blocking = TRUE))

ps_1 <- pred_res_4_FB(pp, N = 5000, B = 100, cores = 10, seed = 1)
save(ps_1, file = file("ps_1a_long.RData", blocking = TRUE))
ps_2 <- pred_res_4_FB(pp, N = 5000, B = 100, cores = 10, seed = 101)
save(ps_2, file = file("ps_2a_long.RData", blocking = TRUE))
ps_3 <- pred_res_4_FB(pp, N = 5000, B = 100, cores = 10, seed = 301)
save(ps_3, file = file("ps_3a_long.RData", blocking = TRUE))
ps_4 <- pred_res_4_FB(pp, N = 5000, B = 100, cores = 10, seed = 401)
save(ps_4, file = file("ps_4a_long.RData", blocking = TRUE))
ps_5 <- pred_res_4_FB(pp, N = 5000, B = 100, cores = 10, seed = 501)
save(ps_5, file = file("ps_5a_long.RData", blocking = TRUE))
ps_6 <- pred_res_4_FB(pp, N = 5000, B = 100, cores = 10, seed = 601)
save(ps_6, file = file("ps_6a_long.RData", blocking = TRUE))
ps_7 <- pred_res_4_FB(pp, N = 5000, B = 100, cores = 10, seed = 701)
save(ps_7, file = file("ps_7a_long.RData", blocking = TRUE))
ps_8 <- pred_res_4_FB(pp, N = 5000, B = 100, cores = 10, seed = 801)
save(ps_8, file = file("ps_8a_long.RData", blocking = TRUE))
ps_9 <- pred_res_4_FB(pp, N = 5000, B = 100, cores = 10, seed = 901)
save(ps_9, file = file("ps_9a_long.RData", blocking = TRUE))
ps_10 <- pred_res_4_FB(pp, N = 5000, B = 100, cores = 10, seed = 1001)
save(ps_10, file = file("ps_10a_long.RData", blocking = TRUE))






