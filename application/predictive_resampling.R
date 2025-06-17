
rm(list = ls())

library(Rcpp)
library(RcppArmadillo)
library(SparseArray)

data_2009 <- read.csv("corvus_2009.csv")
y <- data_2009[,1:4]
temp <- as.numeric(data_2009$Temperature)
hab <- as.numeric(data_2009$Habitat)
sites_2009 <- data_2009$Site

data_2010 <- read.csv("corvus_2010.csv")
y_test <- data_2010[,1:4]
temp_test <- as.numeric(data_2010$Temperature)
hab_test <- as.numeric(data_2010$Habitat)
sites_2010 <- data_2010$Site


# MAD Sequence
Rcpp::sourceCpp("Functions_Corvids.cpp")
source("Functions_Corvids.R")
D <- 4
supp <- 0:50
supp_temp <- 0:2
supp_hab <- 1:5
p_0 <- SVT_SparseArray(dim = c(rep(length(supp), D), length(supp_temp), length(supp_hab)))
p_0_az <- 1 / (length(supp)^D * length(supp_temp) * length(supp_hab))
alpha <- 1
lambda <- 3/4
N_star <- 500

opt <- nlminb(start = c(4, 2, 1.5, 2, 0.5, 0.4),
              lower = c(1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10),
              upper = c(Inf, Inf, Inf, Inf, 1, 1),
              objective = function(h) -preq_4_FB_perm(as.matrix(y), temp, hab, p_0, p_0_az,
                                                      alpha, lambda, N_star, h[1:4], h[5:6], supp, supp_temp, supp_hab,
                                                      10, 10, 1),
              control = list(rel.tol = 1e-3))
save(opt, file = file("corv_opt_a.RData"))
sigma_opt <- opt$par[1:4]
delta_opt <- opt$par[5:6]
pp <- pred_4_FB_perm(as.matrix(y), temp, hab, p_0, p_0_az, alpha, lambda, N_star, sigma_opt, 
                     delta_opt, supp, supp_temp, supp_hab, 10, 10, 12345)

# Predictive resampling
ps <- pred_res_4_FB(pp, N = 1000, B = 1000, cores = 20, seed = 1)
save(ps, file = file("ps_a.RData", blocking = TRUE))

