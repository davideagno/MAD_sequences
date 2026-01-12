
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
Rcpp::sourceCpp("Functions_Corvids_2.cpp")
D <- 4
supp <- 0:50
supp_temp <- 0:2
supp_hab <- 1:5
p_0 <- SVT_SparseArray(dim = c(rep(length(supp), D), length(supp_temp), length(supp_hab)))
p_0_az <- 1 / (length(supp)^D * length(supp_temp) * length(supp_hab))
alpha <- 1
lambda <- 3/4
N_star <- 500

load("corv_opt_a.RData")
sigma_opt <- opt$par[1:4]
delta_opt <- opt$par[5:6]
load("corv_pp.RData")

p_n <- pp$p_n
nz_log <- as.numeric(is_nonzero(p_n))
idx_nz <- c(1:length(nz_log))[which(nz_log == 1)]
p_n_vec <- p_n[idx_nz]
p_az <- pp$p_az

D <- 4
require(crch)
left <- min(supp) - 0.5
interval <- rep(0, D)
for (k in 1:D) {
  extremes <- crch::qtnorm(c(0.025, 0.975), mean = 50, sd = sigma_opt[k], left = left)
  interval[k] <- max(0, (extremes[2] - extremes[1]) / 2)
}

conv_appl <- matrix(nrow = 10000, ncol = 10)
for (i in 1:10) {
  set.seed(i)
  print(i)
  conv_appl[,i] <- pred_res_4_FB_one_conv(p_n_vec, p_az, idx_nz, NROW(y), 10000,
                                          alpha, 0.75, N_star, sigma_opt, delta_opt,
                                          supp, supp_temp, supp_hab, interval)
}


library(ggplot2)
library(latex2exp)
df_conv <- data.frame(conv = c(conv_appl),
                      fs = rep(1:NROW(conv_appl), NCOL(conv_appl)),
                      group = rep(1:NCOL(conv_appl), each = NROW(conv_appl)),
                      title = "Convergence")
ggplot(df_conv) +
  facet_wrap(~ title) +
  geom_line(aes(x = fs, y = conv, group = group), alpha = 0.7) +
  labs(x = "Number of forward samples", y = TeX("$L_1$ distance")) +
  theme_light() +
  theme(legend.text.align = 0,
        strip.text = element_text(size = 12, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top")

