rm(list = ls())

library(parallel)
library(rmp)
library(dplyr)
library(tidyr)
library(ggplot2)
library(latex2exp)

Rcpp::sourceCpp("illustration_functions.cpp")
emp_hpd_lower <- function(x, conf_level = 0.95) robustBLME::hpd(x, conf_level)[1]
emp_hpd_upper <- function(x, conf_level = 0.95) robustBLME::hpd(x, conf_level)[2]
pred_mad <- function(y, p_0, alpha, lambda, N_star, sd_rg, support, nPerm, seed = 1, cores = 5) {
  idx_perm <- y_poss_perm <- matrix(nrow = length(y), ncol = nPerm) 
  set.seed(seed)
  for (p in 1:nPerm) {
    idx_perm[,p] <- sample(1:length(y), size = length(y), replace = FALSE)
    y_poss_perm[,p] <- y[idx_perm[,p]]
  }
  
  out <- mclapply(1:nPerm, mclapply_function_cpp, mc.cores = cores, y = y_poss_perm, 
                  p_0 = p_0, sd_rg = sd_rg, support = support, lambda = lambda, 
                  N_star = N_star, alpha = alpha)
  out <- rowMeans(matrix(unlist(out), nrow = length(support), ncol = nPerm))
  out
}
preq_ll <- function(sd_rg, y, p_0, alpha, lambda, N_star, support, nPerm, seed = 1, cores = 5) {
  idx_perm <- y_poss_perm <- matrix(nrow = length(y), ncol = nPerm) 
  set.seed(seed)
  for (p in 1:nPerm) {
    idx_perm[,p] <- sample(1:length(y), size = length(y), replace = FALSE)
    y_poss_perm[,p] <- y[idx_perm[,p]]
  }
  
  out <- mclapply(1:nPerm, mclapply_function_pq_cpp, mc.cores = cores, y = y_poss_perm,
                  p_0 = p_0, sd_rg = sd_rg, support = support, lambda = lambda, 
                  N_star = N_star, alpha = alpha)
  out <- mean(unlist(out))
  out
}
pred_res <- function(y, p_n, alpha, lambda, N_star, support, N, B, sd_rg, seed = 1, cores = 5) {
  require(parallel)
  m <- length(y)
  
  set.seed(seed)
  out <- mclapply(1:B, mclapply_function_ps_cpp, mc.cores = cores,
                  p_n = p_n, m = m, N = N, sd_rg = sd_rg, lambda = lambda,
                  N_star = N_star, alpha = alpha, support = support)
  out <- matrix(unlist(out), nrow = length(support), ncol = B)
  out
}
mclapply_function_ps_conv_val <- function(b, p_n, m, N, sd_rg, lambda, N_star, alpha, support, val, seed) {
  set.seed(seed*b)
  mclapply_function_ps_conv_val_cpp(val, b, p_n, m, N, sd_rg, lambda, N_star, alpha, support)
}
pred_res_conv_val <- function(val, y, p_n, alpha, lambda, N_star, support, N, B, sd_rg, seed = 1, cores = 5) {
  require(parallel)
  m <- length(y)
  
  out <- mclapply(1:B, mclapply_function_ps_conv_val, mc.cores = cores,
                  p_n = p_n, m = m, N = N, sd_rg = sd_rg, lambda = lambda,
                  N_star = N_star, alpha = alpha, support = support, val = val, seed = seed)
  out <- matrix(unlist(out), nrow = N - m, ncol = B)
  out
}

mclapply_function_ps_conv <- function(b, p_n, m, N, sd_rg, lambda, N_star, alpha, support, seed) {
  set.seed(seed*b)
  mclapply_ps_conv_cpp(p_n, m, N, sd_rg, lambda, N_star, alpha, support)
}
pred_res_conv <- function(y, p_n, alpha, lambda, N_star, support, N, B, sd_rg, seed = 1, cores = 5) {
  require(parallel)
  m <- length(y)
  
  out <- mclapply(1:B, mclapply_function_ps_conv, mc.cores = cores,
                  p_n = p_n, m = m, N = N, sd_rg = sd_rg, lambda = lambda,
                  N_star = N_star, alpha = alpha, support = support, seed = seed)
  out <- matrix(unlist(out), nrow = N - m, ncol = B)
  out
}

pred_mad_cons <- function(y, sample_size, p_0, alpha, lambda, N_star, sd_rg, support, nPerm, seed = 1, cores = 5) {
  idx_perm <- y_poss_perm <- matrix(nrow = length(y), ncol = nPerm) 
  set.seed(seed)
  for (p in 1:nPerm) {
    idx_perm[,p] <- sample(1:length(y), size = length(y), replace = FALSE)
    y_poss_perm[,p] <- y[idx_perm[,p]]
  }
  
  out <- mclapply(1:nPerm, mclapply_function_cons_cpp, mc.cores = cores, y = y_poss_perm, 
                  p_0 = p_0, sd_rg = sd_rg, support = support, lambda = lambda, 
                  N_star = N_star, alpha = alpha, sample_size = sample_size)
  out <- rowMeans(matrix(unlist(out), nrow = length(support), ncol = nPerm))
  out
}

# Data simulation
nn <- 500
set.seed(1)
y <- sample(c(1, 18, 40, 62), size = nn, replace = TRUE, prob = c(3,2,3,2))
for (i in 1:nn) {
  y[i] <- rpois(1, lambda = y[i])
}
supp <- 0:100
p_true <- 0.3*dpois(supp, lambda = 1) + 0.2*dpois(supp, lambda = 18) + 0.3*dpois(supp, lambda = 40) + 0.2*dpois(supp, lambda = 62)
N <- 1e4
B <- 1000
N_star <- 500
p_0 <- rep(1/length(supp), length(supp))
alpha <- 1
PP <- 10


# Canale & Dunson
set.seed(1234)
idx_rgdpm <- 3
system.time(rgdpm <- rmp::rmg(ydis = y, nrep = (idx_rgdpm*2)*B, nb = idx_rgdpm*B, alpha = 1))
pm_rgdpm <- ci_low_rgdpm <- ci_upp_rgdpm <- rep(0, length(supp))
for (i in 1:length(rgdpm$pmf$domain)) {
  if (rgdpm$pmf$domain[i] %in% supp) {
    idx_supp <- which(supp == rgdpm$pmf$domain[i])
    pm_rgdpm[idx_supp] <- rgdpm$pmf$post.pmf[i]
    ci_low_rgdpm[idx_supp] <- rgdpm$pmf$lower.95[i]
    ci_upp_rgdpm[idx_supp] <- rgdpm$pmf$upper.95[i]
  }  
}
plot(rgdpm$mcmc.chains$pmf[,50], type = "l")
plot(rgdpm$mcmc.chains$pmf[,6], type = "l")


# MAD-ada
LL <- 0.75
system.time(opt_2 <- nlminb(start = 1, lower = 1e-10, upper = Inf,
                            objective = function(x) -preq_ll(x, y, p_0, alpha, LL, N_star, supp, PP, cores = 8)))
sd_2 <- opt_2$par
system.time(pred_2 <- pred_mad(y, p_0, alpha, LL, N_star, sd_2, supp, PP, cores = 5))
system.time(ps_2 <- pred_res(y, pred_2, alpha, LL, N_star, supp, N, B, sd_2, cores = 8))
ci_low_2 <- apply(ps_2, 1, emp_hpd_lower); ci_upp_2 <- apply(ps_2, 1, emp_hpd_upper)


# MAD-1
LL <- 1
system.time(opt_1 <- nlminb(start = 1, lower = 1e-10, upper = Inf,
                            objective = function(x) -preq_ll(x, y, p_0, alpha, LL, 1, supp, PP)))
opt_1
sd_1 <- opt_1$par
system.time(pred_1 <- pred_mad(y, p_0, alpha, LL, 1, sd_1, supp, PP))
system.time(ps_1 <- pred_res(y, pred_1, alpha, LL, 1, supp, N, B, sd_1, seed = 1))
ci_low_1 <- apply(ps_1, 1, emp_hpd_lower); ci_upp_1 <- apply(ps_1, 1, emp_hpd_upper); pm_1 <- apply(ps_1, 1, mean)


# DP
LL <- 1
sd_dp <- 1e-10
pred_dp <- pred_mad(y, p_0, alpha, LL, N_star, sd_dp, supp, PP)
ps_dp <- pred_res(y, pred_dp, alpha, LL, N_star, supp, N, B, sd_dp)
ci_low_dp <- apply(ps_dp, 1, emp_hpd_lower); ci_upp_dp <- apply(ps_dp, 1, emp_hpd_upper)



df <- data.frame(supp = rep(supp, 4),
                 pm = c(pm_1, pred_2, pred_dp, pm_rgdpm),
                 ci_low = c(ci_low_1, ci_low_2, ci_low_dp, ci_low_rgdpm),
                 ci_upp = c(ci_upp_1, ci_upp_2, ci_upp_dp, ci_upp_rgdpm),
                 setting = rep(1:4, each = length(supp)),
                 p_true = rep(p_true, 4))
df$setting <- as.factor(df$setting)
levels(df$setting) <- c("MAD-1", "MAD-ada", "DP", "RG-DPM")
ggplot(df) + 
  facet_wrap(~ setting) + 
  geom_segment(aes(x = supp, xend = supp, y = ci_low, yend = ci_upp, col = "z_ci"), alpha = 0.15, size = 2) +
  geom_segment(aes(x = supp, xend = supp, y = 0, yend = pm, col = "pm"), alpha = 1, size = .4) + 
  geom_point(aes(x = supp, y = p_true), col = "forestgreen", size = .5) +
  scale_color_manual(name = "", values = c("pm" = "red", "z_ci" = "blue3"), 
                     labels = c("pm" = TeX("$P_n$"), "z_ci" = "Posterior 0.95 CI")) +
  ylim(c(0,0.185)) + xlab("Support") + ylab("Probability") + xlim(c(0,85)) +
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



# CONVERGENCE of MARTINGALE for P(Y<26)
LL <- 0.75
N_star <- 500
y_val <- 26
ps_val <- pred_res_conv_val(y_val, y, pred_2, alpha, LL, N_star, supp, 
                           N = 20000, B = 100, sd_2, seed = 123)
dim(ps_val)

df_conv_val <- as.data.frame(ps_val) %>%
  mutate(row = 1:n()) %>%
  pivot_longer(cols = -row, names_to = "column", values_to = "value")
df_conv_val$title <- "Convergence for MAD-ada"
df_conv_val$row <- df_conv_val$row + nn

pl_conv_val <- ggplot(df_conv_val) +
  facet_wrap(~ title) +
  geom_line(aes(x = row, y = value, group = column), alpha = 0.2) +
  ylim(0,1) + xlim(c(nn, nn+NROW(ps_val))) +
  labs(x = "Number of forward samples", y = "Probability") +
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
pl_conv_val



# CONVERGENCE of P_n in L_1
LL <- 0.75
N_star <- 500
cc <- pred_res_conv(y, pred_2, alpha, LL, N_star, supp, N = 20000, B = 100, 
                    sd_rg = sd_2, seed = 123)
plot(cc[,1], type = "l"); abline(v = N, lty = 2)
df_conv <- as.data.frame(cc) %>%
  mutate(row = 1:n()) %>%
  pivot_longer(cols = -row, names_to = "column", values_to = "value")
df_conv$title <- "Convergence for MAD-ada"
df_conv$row <- df_conv$row + nn
pl_conv <- ggplot(df_conv) +
  facet_wrap(~ title) +
  geom_line(aes(x = row, y = value, group = column), alpha = 0.2) +
  ylim(0,1) + xlim(c(nn, nn+NROW(cc))) +
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
pl_conv
# ggsave(pl_conv, filename = "pl_ill_conv.pdf", width = 9, height = 6)



# MARTINGALE POSTERIOR for P(Y<26)
y_val <- 26
ps_mar_1 <- apply(ps_1, 2, function(x) sum(x[1:y_val]))
ps_mar_2 <- apply(ps_2, 2, function(x) sum(x[1:y_val]))
ps_mar_dp <- apply(ps_dp, 2, function(x) sum(x[1:y_val]))
df <- data.frame(p1 = ps_mar_1, p2 = ps_mar_2, dp = ps_mar_dp, title = "Martingale posterior for P{Y < 26}")
pl_mart <- ggplot(df) + 
  facet_wrap(~ title) +
  stat_density(aes(x = dp, col = "dp"), linetype = 2, geom = "line", position = "identity", adjust = 2.5) +
  stat_density(aes(x = p1, col = "1"), linetype = 1, geom = "line", position = "identity", adjust = 2.5) +
  stat_density(aes(x = p2, col = "2"), linetype = 1, geom = "line", position = "identity", adjust = 2.5) +
  scale_color_manual(name = "lambda",
                     values = c("dp" = "grey30",
                                "2" = "darkorange",
                                "1" = "blue2"),
                     labels = c("dp" = "DP", 
                                "2" = "MAD-ada",
                                "1" = "MAD-1"
                     )) +
  xlim(c(0,1)) + xlab("Probability") +
  ylab("Density") +
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
pl_mart

