
rm(list = ls())

library(mvtnorm)
library(truncnorm)
library(fastmatrix)
library(TeachingDemos)
library(ggplot2)
library(coda)
library(latex2exp)

library(predtools)
library(magrittr)
library(dplyr)
library(plotROC)
library(pROC)

library(parallel)

emp_hpd_lower <- function(x, conf_level = 0.95) robustBLME::hpd(x, conf_level)[1]
emp_hpd_upper <- function(x, conf_level = 0.95) robustBLME::hpd(x, conf_level)[2]
rho <- function(n, lambda, N_star) {
  lambda + (1 - lambda)*exp(-(1/N_star)*n)
}
drg <- function(x, mean, sd) {
  D <- length(c(x))
  out <- rep(NA, D)
  for (i in 1:D) {
    out[i] <- pnorm(x[i]+0.5, mean = mean, sd = sd) - pnorm(x[i]-0.5, mean = mean, sd = sd)
  }
  out <- out / sum(out)
  out[which(out <= 1e-250)] <- 1e-250
  out
}
mhk_rg <- function(p_n, y, sd_rg, support) {
  out <- rep(NA, length(p_n))
  idx_y <- which(support == y)
  prob_acc <- pmin(p_n / p_n[idx_y], 1)
  k0 <- drg(support, mean = y, sd = sd_rg)
  out <- prob_acc * k0
  rej_prob <- 1 - sum(out)
  out[idx_y] <- out[idx_y] + rej_prob
  out
}

pred_rg_rho <- function(y, p_0, alpha, lambda, N_star, sd_rg, support, nPerm, seed = 1, cores = 5) {
  idx_perm <- matrix(nrow = length(y), ncol = nPerm) 
  set.seed(seed)
  for (p in 1:nPerm) {
    idx_perm[,p] <- sample(1:length(y), size = length(y), replace = FALSE)
  }
  
  mclapply_function <- function(i) {
    y_perm <- y[idx_perm[,i]]
    p_n <- p_0
    for (n in 1:length(y_perm)) {
      ker <- mhk_rg(p_n, y_perm[n], sd_rg, support)
      rho_n <- rho(n, lambda, N_star)
      w_n <- 1 / ((alpha + n)^rho_n)
      p_n <- (1 - w_n)*p_n + w_n*ker
    }
    p_n
  }
  
  out <- mclapply(1:nPerm, mclapply_function, mc.cores = cores)
  out <- rowMeans(matrix(unlist(out), nrow = length(support), ncol = nPerm))
  out
}
preq_logl_rg_rho <- function(sd_rg, y, p_0, alpha, lambda, N_star, support, nPerm, seed = 1, cores = 5) {
  idx_perm <- matrix(nrow = length(y), ncol = nPerm) 
  set.seed(seed)
  for (p in 1:nPerm) {
    idx_perm[,p] <- sample(1:length(y), size = length(y), replace = FALSE)
  }
  
  mclapply_function <- function(i) {
    pll <- 0
    y_perm <- y[idx_perm[,i]]
    p_n <- p_0
    for (n in 1:length(y_perm)) {
      y_n <- y_perm[n]
      pll <- pll + log(p_n[which(support == y_n)])
      ker <- mhk_rg(p_n, y_n, sd_rg, support)
      rho_n <- rho(n, lambda, N_star)
      w_n <- 1 / ((alpha + n)^rho_n)
      p_n <- (1 - w_n)*p_n + w_n*ker
    }
    pll
  }
  
  out <- mclapply(1:nPerm, mclapply_function, mc.cores = cores)
  mean(unlist(out))
}
pred_res_rg_rho <- function(y, p_n, alpha, lambda, N_star, support, N, B, sd_rg, seed = 1, cores = 5) {
  require(parallel)
  m <- length(y)
  
  mclapply_function <- function(i) {
    prob <- p_n
    for (n in (m+1):N) {
      y_past <- sample(support, 1, prob = prob)
      ker <- mhk_rg(p_n, y_past, sd_rg, support)
      rho_n <- rho(n, lambda, N_star)
      w_n <- 1 / ((alpha + n)^rho_n)
      prob <- (1 - w_n)*prob + w_n*ker
    }
    prob
  }
  
  out <- mclapply(1:B, mclapply_function, mc.cores = cores)
  out <- matrix(unlist(out), nrow = length(support), ncol = B)
  out
}

pred_rg_fixed <- function(y, p_0, alpha, lambda, sd_rg, support, nPerm, seed = 1, cores = 5) {
  idx_perm <- matrix(nrow = length(y), ncol = nPerm) 
  set.seed(seed)
  for (p in 1:nPerm) {
    idx_perm[,p] <- sample(1:length(y), size = length(y), replace = FALSE)
  }
  
  mclapply_function <- function(i) {
    y_perm <- y[idx_perm[,i]]
    p_n <- p_0
    for (n in 1:length(y_perm)) {
      ker <- mhk_rg(p_n, y_perm[n], sd_rg, support)
      w_n <- 1 / ((alpha + n)^lambda)
      p_n <- (1 - w_n)*p_n + w_n*ker
    }
    p_n
  }
  
  out <- mclapply(1:nPerm, mclapply_function, mc.cores = cores)
  out <- rowMeans(matrix(unlist(out), nrow = length(support), ncol = nPerm))
  out
}
preq_logl_rg_fixed <- function(sd_rg, y, p_0, alpha, lambda, support, nPerm, seed = 1, cores = 5) {
  idx_perm <- matrix(nrow = length(y), ncol = nPerm) 
  set.seed(seed)
  for (p in 1:nPerm) {
    idx_perm[,p] <- sample(1:length(y), size = length(y), replace = FALSE)
  }
  
  mclapply_function <- function(i) {
    pll <- 0
    y_perm <- y[idx_perm[,i]]
    p_n <- p_0
    for (n in 1:length(y_perm)) {
      y_n <- y_perm[n]
      pll <- pll + log(p_n[which(support == y_n)])
      ker <- mhk_rg(p_n, y_n, sd_rg, support)
      w_n <- 1 / ((alpha + n)^lambda)
      p_n <- (1 - w_n)*p_n + w_n*ker
    }
    pll
  }
  
  out <- mclapply(1:nPerm, mclapply_function, mc.cores = cores)
  mean(unlist(out))
}
pred_res_rg_fixed <- function(y, p_n, alpha, support, N, B, sd_rg, lambda, seed = 1, cores = 5) {
  require(parallel)
  m <- length(y)
  
  mclapply_function <- function(i) {
    prob <- p_n
    for (n in (m+1):N) {
      y_past <- sample(support, 1, prob = prob)
      ker <- mhk_rg(p_n, y_past, sd_rg, support)
      w_n <- 1 / ((alpha + n)^lambda)
      prob <- (1 - w_n)*prob + w_n*ker
    }
    prob
  }
  
  out <- mclapply(1:B, mclapply_function, mc.cores = cores)
  out <- matrix(unlist(out), nrow = length(support), ncol = B)
  out
}

# Data simulation
nn <- 50
set.seed(12)
y <- rbinom(nn, size = 1, prob = 0.6)
for (i in 1:nn) {
  if (y[i] == 0) {
    y[i] <- rpois(1, lambda = 25)
  } else {
    y[i] <- rpois(1, lambda = 60)
  }
}
supp <- 0:100
N <- 1000
B <- 500
p_0 <- rep(1/length(supp), length(supp))
alpha <- 1
p_true <- 0.4*dpois(supp, lambda = 25) + 0.6*dpois(supp, lambda = 60)
PP <- 10


# MAD-1
LL <- 1
opt_1 <- nlminb(start = 1, lower = 1e-10, upper = Inf,
                objective = function(x) -preq_logl_rg_fixed(x, y, p_0, alpha, LL, supp, PP))
opt_1
sd_1 <- opt_1$par
p_n_1 <- pred_rg_fixed(y, p_0, alpha, LL, sd_1, supp, PP)
ps_1 <- pred_res_rg_fixed(y, p_n_1, alpha, supp, N, B, sd_1, LL)
ci_lower_1 <- apply(ps_1, 1, emp_hpd_lower)
ci_upper_1 <- apply(ps_1, 1, emp_hpd_upper)


# MAD-ada
LL <- 3/4
NN <- 500
opt_1_34 <- nlminb(start = 1, lower = 1e-10, upper = Inf,
                   objective = function(x) -preq_logl_rg_rho(x, y, p_0, alpha, LL, NN, supp, PP))
opt_1_34
sd_1_34 <- opt_1_34$par
p_n_1_34 <- pred_rg_rho(y, p_0, alpha, LL, NN, sd_1_34, supp, PP)
ps_1_34 <- pred_res_rg_rho(y, p_n_1_34, alpha, LL, NN, supp, N, B, sd_1_34)
ci_lower_1_34 <- apply(ps_1_34, 1, emp_hpd_lower)
ci_upper_1_34 <- apply(ps_1_34, 1, emp_hpd_upper)


# DP
LL <- 1
sd_dp <- 1e-10
p_n_dp <- pred_rg_fixed(y, p_0, alpha, LL, sd_dp, supp, PP)
ps_dp <- pred_res_rg_fixed(y, p_n_dp, alpha, supp, N, B, sd_dp, LL)
ci_lower_dp <- apply(ps_dp, 1, emp_hpd_lower)
ci_upper_dp <- apply(ps_dp, 1, emp_hpd_upper)


# Martingale posterior for Pr{Y <= 25}
pr0_1 <- pr0_1_34 <- pr0_dp <- rep(NA, B)
for (b in 1:B) {
  pr0_1[b] <- sum(ps_1[1:26, b])
  pr0_1_34[b] <- sum(ps_1_34[1:26, b])
  pr0_dp[b] <- sum(ps_dp[1:26, b])
}



df <- data.frame(supp = rep(supp, 3),
                 p_mad = c(p_n_1, p_n_1_34, p_n_dp),
                 ci_lower = c(ci_lower_1, ci_lower_1_34, ci_lower_dp),
                 ci_upper = c(ci_upper_1, ci_upper_1_34, ci_upper_dp),
                 setting = rep(c("MAD-1", "MAD-ada", "DP"), each = length(supp)),
                 p_true = rep(p_true, 3))

df_1 <- data.frame(supp = supp, p_mad = p_n_1, ci_lower = ci_lower_1, ci_upper = ci_upper_1, title = "MAD-1")
pl_1 <- ggplot(df_1) + 
  facet_wrap(~ title) +
  geom_segment(aes(x = supp, xend = supp, y = ci_lower, yend = ci_upper, col = "z_ci"), alpha = 0.15, size = 2) +
  geom_segment(aes(x = supp, xend = supp, y = 0, yend = p_mad, col = "pm"), alpha = 1, size = .4) +
  geom_point(aes(x = supp, y = p_true), col = "forestgreen", size = .5) +
  scale_color_manual(name = "", values = c("pm" = "red", "z_ci" = "blue3"), labels = c("pm" = TeX("$P_n$"),
                                                                                       "z_ci" = "Posterior 0.95 CI")) +
  labs(x = "Support", y = "Probability") +
  ylim(c(0, 0.075)) +
  theme_light() +
  theme(legend.position = "top",
        strip.text = element_text(size = 11, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 11),
        axis.text.y = element_text(size = 5))
pl_1

df_1_34 <- data.frame(supp = supp, p_mad = p_n_1_34, ci_lower = ci_lower_1_34, ci_upper = ci_upper_1_34, title = "MAD-ada")
pl_1_34 <- ggplot(df_1_34) + 
  facet_wrap(~ title) +
  geom_segment(aes(x = supp, xend = supp, y = ci_lower, yend = ci_upper, col = "z_ci"), alpha = 0.15, size = 2) +
  geom_segment(aes(x = supp, xend = supp, y = 0, yend = p_mad, col = "pm"), alpha = 1, size = .4) +
  geom_point(aes(x = supp, y = p_true), col = "forestgreen", size = .5) +
  scale_color_manual(name = "", values = c("pm" = "red", "z_ci" = "blue3"), labels = c("pm" = TeX("$P_n$"),
                                                                                       "z_ci" = "Posterior 0.95 CI")) +
  labs(x = "Support", y = "Probability") +
  ylim(c(0, 0.075)) +
  theme_light() +
  theme(legend.position = "top",
        strip.text = element_text(size = 11, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 11),
        axis.text.y = element_text(size = 5))
pl_1_34

df_dp <- data.frame(supp = supp, p_mad = p_n_dp, ci_lower = ci_lower_dp, ci_upper = ci_upper_dp, title = "DP")
pl_dp <- ggplot(df_dp) + 
  facet_wrap(~ title) +
  geom_segment(aes(x = supp, xend = supp, y = ci_lower, yend = ci_upper, col = "z_ci"), alpha = 0.15, size = 2) +
  geom_segment(aes(x = supp, xend = supp, y = 0, yend = p_mad, col = "pm"), alpha = 1, size = .4) +
  geom_point(aes(x = supp, y = p_true), col = "forestgreen", size = .5) +
  scale_color_manual(name = "", values = c("pm" = "red", "z_ci" = "blue3"), labels = c("pm" = TeX("$P_n$"),
                                                                                       "z_ci" = "Posterior 0.95 CI")) +
  labs(x = "Support", y = "Probability") +
  ylim(c(0, 0.15)) +
  theme_light() +
  theme(legend.position = "top",
        strip.text = element_text(size = 11, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 11),
        axis.text.y = element_text(size = 5))
pl_dp


# Prob{Y <= 25}
pr0_1 <- pr0_1_34 <- pr0_dp <- rep(NA, B)
for (b in 1:B) {
  pr0_1[b] <- sum(ps_1[1:26, b])
  pr0_1_34[b] <- sum(ps_1_34[1:26, b])
  pr0_dp[b] <- sum(ps_dp[1:26, b])
}

df_pr <- data.frame(p1 = pr0_1, p2 = pr0_1_34, p3 = pr0_dp, title = "Martingale posterior for Pr{Y < 26}")
pl_pr <- ggplot(df_pr) + 
  facet_wrap(~ title) +
  stat_density(aes(x = p3, col = "p3"), linetype = 1,
               geom = "line", position = "identity", adjust = 3) +
  stat_density(aes(x = p1, col = "p1"), linetype = 1, 
               geom = "line", position = "identity", adjust = 3) +
  stat_density(aes(x = p2, col = "p2"), linetype = 1,
               geom = "line", position = "identity", adjust = 3) +
  geom_vline(xintercept = sum(p_true[1:26]), linetype = 3) +
  scale_color_manual(name = "Habitat",
                     values = c("p1" = "forestgreen",
                                "p2" = "red2",
                                "p3" = "blue2"),
                     labels = c("p1" = "MAD-1", 
                                "p2" = "MAD-ada",
                                "p3" = "DP"
                     )) +
  guides(colour = guide_legend(nrow = 1)) +
  xlim(c(0,0.6)) + xlab("Probability") +
  ylab("Density") +
  theme_light() +
  theme(legend.text.align = 0,
        strip.text = element_text(size = 11, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 11),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top")
pl_pr

library(ggpubr)
pl_ill <- ggarrange(pl_1, pl_1_34, pl_dp, pl_pr, common.legend = TRUE)
pl_ill
# ggsave(pl_ill, filename = "pl_ill.pdf", width = 12, height = 10)


# Figure for slide
df_mad <- data.frame(supp = supp, p_mad = p_n_1, ci_lower = ci_lower_1, ci_upper = ci_upper_1, title = "MAD")
pl_mad <- ggplot(df_mad) + 
  facet_wrap(~ title) +
  geom_segment(aes(x = supp, xend = supp, y = ci_lower, yend = ci_upper, col = "z_ci"), alpha = 0.15, size = 2) +
  geom_segment(aes(x = supp, xend = supp, y = 0, yend = p_mad, col = "pm"), alpha = 1, size = .4) +
  geom_point(aes(x = supp, y = p_true), col = "forestgreen", size = .5) +
  scale_color_manual(name = "", values = c("pm" = "red", "z_ci" = "blue3"), labels = c("pm" = TeX("$P_n$"),
                                                                                       "z_ci" = "Posterior 0.95 CI")) +
  labs(x = "Support", y = "Probability") +
  ylim(c(0, 0.075)) +
  theme_light() +
  theme(legend.position = "top",
        strip.text = element_text(size = 11, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 11),
        axis.text.y = element_text(size = 5))
pl_mad

df_pr <- data.frame(p1 = pr0_1, p2 = pr0_1_34, p3 = pr0_dp, title = "Martingale posterior for Pr{Y < 26}")
pl_pr_2 <- ggplot(df_pr) + 
  facet_wrap(~ title) +
  stat_density(aes(x = p3, col = "p3"), linetype = 1,
               geom = "line", position = "identity", adjust = 3) +
  stat_density(aes(x = p1, col = "p1"), linetype = 1, 
               geom = "line", position = "identity", adjust = 3) +
  geom_vline(xintercept = sum(p_true[1:26]), linetype = 3) +
  scale_color_manual(name = "Habitat",
                     values = c("p1" = "red2",
                                "p3" = "blue2"),
                     labels = c("p1" = "MAD",
                                "p3" = "DP"
                     )) +
  guides(colour = guide_legend(nrow = 1)) +
  xlim(c(0,0.5)) + xlab("Probability") +
  ylab("Density") +
  theme_light() +
  theme(legend.text.align = 0,
        strip.text = element_text(size = 11, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 11),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top")
pl_pr_2

pl_ill_slide <- ggarrange(pl_mad, pl_dp, pl_pr_2, nrow = 1)
pl_ill_slide
# ggsave(pl_ill_slide, filename = "pl_ill_slide.pdf", width = 18, height = 6)




# EMPIRICAL CONSISTENTCY
hellinger <- function(p, q) sqrt(sum((sqrt(p) - sqrt(q))^2)) / sqrt(2)

n_val <- c(seq(0, 1000, by = 10),
           seq(1050, 2000, by = 50),
           seq(2100, 5000, by = 100),
           seq(5250, 10000, by = 150))
h_dp <- h_mad_1 <- h_mad_ada <- rep(NA, length(n_val)-1)

n <- n_val[2]
idx <- 1
extra_n <- n_val[2]
set.seed(2)
y <- rbinom(extra_n, size = 1, prob = 0.6)
for (i in 1:n) {
  if (y[i] == 0) {
    y[i] <- rpois(1, lambda = 25)
  } else {
    y[i] <- rpois(1, lambda = 60)
  }
}
p_dp <- pred_rg_fixed(y, p_0, alpha, 1, sd_dp, supp, 10)
h_dp[idx] <- hellinger(p_dp, p_true)
p_mad_1 <- pred_rg_fixed(y, p_0, alpha, 1, sd_1, supp, 10)
h_mad_1[idx] <- hellinger(p_mad_1, p_true)
p_mad_ada <- pred_rg_rho(y, p_0, alpha, 3/4, 500, sd_1_34, supp, 10)
h_mad_ada[idx] <- hellinger(p_mad_ada, p_true)
for (j in 3:length(n_val)) {
  n <- n_val[j]
  idx <- j - 1
  extra_n <- n - n_val[idx]
  print(idx)
  
  set.seed(j)
  y_new <- rbinom(n = extra_n, size = 1, prob = 0.6)
  for (i in 1:extra_n) {
    if (y_new[i] == 0) {
      y_new[i] <- rpois(1, lambda = 25)
    } else {
      y_new[i] <- rpois(1, lambda = 60)
    }
  }
  y <- c(y, y_new)
  
  p_dp <- pred_rg_fixed(y_new, p_dp, alpha + n, 1, sd_dp, supp, 10)
  h_dp[idx] <- hellinger(p_dp, p_true)
  
  p_mad_1 <- pred_rg_fixed(y_new, p_mad_1, alpha + n, 1, sd_1, supp, 1)
  h_mad_1[idx] <- hellinger(p_mad_1, p_true)
  
  p_mad_ada <- pred_rg_rho(y_new, p_mad_ada, alpha+n, 3/4, 500, sd_1_34, supp, 1)
  h_mad_ada[idx] <- hellinger(p_mad_ada, p_true)
}

df_h <- data.frame(n_val = n_val[-1],
                   mad_1 = h_mad_1,
                   mad_ada = h_mad_ada,
                   dp = h_dp)
cons_plot <- ggplot(df_h) + 
  geom_line(aes(x = n_val, y = dp, colour = "dp")) +
  geom_line(aes(x = n_val, y = mad_1, colour = "mad_1")) +
  geom_line(aes(x = n_val, y = mad_ada, colour = "mad_ada")) +
  scale_color_manual(name = "", values = c("mad_ada" = "red", "mad_1" = "deepskyblue", "dp" = "black"), 
                     labels = c("mad_ada" = "MAD-ada", "mad_1" = "MAD-1", "dp" = "DP")) +
  labs(x = "Sample size", y = "Hellinger distance with the true data generator") +
  theme_light() +
  theme(legend.text.align = 0,
        strip.text = element_text(size = 11, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top")
cons_plot
# ggsave(cons_plot, filename = "emp_cons.pdf")
