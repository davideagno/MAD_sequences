
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
nn <- 20
set.seed(123)
y <- rpois(nn, lambda = 15)
supp <- 0:30
N <- 1000
B <- 500
p_0 <- rep(1/length(supp), length(supp))
alpha <- 1
p_true <- dpois(supp, lambda = 15)
PP <- 10


# MAD-1
LL <- 1
opt_1 <- nlminb(start = 1, lower = 1e-10, upper = Inf,
                objective = function(x) -preq_logl_rg_fixed(x, y, p_0, alpha, LL, supp, PP))
opt_1
sd_1 <- opt_1$par
p_n_1 <- pred_rg_fixed(y, p_0, alpha, LL, sd_1, supp, PP)


# DP
LL <- 1
sd_dp <- 1e-10
p_n_dp <- pred_rg_fixed(y, p_0, alpha, LL, sd_dp, supp, PP)


freq_obs <- rep(NA, length(supp))
for (i in 1:length(supp)) freq_obs[i] <- sum(y == supp[i])
freq_obs <- freq_obs / sum(freq_obs)


df <- data.frame(supp = rep(supp, 2),
                 p_mad = c(p_n_1, p_n_dp),
                 setting = rep(c("MAD", "DP"), each = length(supp)),
                 freq_obs = rep(freq_obs, 2))
ggplot(df) +
  facet_wrap(~ setting) +
  geom_segment(aes(x = supp, xend = supp, y = 0, yend = p_mad, linetype = "pm"), col = "red", alpha = 1, size = .4) +
  geom_segment(aes(x = supp+0.25, xend = supp+0.25, y = 0, yend = freq_obs, , linetype = "obs"), col = "gray60", alpha = 1, size = .4) +
  scale_linetype_manual(name = "", values = c("pm" = "solid", "obs" = "dashed"), labels = c("pm" = "Posterior mean", "obs" = "Observed frequency")) + 
  scale_color_manual(name = "", values = c("pm" = "red", "obs" = "gray10"), labels = c("pm" = "Posterior mean", "obs" = "Observed frequency")) +
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
