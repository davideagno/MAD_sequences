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

pred_cop <- function(y, p_0, rho, support, nPerm, seed = 1, cores = 5) {
  idx_perm <- matrix(nrow = length(y), ncol = nPerm) 
  set.seed(seed)
  for (p in 1:nPerm) {
    idx_perm[,p] <- sample(1:length(y), size = length(y), replace = FALSE)
  }
  
  mclapply_function <- function(i) {
    y_perm <- y[idx_perm[,i]]
    p_n <- p_0
    for (n in 1:length(y_perm)) {
      y_n <- y_perm[n]
      w_n <- (2 - 1/n) * (n + 1)^(-1)
      p_n <- (1 - w_n + w_n*(1 - rho)) * p_n
      p_n[which(support == y_n)] <- p_n[which(support == y_n)] + w_n*rho
    }
    p_n
  }
  
  out <- mclapply(1:nPerm, mclapply_function, mc.cores = cores)
  out <- rowMeans(matrix(unlist(out), nrow = length(support), ncol = nPerm))
  out
}
preq_logl_cop <- function(rho, y, p_0, support, nPerm, seed = 1, cores = 5) {
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
      w_n <- (2 - 1/n) * (n + 1)^(-1)
      pll <- pll + log(p_n[which(support == y_n)])
      p_n <- (1 - w_n + w_n*(1 - rho)) * p_n
      p_n[which(support == y_n)] <- p_n[which(support == y_n)] + w_n*rho
    }
    pll
  }
  
  out <- mclapply(1:nPerm, mclapply_function, mc.cores = cores)
  mean(unlist(out))
}
pred_res_cop <- function(y, p_n, support, N, B, rho, seed = 1, cores = 5) {
  require(parallel)
  m <- length(y)
  
  mclapply_function <- function(i) {
    prob <- p_n
    for (n in (m+1):N) {
      y_past <- sample(support, 1, prob = prob)
      w_n <- (2 - 1/n) * (n + 1)^(-1)
      prob <- (1 - w_n + w_n*(1 - rho)) * prob
      prob[which(support == y_past)] <- prob[which(support == y_past)] + w_n*rho
    }
    prob
  }
  
  out <- mclapply(1:B, mclapply_function, mc.cores = cores)
  out <- matrix(unlist(out), nrow = length(support), ncol = B)
  out
}

drg <- function(x, mean, sd) {
  D <- length(c(x))
  if (D > 1) {
    out <- rep(NA, D)
    for (i in 1:D) {
      out[i] <- pnorm(x[i]+0.5, mean = mean, sd = sd) - pnorm(x[i]-0.5, mean = mean, sd = sd)
    }
    out <- out / sum(out)
  }
  
  if (D == 1) {
    out <- pnorm(x+0.5, mean = mean, sd = sd) - pnorm(x-0.5, mean = mean, sd = sd)
  }
  
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
pred_rg_dpm <- function(y, p_0, sd_rg, support, nPerm, seed = 1, cores = 5) {
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
      w_n <- (2 - 1/n) * (n + 1)^(-1)
      p_n <- (1 - w_n)*p_n + w_n*ker
    }
    p_n
  }
  
  out <- mclapply(1:nPerm, mclapply_function, mc.cores = cores)
  out <- rowMeans(matrix(unlist(out), nrow = length(support), ncol = nPerm))
  out
}
preq_logl_rg_dpm <- function(sd_rg, y, p_0, support, nPerm, seed = 1, cores = 5) {
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
      w_n <- (2 - 1/n) * (n + 1)^(-1)
      p_n <- (1 - w_n)*p_n + w_n*ker
    }
    pll
  }
  
  out <- mclapply(1:nPerm, mclapply_function, mc.cores = cores)
  mean(unlist(out))
}
pred_res_rg_dpm <- function(y, p_n, support, N, B, sd_rg, seed = 1, cores = 5) {
  require(parallel)
  m <- length(y)
  
  mclapply_function <- function(i) {
    prob <- p_n
    for (n in (m+1):N) {
      y_past <- sample(support, 1, prob = prob)
      ker <- mhk_rg(p_n, y_past, sd_rg, support)
      w_n <- (2 - 1/n) * (n + 1)^(-1)
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


# MAD
opt_mad <- nlminb(start = 1, lower = 1e-10, upper = Inf,
                  objective = function(x) -preq_logl_rg_dpm(x, y, p_0, supp, PP))
opt_mad
sd_mad <- opt_mad$par
p_n_mad <- pred_rg_dpm(y, p_0, sd_mad, supp, PP)
ps_mad <- pred_res_rg_dpm(y, p_n_mad, supp, N, B, sd_mad)
ci_lower_mad <- apply(ps_mad, 1, emp_hpd_lower)
ci_upper_mad <- apply(ps_mad, 1, emp_hpd_upper)


# COPULA - \rho = 0.3
rho_cop_03 <- 0.3
p_n_cop_03 <- pred_cop(y, p_0, rho_cop_03, supp, PP)
ps_cop_03 <- pred_res_cop(y, p_n_cop_03, supp, N, B, rho_cop_03)
ci_lower_cop_03 <- apply(ps_cop_03, 1, emp_hpd_lower)
ci_upper_cop_03 <- apply(ps_cop_03, 1, emp_hpd_upper)


# COPULA - \rho preqll
opt_cop <- nlminb(start = 1, lower = 1e-10, upper = 1-1e-10,
                  objective = function(x) -preq_logl_cop(x, y, p_0, supp, PP))
opt_cop
rho_cop <- opt_cop$par
p_n_cop <- pred_cop(y, p_0, rho_cop, supp, PP)
ps_cop <- pred_res_cop(y, p_n_cop, supp, N, B, rho_cop)
ci_lower_cop <- apply(ps_cop, 1, emp_hpd_lower)
ci_upper_cop <- apply(ps_cop, 1, emp_hpd_upper)


pl_max <- max(ci_upper_mad, ci_upper_cop, ci_upper_cop_03)
df_mad <- data.frame(supp = supp, p_mad = p_n_mad, ci_lower = ci_lower_mad, ci_upper = ci_upper_mad, title = "MAD")
pl_mad <- ggplot(df_mad) + 
  facet_wrap(~ title) +
  geom_segment(aes(x = supp, xend = supp, y = ci_lower, yend = ci_upper, col = "z_ci"), alpha = 0.15, size = 2) +
  geom_segment(aes(x = supp, xend = supp, y = 0, yend = p_mad, col = "pm"), alpha = 1, size = .4) +
  geom_point(aes(x = supp, y = p_true), col = "forestgreen", size = .5) +
  scale_color_manual(name = "", values = c("pm" = "red", "z_ci" = "blue3"), labels = c("pm" = TeX("$P_n$"),
                                                                                       "z_ci" = "Posterior 0.95 CI")) +
  labs(x = "Support", y = "Probability") +
  ylim(c(0, pl_max)) +
  theme_light() +
  theme(legend.position = "top",
        strip.text = element_text(size = 13, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 13),
        axis.text.y = element_text(size = 5))
pl_mad


df_cop_03 <- data.frame(supp = supp, p_cop_03 = p_n_cop_03, ci_lower = ci_lower_cop_03, ci_upper = ci_upper_cop_03, title = "COP-2")
pl_cop_03 <- ggplot(df_cop_03) + 
  facet_wrap(~ title) +
  geom_segment(aes(x = supp, xend = supp, y = ci_lower, yend = ci_upper, col = "z_ci"), alpha = 0.15, size = 2) +
  geom_segment(aes(x = supp, xend = supp, y = 0, yend = p_cop_03, col = "pm"), alpha = 1, size = .4) +
  geom_point(aes(x = supp, y = p_true), col = "forestgreen", size = .5) +
  scale_color_manual(name = "", values = c("pm" = "red", "z_ci" = "blue3"), labels = c("pm" = TeX("$P_n$"),
                                                                                       "z_ci" = "Posterior 0.95 CI")) +
  labs(x = "Support", y = "Probability") +
  ylim(c(0, pl_max)) +
  theme_light() +
  theme(legend.position = "top",
        strip.text = element_text(size = 13, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 13),
        axis.text.y = element_text(size = 5))
pl_cop_03


df_cop <- data.frame(supp = supp, p_cop = p_n_cop, ci_lower = ci_lower_cop, ci_upper = ci_upper_cop, title = "COP-1")
pl_cop <- ggplot(df_cop) + 
  facet_wrap(~ title) +
  geom_segment(aes(x = supp, xend = supp, y = ci_lower, yend = ci_upper, col = "z_ci"), alpha = 0.15, size = 2) +
  geom_segment(aes(x = supp, xend = supp, y = 0, yend = p_cop, col = "pm"), alpha = 1, size = .4) +
  geom_point(aes(x = supp, y = p_true), col = "forestgreen", size = .5) +
  scale_color_manual(name = "", values = c("pm" = "red", "z_ci" = "blue3"), labels = c("pm" = TeX("$P_n$"),
                                                                                       "z_ci" = "Posterior 0.95 CI")) +
  labs(x = "Support", y = "Probability") +
  ylim(c(0, pl_max)) +
  theme_light() +
  theme(legend.position = "top",
        strip.text = element_text(size = 13, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 13),
        axis.text.y = element_text(size = 5))
pl_cop

# Figure S1 of the Supplement
library(ggpubr)
pl_ill_cop <- ggarrange(pl_mad, pl_cop, pl_cop_03, common.legend = TRUE, nrow = 1)
pl_ill_cop
# ggsave(pl_ill_cop, filename = "pl_ill_cop.pdf", width = 16, height = 6)
