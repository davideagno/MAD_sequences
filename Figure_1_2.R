
rm(list = ls())

library(TeachingDemos)
library(ggplot2)
library(tmvtnorm)
library(latex2exp)
library(ggthemes)

rn <- function(alpha, n) {
  out <- (n+1):2e6
  out <- out + alpha
  1 / sum(1 / out ^ 2)
} 
drg <- function(x, mean, sd, left = -0.5) {
  require(tmvtnorm)
  D <- length(c(x))
  out <- rep(NA, D)
  for (i in 1:D) {
    out[i] <- ptmvnorm(lowerx = x[i]-0.5, upperx = x[i]+0.5, mean = mean, sigma = sd^2, lower = left)[1]
  }
  out
}
drg_rev <- function(x, mean, sd, left = -0.5) {
  require(tmvtnorm)
  D <- length(c(mean))
  out <- rep(NA, D)
  for (i in 1:D) {
    out[i] <- ptmvnorm(lowerx = x-0.5, upperx = x+0.5, mean = mean[i], sigma = sd^2, lower = left)[1]
  }
  out
}
var_p1 <- function(alpha, sd, x, p_0, support) {
  k0_2 <- sum(drg_rev(x, support, sd)^2)
  num <- p_0*(k0_2 - p_0)
  den <- (alpha + 1)^2
  num / den
}

mhk <- function(p_n, y, sd, support) {
  out <- rep(NA, length(support))
  idx_obs <- y + 1
  log_k0 <- log(drg(support, y, sd))
  log_k0_n <- log(drg_rev(y, support, sd))
  prob_acc <- exp(log(p_n) + log_k0_n - log_k0 - log(p_n[idx_obs]))
  prob_acc[is.nan(prob_acc)] <- 0
  prob_acc <- pmin(prob_acc, 1)
  out <- prob_acc * exp(log_k0)
  out[idx_obs] <- out[idx_obs] + 1 - sum(out)
  out
}
pred <- function(y, alpha, p_0, sd, support) {
  p_n <- p_0
  for (n in 1:length(y)) {
    ker <- mhk(p_n, y[n], sd, support)
    w_n <- 1 / (alpha + n)
    p_n <- (1 - w_n) * p_n + w_n * ker
  }
  p_n
}
pred_perm <- function(y, alpha, p_0, sd, support, nPerm, cores, seed = 1) {
  idx <- matrix(nrow = length(y), ncol = nPerm)
  set.seed(seed)
  for (i in 1:nPerm) {
    idx[,i] <- sample(1:length(y), size = length(y), replace = FALSE)
  }
  
  mclapply_function <- function(b) {
    pred(y[idx[,b]], alpha, p_0, sd, support)
  }
  out <- mclapply(1:nPerm, mclapply_function, mc.cores = cores)
  out <- matrix(unlist(out), ncol = nPerm)
  rowMeans(out)
}
var_pn <- function(sd, x, alpha, p_n, N, support) {
  num <- sum((mhk(p_n, x, sd, support)^2) * p_n) - (p_n[x+1])^2
  den <- (alpha + N + 1)^2
  num / den
}
cov_pn <- function(sd, x1, x2, alpha, p_n, N, support) {
  num <- sum(mhk(p_n, x1, sd, support) * mhk(p_n, x2, sd, support) * p_n) - (p_n[x1+1]*p_n[x2+1])
  den <- (alpha + N + 1)^2
  num / den
}
cor_pn <- function(sd, x1, x2, alpha, p_n, N, support) {
  num <- sum(mhk(p_n, x1, sd, support) * mhk(p_n, x2, sd, support) * p_n) - (p_n[x1+1]*p_n[x2+1])
  den_1 <- sqrt(sum((mhk(p_n, x1, sd, support)^2) * p_n) - (p_n[x1+1])^2)
  den_2 <- sqrt(sum((mhk(p_n, x2, sd, support)^2) * p_n) - (p_n[x2+1])^2)
  num / (den_1 * den_2)
}

supp <- 0:100
alpha <- 1
p_0 <- dpois(supp, lambda = 30) + dpois(supp, lambda = 55)
p_0 <- p_0 / sum(p_0)
plot(supp, p_0, type = "h")

# Kernels
y_obs <- 40
sd_grid <- c(0.5, 2, 10, 50)
kns <- matrix(NA, nrow = length(supp), ncol = length(sd_grid))
for (i in 1:length(sd_grid)) {
  kns[,i] <- mhk(p_0, y_obs, sd_grid[i], supp)
}
df_kk <- data.frame(supp = rep(supp, length(sd_grid)),
                    kernel = c(kns),
                    k_0 = c(drg(supp, mean = y_obs, sd = sd_grid[1]),
                            drg(supp, mean = y_obs, sd = sd_grid[2]),
                            drg(supp, mean = y_obs, sd = sd_grid[3]),
                            drg(supp, mean = y_obs, sd = sd_grid[4])),
                    p_0 = rep(p_0, length(sd_grid)),
                    sigma = rep(sd_grid, each = length(supp)))
kk <- ggplot(df_kk) +
  facet_wrap(~ sigma, nrow = 2, labeller = label_bquote(sigma == .(sigma)) ) +
  geom_segment(aes(x = supp+0.2, xend = supp+0.2, y = 0, yend = k_0, color = "k0")) +
  geom_segment(aes(x = supp-0.2, xend = supp-0.2, y = 0, yend = p_0, color = "pred")) +
  geom_segment(aes(x = supp, xend = supp, y = 0, yend = kernel, color = "ker")) +
  scale_color_tableau(name = "",
                      palette = "Color Blind",
                      labels = c("k0" = TeX(r'($k_{*}(y \;|\; y_n,\sigma)$)'), 
                                 "ker" = TeX(r'($k_{n-1}(y \;|\; y_n,\sigma)$)'),
                                 "pred" = TeX(r'($p_{n-1}(y)$)'))) +
  labs(x = "Support", y = "Probability") +
  theme_light() +
  theme(legend.position = "top",
        strip.text = element_text(size = 10, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"))
kk
# ggsave(kk, filename = "ker_plot.pdf", width = 8, height = 7)


kk_slide <- ggplot(df_kk) +
  facet_wrap(~ sigma, nrow = 1, labeller = label_bquote(sigma == .(sigma)) ) +
  geom_segment(aes(x = supp+0.2, xend = supp+0.2, y = 0, yend = k_0, color = "k0")) +
  geom_segment(aes(x = supp-0.2, xend = supp-0.2, y = 0, yend = p_0, color = "pred")) +
  geom_segment(aes(x = supp, xend = supp, y = 0, yend = kernel, color = "ker")) +
  scale_color_tableau(name = "",
                      palette = "Color Blind",
                      labels = c("k0" = TeX(r'($k_{*}(y \;|\; y_n,\sigma)$)'), 
                                 "ker" = TeX(r'($k_{n-1}(y \;|\; y_n,\sigma)$)'),
                                 "pred" = TeX(r'($p_{n-1}(y)$)'))) +
  labs(x = "Support", y = "Probability") +
  theme_light() +
  theme(legend.position = "top",
        strip.text = element_text(size = 10, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"))
kk_slide
# ggsave(kk_slide, filename = "ker_plot_slide.pdf", width = 16, height = 4)



# Variance
show_col(tableau_color_pal("Color Blind")(10))
nn <- 50
sd_val <- seq(1e-100, 40, length = 200)
vv_20 <- vv_40 <- vv_60 <- vv_80 <- rep(NA, length(sd_val))
plot(supp, p_0, type = "h")
for (i in 1:length(sd_val)) {
  vv_20[i] <- var_pn(sd_val[i], 20, alpha, p_0, nn, supp)
  vv_40[i] <- var_pn(sd_val[i], 40, alpha, p_0, nn, supp)
  vv_60[i] <- var_pn(sd_val[i], 60, alpha, p_0, nn, supp)
  vv_80[i] <- var_pn(sd_val[i], 80, alpha, p_0, nn, supp)
}
df_vv <- data.frame(sd = sd_val, var_20 = vv_20, var_40 = vv_40, var_60 = vv_60,
                    var_80 = vv_80, title = "Induced variance")
vv <- ggplot(df_vv) + 
  facet_wrap(~ title) +
  geom_line(aes(x = sd, y = var_80, color = "80")) +
  geom_line(aes(x = sd, y = var_60, color = "60")) +
  geom_line(aes(x = sd, y = var_40, color = "40")) +
  geom_line(aes(x = sd, y = var_20, color = "20")) +
  scale_color_manual(name = "",
                     values = c("20" = "#1170aa",
                                "40" = "#a3cce9",
                                "60" = "#a3acb9",
                                "80" = "#5fa2ce"),
                     labels = c("20" = TeX(r'($p_n(20)$)'),
                                "40" = TeX(r'($p_n(40)$)'),
                                "60" = TeX(r'($p_n(60)$)'),
                                "80" = TeX(r'($p_n(80)$)')
                     )) +
  xlab(expression(sigma)) + ylab("Variance") +
  theme_light() +
  theme(legend.position = "top",
        axis.text.y = element_text(size = 4.5),
        axis.title.y = element_text(size = 9),
        strip.text = element_text(size = 11, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"))
vv


# Induced correlation
pp <- p_0
cc_41 <- cc_60 <- cc_10 <- rep(NA, length(sd_val))
for (i in 1:length(sd_val)) {
  cc_41[i] <- cor_pn(sd_val[i], 40, 41, alpha, pp, nn, supp)
  cc_60[i] <- cor_pn(sd_val[i], 40, 60, alpha, pp, nn, supp)
  cc_10[i] <- cor_pn(sd_val[i], 40, 10, alpha, pp, nn, supp)
}
df_cc <- data.frame(sd = sd_val, cor_41 = cc_41, cor_60 = cc_60, cor_10 = cc_10,
                    title = "Induced correlation")
cc <- ggplot(df_cc) + 
  facet_wrap(~ title) +
  geom_line(aes(x = sd, y = cor_41, color = "41")) +
  geom_line(aes(x = sd, y = cor_60, color = "60")) +
  geom_line(aes(x = sd, y = cor_10, color = "10")) +
  scale_color_manual(name = "",
                     values = c("10" = "#ffbc79",
                                "41" = "#fc7d0b",
                                "60" = "#57606c"),
                     labels = c("10" = TeX(r'($p_{n}(40), p_{n}(10)$)'),
                                "41" = TeX(r'($p_{n}(40), p_{n}(41)$)'),
                                "60" = TeX(r'($p_{n}(40), p_{n}(60)$)')
                     )) +
  ylim(c(-.5, 1)) +
  xlab(expression(sigma)) + ylab("Correlation") +
  theme_light() +
  theme(legend.position = "top",
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 9),
        strip.text = element_text(size = 11, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"))
cc

cv <- ggpubr::ggarrange(vv, cc, ncol = 2)
cv
# ggsave(cv, filename = "cor_var_plot.pdf", width = 11, height = 6)
