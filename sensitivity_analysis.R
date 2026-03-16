rm(list = ls())
setwd("/Users/davideagnoletto/Dropbox/sim_martingales/review")

library(parallel)
library(rmp)
library(dplyr)
library(tidyr)
library(ggplot2)
library(latex2exp)

Rcpp::sourceCpp("illustration_2.cpp")
source("illustration_functions.R")

# Data simulation
nn <- 150
set.seed(1234)
y <- sample(c(2,20), size = nn, replace = TRUE, prob = c(1,1))
for (i in 1:nn) y[i] <- rpois(1, lambda = y[i])
supp <- 0:50
p_true <- 0.5*dpois(supp, lambda = 2) + 0.5*dpois(supp, lambda = 20)
plot(supp, p_true, type = "h")
N <- 2e4
B <- 1000
p_0 <- rep(1/length(supp), length(supp))
alpha <- 1
PP <- 10
y_val <- 9

N_val <- 50000

LL <- 1; N_star <- 1
opt <- nlminb(start = 1, lower = 1e-10, upper = Inf,
              objective = function(x) -preq_ll(x, y, p_0, alpha, LL, N_star, supp, PP))
sd <- opt$par
pred <- pred_mad(y, p_0, alpha, LL, N_star, sd, supp, PP)
ps <- pred_res(y, pred, alpha, LL, N_star, supp, N, B, sd)
ps_th_1 <- apply(ps, 2, function(x) sum(x[1:y_val]))
ps_conv_1 <- pred_res_conv_val(y_val, y, pred, alpha, LL, N_star, supp, N = N_val, B = 10, 
                           sd_rg = sd, seed = 12345)



LL <- 3/4; N_star <- 20
opt <- nlminb(start = 1, lower = 1e-10, upper = Inf,
              objective = function(x) -preq_ll(x, y, p_0, alpha, LL, N_star, supp, PP))
sd <- opt$par
pred <- pred_mad(y, p_0, alpha, LL, N_star, sd, supp, PP)
ps <- pred_res(y, pred, alpha, LL, N_star, supp, N, B, sd)
ps_th_34_20 <- apply(ps, 2, function(x) sum(x[1:y_val]))
ps_conv_34_20 <- pred_res_conv_val(y_val, y, pred, alpha, LL, N_star, supp, N = N_val, B = 10, 
                               sd_rg = sd, seed = 12345)

LL <- 2/3; N_star <- 20
opt <- nlminb(start = 1, lower = 1e-10, upper = Inf,
              objective = function(x) -preq_ll(x, y, p_0, alpha, LL, N_star, supp, PP))
sd <- opt$par
pred <- pred_mad(y, p_0, alpha, LL, N_star, sd, supp, PP)
ps <- pred_res(y, pred, alpha, LL, N_star, supp, N, B, sd)
ps_th_23_20 <- apply(ps, 2, function(x) sum(x[1:y_val]))
ps_conv_23_20 <- pred_res_conv_val(y_val, y, pred, alpha, LL, N_star, supp, N = N_val, B = 10, 
                                   sd_rg = sd, seed = 12345)

LL <- 0.6; N_star <- 20
opt <- nlminb(start = 1, lower = 1e-10, upper = Inf,
              objective = function(x) -preq_ll(x, y, p_0, alpha, LL, N_star, supp, PP))
sd <- opt$par
pred <- pred_mad(y, p_0, alpha, LL, N_star, sd, supp, PP)
ps <- pred_res(y, pred, alpha, LL, N_star, supp, N, B, sd)
ps_th_6_20 <- apply(ps, 2, function(x) sum(x[1:y_val]))
ps_conv_6_20 <- pred_res_conv_val(y_val, y, pred, alpha, LL, N_star, supp, N = N_val, B = 10, 
                                   sd_rg = sd, seed = 12345)

LL <- 3/4; N_star <- 500
opt <- nlminb(start = 1, lower = 1e-10, upper = Inf,
              objective = function(x) -preq_ll(x, y, p_0, alpha, LL, N_star, supp, PP))
sd <- opt$par
pred <- pred_mad(y, p_0, alpha, LL, N_star, sd, supp, PP)
ps <- pred_res(y, pred, alpha, LL, N_star, supp, N, B, sd)
ps_th_34_500 <- apply(ps, 2, function(x) sum(x[1:y_val]))
ps_conv_34_500 <- pred_res_conv_val(y_val, y, pred, alpha, LL, N_star, supp, N = N_val, B = 10, 
                                   sd_rg = sd, seed = 12345)

LL <- 2/3; N_star <- 500
opt <- nlminb(start = 1, lower = 1e-10, upper = Inf,
              objective = function(x) -preq_ll(x, y, p_0, alpha, LL, N_star, supp, PP))
sd <- opt$par
pred <- pred_mad(y, p_0, alpha, LL, N_star, sd, supp, PP)
ps <- pred_res(y, pred, alpha, LL, N_star, supp, N, B, sd)
ps_th_23_500 <- apply(ps, 2, function(x) sum(x[1:y_val]))
ps_conv_23_500 <- pred_res_conv_val(y_val, y, pred, alpha, LL, N_star, supp, N = N_val, B = 10, 
                                    sd_rg = sd, seed = 12345)

LL <- 0.6; N_star <- 500
opt <- nlminb(start = 1, lower = 1e-10, upper = Inf,
              objective = function(x) -preq_ll(x, y, p_0, alpha, LL, N_star, supp, PP))
sd <- opt$par
pred <- pred_mad(y, p_0, alpha, LL, N_star, sd, supp, PP)
ps <- pred_res(y, pred, alpha, LL, N_star, supp, N, B, sd)
ps_th_6_500 <- apply(ps, 2, function(x) sum(x[1:y_val]))
ps_conv_6_500 <- pred_res_conv_val(y_val, y, pred, alpha, LL, N_star, supp, N = N_val, B = 10, 
                                    sd_rg = sd, seed = 12345)

LL <- 3/4; N_star <- 1000
opt <- nlminb(start = 1, lower = 1e-10, upper = Inf,
              objective = function(x) -preq_ll(x, y, p_0, alpha, LL, N_star, supp, PP))
sd <- opt$par
pred <- pred_mad(y, p_0, alpha, LL, N_star, sd, supp, PP)
ps <- pred_res(y, pred, alpha, LL, N_star, supp, N, B, sd)
ps_th_34_1000 <- apply(ps, 2, function(x) sum(x[1:y_val]))
ps_conv_34_1000 <- pred_res_conv_val(y_val, y, pred, alpha, LL, N_star, supp, N = N_val, B = 10, 
                                    sd_rg = sd, seed = 12345)

LL <- 2/3; N_star <- 1000
opt <- nlminb(start = 1, lower = 1e-10, upper = Inf,
              objective = function(x) -preq_ll(x, y, p_0, alpha, LL, N_star, supp, PP))
sd <- opt$par
pred <- pred_mad(y, p_0, alpha, LL, N_star, sd, supp, PP)
ps <- pred_res(y, pred, alpha, LL, N_star, supp, N, B, sd)
ps_th_23_1000 <- apply(ps, 2, function(x) sum(x[1:y_val]))
ps_conv_23_1000 <- pred_res_conv_val(y_val, y, pred, alpha, LL, N_star, supp, N = N_val, B = 10, 
                                     sd_rg = sd, seed = 12345)

LL <- 0.6; N_star <- 1000
opt <- nlminb(start = 1, lower = 1e-10, upper = Inf,
              objective = function(x) -preq_ll(x, y, p_0, alpha, LL, N_star, supp, PP))
sd <- opt$par
pred <- pred_mad(y, p_0, alpha, LL, N_star, sd, supp, PP)
ps <- pred_res(y, pred, alpha, LL, N_star, supp, N, B, sd)
ps_th_6_1000 <- apply(ps, 2, function(x) sum(x[1:y_val]))
ps_conv_6_1000 <- pred_res_conv_val(y_val, y, pred, alpha, LL, N_star, supp, N = N_val, B = 10, 
                                     sd_rg = sd, seed = 12345)

# Figure S3 of the Supplement
df_var <- data.frame(p1 = rep(ps_th_1, 3),
                     p34 = c(ps_th_34_1000, ps_th_34_500, ps_th_34_20),
                     p23 = c(ps_th_23_1000, ps_th_23_500, ps_th_23_20),
                     p6 = c(ps_th_6_1000, ps_th_6_500, ps_th_6_20),
                     title = rep(c(1000, 500, 20), each = length(ps_th_1)))
plot_var <- ggplot(df_var) + 
  facet_wrap(~ title, labeller = label_bquote(N["*"] == .(title))) +
  stat_density(aes(x = p1, col = "1"), linetype = 1, geom = "line", position = "identity", adjust = 2) +
  stat_density(aes(x = p34, col = "2"), linetype = 1, geom = "line", position = "identity", adjust = 2.5) +
  stat_density(aes(x = p23, col = "3"), linetype = 1, geom = "line", position = "identity", adjust = 2.5) +
  stat_density(aes(x = p6, col = "4"), linetype = 1, geom = "line", position = "identity", adjust = 2.5) +
  scale_color_manual(name = "lambda",
                     values = c("1" = "#1170aa",
                                "2" = "#fc7d0b",
                                "3" = "#57606c",
                                "4" = "#a3cce9"),
                     labels = c("1" = TeX(r'($\lambda = 1$)'), 
                                "2" = TeX(r'($\lambda = 3/4$)'),
                                "3" = TeX(r'($\lambda = 2/3$)'),
                                "4" = TeX(r'($\lambda = 0.6$)')
                     )) +
  xlab("Probability") + ylab("Density") +
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
plot_var



df_conv <- data.frame(conv_1 = rep(ps_conv_1[,1], 3),
                      conv_34 = c(ps_conv_34_1000[,1], ps_conv_34_500[,1], ps_conv_34_20[,1]),
                      conv_23 = c(ps_conv_23_1000[,1], ps_conv_23_500[,1], ps_conv_23_20[,1]),
                      conv_6 = c(ps_conv_6_1000[,1], ps_conv_6_500[,1], ps_conv_6_20[,1]),
                      fow_sam = rep(c(1:NROW(ps_conv_1)), 3),
                      title = c(rep(c(1000, 500, 20), each = NROW(ps_conv_1))))
plot_conv <- ggplot(df_conv) + 
  facet_wrap(~ title, labeller = label_bquote(N["*"] == .(title))) +
  geom_line(aes(x = fow_sam, y = conv_6, color = "c6")) +
  geom_line(aes(x = fow_sam, y = conv_23, color = "c23")) +
  geom_line(aes(x = fow_sam, y = conv_34, color = "c34")) +
  geom_line(aes(x = fow_sam, y = conv_1, color = "c1")) +
  scale_color_manual(values = c("c1" = "#1170aa", 
                                "c23" = "#57606c", 
                                "c34" = "#fc7d0b", 
                                "c6" = "#a3cce9"),
                     labels = c("c1" = TeX(r'($\lambda = 1$)'), 
                                "c34" = TeX(r'($\lambda = 3/4$)'), 
                                "c23" = TeX(r'($\lambda = 2/3$)'),
                                "c6" = TeX(r'($\lambda = 0.6$)'))) +
  xlab("Number of forward samples") + ylab("Probability") +
  theme_light() + xlim(0,20000) +
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
plot_conv


pl_sens <- ggarrange(plot_conv, plot_var, common.legend = TRUE, nrow = 2)
pl_sens
# ggsave(pl_sens, filename = "pl_sens.pdf", width = 11, height = 9)

