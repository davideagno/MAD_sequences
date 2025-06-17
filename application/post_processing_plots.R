
rm(list = ls())

library(Rcpp)
library(RcppArmadillo)
library(SparseArray)
library(ggplot2)
library(ggcorrplot)
library(latex2exp)

load("mad_object_a.RData")
p_n <- pp$p_n
nz_log <- as.numeric(is_nonzero(p_n))
idx_nz <- c(1:length(nz_log))[which(nz_log == 1)]
load("pred_a.RData")



# CI predictive conditional - CROWS across TEMPERATURES
ps_cc_t1 <- ps_cc_t2 <- ps_cc_t3 <- matrix(nrow = NCOL(ps), ncol = 51)
for (s in 1:NCOL(ps)) {
  print(s)
  post <- SVT_SparseArray(dim = dim(p_n))
  post[idx_nz] <- ps[,s]

  pc_1 <- pc_2 <- pc_3 <- rep(NA, 51)
  for (j in 1:51) {
    pc_1[j] <- sum(post[j,,,,1,])
    pc_2[j] <- sum(post[j,,,,2,])
    pc_3[j] <- sum(post[j,,,,3,])
  }
  ps_cc_t1[s,] <- pc_1 / sum(pc_1)
  ps_cc_t2[s,] <- pc_2 / sum(pc_2)
  ps_cc_t3[s,] <- pc_3 / sum(pc_3)
}
# save(ps_cc_t1, ps_cc_t2, ps_cc_t3, file = "pred_res_crows_ax_temperatures.RData")
# load("pred_res_crows_ax_temperatures.RData")

ci_cc_t1 <- ci_cc_t2 <- ci_cc_t3 <- matrix(nrow = 2, ncol = 51)
for (i in 1:51) {
  ci_cc_t1[,i] <- robustBLME::hpd(ps_cc_t1[,i])
  ci_cc_t2[,i] <- robustBLME::hpd(ps_cc_t2[,i])
  ci_cc_t3[,i] <- robustBLME::hpd(ps_cc_t3[,i])
}

pc_cc_t1 <- pred_cond_mad_temp(pp$p_n, pp$p_az, 1, 51)[,1]
pc_cc_t2 <- pred_cond_mad_temp(pp$p_n, pp$p_az, 2, 51)[,1]
pc_cc_t3 <- pred_cond_mad_temp(pp$p_n, pp$p_az, 3, 51)[,1]

pc_cc_t1 <- apply(ps_cc_t1, 2, mean)
pc_cc_t2 <- apply(ps_cc_t2, 2, mean)
pc_cc_t3 <- apply(ps_cc_t3, 2, mean)

df_cct <- data.frame(support = rep(0:50, 3),
                     pred = c(pc_cc_t1, pc_cc_t2, pc_cc_t3),
                     ci_lower = c(ci_cc_t1[1,], ci_cc_t2[1,], ci_cc_t3[1,]),
                     ci_upper = c(ci_cc_t1[2,], ci_cc_t2[2,], ci_cc_t3[2,]),
                     temp = rep(c("Cold temperature", "Mild temperature", "Warm temperature"), each = 51))
pl_cct <- ggplot(df_cct) + 
  facet_wrap(~ temp) +
  geom_segment(aes(x = support, xend = support, y = ci_lower, yend = ci_upper, col = "z_ci"), alpha = 0.15, size = 4) +
  geom_segment(aes(x = support, xend = support, y = 0, yend = pred, col = "pm"), alpha = 1, size = 1) +
  scale_color_manual(name = "", values = c("pm" = "red", "z_ci" = "blue3"), labels = c("pm" = TeX("$P_n$"),
                                                                                       "z_ci" = "0.95-credible interval")) +
  labs(x = "Number of carrion crows per site", y = "Probability") +
  xlim(c(0, 25)) +
  theme_light() +
  theme(legend.text.align = 0,
        strip.text = element_text(size = 11, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top")
pl_cct
# ggsave(pl_cct, filename = "plot_cct.pdf", width = 12, height = 6)



# CORVIDS in URBAN HABITAT
ps_urb <- matrix(nrow = NCOL(ps), ncol = 4)
for (i in 1:NCOL(ps)) {
  print(i)
  post <- SVT_SparseArray(dim = dim(p_n))
  post[idx_nz] <- ps[,i]

  pc_1 <- pc_2 <- pc_3 <- pc_4 <- rep(NA, 51)
  for (j in 1:51) {
    pc_1[j] <- sum(post[j,,,,,4])
    pc_2[j] <- sum(post[,j,,,,4])
    pc_3[j] <- sum(post[,,j,,,4])
    pc_4[j] <- sum(post[,,,j,,4])
  }
  pc_1 <- pc_1 / sum(pc_1)
  pc_2 <- pc_2 / sum(pc_2)
  pc_3 <- pc_3 / sum(pc_3)
  pc_4 <- pc_4 / sum(pc_4)

  ps_urb[i,1] <- sum(0:50 * pc_1)
  ps_urb[i,2] <- sum(0:50 * pc_2)
  ps_urb[i,3] <- sum(0:50 * pc_3)
  ps_urb[i,4] <- sum(0:50 * pc_4)
}
# save(ps_urb, file = "pred_res_corvids_in_urban.RData")
# load("pred_res_corvids_in_urban.RData")
colMeans(ps_urb)
round(robustBLME::hpd(ps_urb[,1]), 2)
round(robustBLME::hpd(ps_urb[,2]), 2)
round(robustBLME::hpd(ps_urb[,3]), 2)
round(robustBLME::hpd(ps_urb[,4]), 2)

library(ggplot2)
library(ggpubr)
library(ggthemes)

df_urb <- data.frame(c1 = ps_urb[,1], c2 = ps_urb[,2], c3 = ps_urb[,3],
                     c4 = ps_urb[,4], title = "Corvids in Urban habitat")
pl_urb <- ggplot(df_urb) + 
  facet_wrap(~ title) +
  stat_density(aes(x = c1, col = "crow"),
               geom = "line", position = "identity", adjust = 1.8) +
  stat_density(aes(x = c2, col = "magpie"),
               geom = "line", position = "identity", adjust = 1.8) +
  stat_density(aes(x = c3, col = "raven"),
               geom = "line", position = "identity", adjust = 3) +
  stat_density(aes(x = c4, col = "jay"),
               geom = "line", position = "identity", adjust = 3) +
  scale_color_tableau( name = "", palette = "Color Blind",
                       labels = c("crow" = "Crow", 
                                  "magpie" = "Magpie", 
                                  "raven" = "Raven",
                                  "jay" = "Jay"
                       )) +
  xlim(c(0, 20)) + xlab("Expected number of corvids per location") +
  ylab("Density") +
  theme_light() +
  theme(legend.text.align = 0,
        strip.text = element_text(size = 11, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top")
pl_urb



pl_appl_slide_2 <- ggarrange(pl_urb, pl_cct_cold)
pl_appl_slide_2
# ggsave(pl_appl_slide_2, filename = "pl_appl_slide_2.pdf", width = 12, height = 6)




# MAGPIES across HABITATS
ps_emh <- matrix(nrow = NCOL(ps), ncol = 5)
for (i in 1:NCOL(ps)) {
  print(i)
  post <- SVT_SparseArray(dim = dim(p_n))
  post[idx_nz] <- ps[,i]

  pc_1 <- pc_2 <- pc_3 <- pc_4 <- pc_5 <- rep(NA, 51)
  for (j in 1:51) {
    pc_1[j] <- sum(post[,j,,,,1])
    pc_2[j] <- sum(post[,j,,,,2])
    pc_3[j] <- sum(post[,j,,,,3])
    pc_4[j] <- sum(post[,j,,,,4])
    pc_5[j] <- sum(post[,j,,,,5])
  }
  pc_1 <- pc_1 / sum(pc_1)
  pc_2 <- pc_2 / sum(pc_2)
  pc_3 <- pc_3 / sum(pc_3)
  pc_4 <- pc_4 / sum(pc_4)
  pc_5 <- pc_5 / sum(pc_5)

  ps_emh[i,1] <- sum(0:50 * pc_1)
  ps_emh[i,2] <- sum(0:50 * pc_2)
  ps_emh[i,3] <- sum(0:50 * pc_3)
  ps_emh[i,4] <- sum(0:50 * pc_4)
  ps_emh[i,5] <- sum(0:50 * pc_5)
}
# save(ps_emh, file = "pred_res_magpies_ax_habitats.RData")
# load("pred_res_magpies_ax_habitats.RData")
colMeans(ps_emh)
round(robustBLME::hpd(ps_emh[,1]), 2)
round(robustBLME::hpd(ps_emh[,2]), 2)
round(robustBLME::hpd(ps_emh[,3]), 2)
round(robustBLME::hpd(ps_emh[,4]), 2)
round(robustBLME::hpd(ps_emh[,5]), 2)

df_emh <- data.frame(t1 = ps_emh[,1], t2 = ps_emh[,2], t3 = ps_emh[,3], t4 = ps_emh[,4], t5 = ps_emh[,5],
                     title = "Eurasian magpies across habitat")
pl_emh <- ggplot(df_emh) + 
  facet_wrap(~ title) +
  stat_density(aes(x = t1, col = "h1"), 
               geom = "line", position = "identity", adjust = 2.5) +
  stat_density(aes(x = t2, col = "h2"), 
               geom = "line", position = "identity", adjust = 4.5) +
  stat_density(aes(x = t3, col = "h3"), 
               geom = "line", position = "identity", adjust = 2.5) +
  stat_density(aes(x = t4, col = "h4"), 
               geom = "line", position = "identity", adjust = 2) +
  stat_density(aes(x = t5, col = "h5"), 
               geom = "line", position = "identity", adjust = 2.5) +
  scale_color_tableau( name = "", palette = "Color Blind",
                       labels = c("h1" = "Broadleaved", 
                                  "h2" = "Coniferous",
                                  "h3" = "Open", 
                                  "h4" = "Urban", 
                                  "h5" = "Wetlands"
                       )) +
  xlim(c(0,15)) + xlab("Expected number of magpies per location") +
  ylab("Density") +
  theme_light() +
  theme(legend.text.align = 0,
        strip.text = element_text(size = 11, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 10),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top")
pl_emh
# ggsave(pl_emh, filename = "pl_supp.pdf", width = 6, height = 6)



# PS for correlation
load("ps_correlation.RData")
robustBLME::hpd(cor_cc_em)
df_cor <- data.frame(c1 = cor_cc_em, c2 = cor_cc_cr, c3 = cor_cc_ej, c4 = cor_em_cr, 
                     c5 = cor_em_ej, c6 = cor_cr_ej, title = "Correlation between species")
pl_cor <- ggplot(df_cor) + 
  facet_wrap(~ title) +
  stat_density(aes(x = c1, col = "c1"), linetype = 1, 
               geom = "line", position = "identity", adjust = 3) +
  stat_density(aes(x = c2, col = "c2"), linetype = 1,
               geom = "line", position = "identity", adjust = 3) +
  stat_density(aes(x = c3, col = "c3"), linetype = 1,
               geom = "line", position = "identity", adjust = 3) +
  stat_density(aes(x = c4, col = "c4"), linetype = 1,
               geom = "line", position = "identity", adjust = 3) +
  stat_density(aes(x = c5, col = "c5"), linetype = 1,
               geom = "line", position = "identity", adjust = 3) +
  stat_density(aes(x = c6, col = "c6"), linetype = 1,
               geom = "line", position = "identity", adjust = 3) +
  scale_color_tableau( name = "", palette = "Color Blind",
                       labels = c("c1" = "CC-EM", 
                                  "c2" = "CC-CR",
                                  "c3" = "CC-EJ", 
                                  "c4" = "EM-CR", 
                                  "c5" = "EM-EJ",
                                  "c6" = "CR-EJ"
                       )) +
  guides(colour = guide_legend(nrow = 1)) +
  xlim(c(-1, 1)) + xlab("Correlation") +
  ylab("Density") +
  theme_light() +
  theme(legend.text.align = 0,
        strip.text = element_text(size = 11, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 10),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top")
pl_cor

library(ggpubr)
pl_main <- ggarrange(pl_cor, pl_urb)
pl_main
# ggsave(pl_main, filename = "pl_main.pdf", width = 11, height = 6)


# COVARIATES MAP
library(sf)
library(giscoR)
library(dplyr)
dat <- read.csv("sites_loc_09.csv")
fin_map <- gisco_get_countries(resolution = "10", country = "FIN",
                               epsg = "3035") %>% mutate(res = "20M")
df_map_temp <- data.frame(x_site = dat$x_site * 1000,
                          y_site = dat$y_site * 1000,
                          hab = dat$hab_site,
                          temp = dat$temp_site,
                          border = rep("gray30", NROW(dat)),
                          lab = rep("Temperature", NROW(dat)))
df_map_temp$lab <- as.factor(df_map_temp$lab)
map_temp <- ggplot(fin_map) +
  geom_sf(fill = "ghostwhite") +
  facet_wrap(~ lab) +
  geom_point(data = df_map_temp, aes(x = x_site, y = y_site, fill = temp), 
             pch = 23, size = 4, colour = "gray30") + 
  scale_fill_manual(name = "",
                    values = c("Cold" = "blue",
                               "Mild" = "grey50",
                               "Warm" = "red3")) +
  scale_size_area(breaks = c(2,4,6), label = c(1,2,3)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_light() +
  theme(strip.text = element_text(size = 12, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 7),
        legend.key.size = unit(11, 'pt'),
        legend.title = element_text(size = 10),
        legend.position = "top")

df_map_hab <- data.frame(x_site = dat$x_site * 1000,
                         y_site = dat$y_site * 1000,
                         hab = dat$hab_site,
                         temp = dat$temp_site,
                         border = rep("gray30", NROW(dat)),
                         lab = rep("Habitat", NROW(dat)))
df_map_hab$lab <- as.factor(df_map_hab$lab)
map_hab <- ggplot(fin_map) +
  geom_sf(fill = "ghostwhite") +
  facet_wrap(~ lab) +
  geom_point(data = df_map_hab, aes(x = x_site, y = y_site, fill = hab), 
             pch = 23, size = 4, colour = "gray30") + 
  scale_fill_tableau( name = "", palette = "Color Blind") +
  scale_size_area(breaks = c(2,4,6), label = c(1,2,3)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_light() +
  theme(strip.text = element_text(size = 12, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 7),
        legend.key.size = unit(11, 'pt'),
        legend.title = element_text(size = 10),
        legend.position = "top")

map_cov <- ggpubr::ggarrange(map_temp, map_hab, ncol = 2)
map_cov
# ggsave(map_cov, filename = "map_cov.pdf", width = 9, height = 8)



















