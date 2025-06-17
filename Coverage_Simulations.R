rm(list = ls())
setwd("/Users/davideagnoletto/Dropbox/sim_martingales/code")

load("/regression/coverage_regr.RData")
cv_regr <- cv_df
cv_regr$type <- rep(c("Regression"), NROW(cv_regr))
load("/classification/coverage_class.RData")
cv_class <- cv_df
cv_class$type <- rep(c("Classification"), NROW(cv_class))

library(ggplot2)
library(ggthemes)
library(ggpubr)

pl_regr <- ggplot(cv_regr) +
  facet_wrap(~ type) +
  geom_boxplot(aes(x = n, y = cv, fill = model), outlier.size = .1) + 
  scale_fill_tableau( name = "Model", palette = "Color Blind",
                      type = "regular", direction = 1,
                      labels = c("1" = "GLM",
                                 "2" = "BART",
                                 "3" = "RF",
                                 "4" = "DP",
                                 "5" = "MAD-1",
                                 "6" = "MAD-2/3",
                                 "7" = "MAD-DPM",
                                 "8" = "MAD-ada")) +
  guides(fill = guide_legend(nrow = 1)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red2") +
  ylim(c(0,1)) + xlab("Sample size") + ylab("Coverage") +
  theme_light() +
  theme(legend.text.align = 0,
        strip.text = element_text(size = 11, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 11),
        legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top")
pl_regr

pl_class <- ggplot(cv_class) +
  facet_wrap(~ type) +
  geom_boxplot(aes(x = n, y = cv, fill = model), outlier.size = .1) + 
  scale_fill_tableau( name = "Model", palette = "Color Blind",
                      type = "regular", direction = 1,
                      labels = c("1" = "GLM",
                                 "2" = "BART",
                                 "3" = "RF",
                                 "4" = "DP",
                                 "5" = "MAD-1",
                                 "6" = "MAD-2/3",
                                 "7" = "MAD-DPM",
                                 "8" = "MAD-ada")) +
  guides(fill = guide_legend(nrow = 1)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red2") +
  ylim(c(0,1)) + xlab("Sample size") + ylab("Coverage") +
  theme_light() +
  theme(legend.text.align = 0,
        strip.text = element_text(size = 11, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 11),
        legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top")
pl_class

pl_cov <- ggarrange(pl_regr, pl_class, common.legend = TRUE)
pl_cov
# ggsave(pl_cov, filename = "pl_cov_sim.pdf", width = 11, height = 6)



