
rm(list = ls())
setwd("/Users/davideagnoletto/Dropbox/sim_martingales/code/application")

library(Rcpp)
library(RcppArmadillo)
library(SparseArray)
library(ggplot2)
library(ggcorrplot)

library(Hmsc)
library(gjam)
library(gllvm)

# DATA
data_2009 <- read.csv("corvus_2009.csv")
y <- data_2009[,1:4]
temp <- as.numeric(data_2009$Temperature)
hab <- as.numeric(data_2009$Habitat)
sites_2009 <- data_2009$Site

# TEST SET (Year 2010)
data_2010 <- read.csv("corvus_2010.csv")
y_test <- data_2010[,1:4]
temp_test <- as.numeric(data_2010$Temperature)
hab_test <- as.numeric(data_2010$Habitat)
sites_2010 <- data_2010$Site

# MAD adaptive
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
                                                      10, 5, 1),
              control = list(rel.tol = 1e-3))
opt
sigma_opt <- opt$par[1:4]
delta_opt <- opt$par[5:6]
system.time(pp <- pred_4_FB_perm(as.matrix(y), temp, hab, p_0, p_0_az, alpha, lambda, N_star, sigma_opt, 
                                 delta_opt, supp, supp_temp, supp_hab, 10, 5, 12345))
save(pp, file = file("mad_object_a.RData", blocking = TRUE))
mad_exp <- predict_mad(pp, temp_test, hab_test)
mean((as.matrix(y_test) - mad_exp)^2)


# GJAM (Generalized Joint Attribute Models)
Xmm <- as.data.frame(model.matrix(~ as.factor(temp) + as.factor(hab)))
colnames(Xmm) <- c("intercept", "temp1", "temp2", "hab2", "hab3", "hab4", "hab5")
form <- as.formula(~ temp1 + temp2 + hab2 + hab3 + hab4 + hab5)
ml <- list(ng = 2000, burnin = 1000, typeNames = 'DA')
set.seed(12345)
gg <- gjam(form, xdata = Xmm, ydata = y, modelList = ml)

Xmm_test <- as.data.frame(model.matrix(~ as.factor(temp_test) + as.factor(hab_test)))
Xmm_test <- cbind(Xmm_test, rep(0, NROW(y_test)))
colnames(Xmm_test) <- c("intercept", "temp1", "temp2", "hab2", "hab3", "hab4", "hab5")
newdata <- list(xdata = Xmm_test, nsim = 200)
pj <- gjamPredict(gg, newdata = newdata)
gjam_exp <- matrix(colMeans(pj$ychains), ncol = D)
mean((as.matrix(y_test) - gjam_exp)^2)


# HMSC (Hierarchical Modeling of Species Communities)
studyDesign <- data.frame(sample = as.factor(sites_2009))
rL <- HmscRandomLevel(units = studyDesign$sample)
dd <- data.frame(temp = as.factor(temp), hab = as.factor(hab))
mm <- Hmsc(Y = y, XData = dd, XFormula = ~ temp + hab, distr = "poisson",
           studyDesign = studyDesign, ranLevels = list(sample = rL))
set.seed(12345)
mm <-  sampleMcmc(mm, thin = 1, samples = 3000, transient = 1000,
                  nChains = 4, nParallel = 4, verbose = TRUE)

studyDesign_test <- data.frame(sample = as.factor(sites_2010))
rL_test <- HmscRandomLevel(units = studyDesign_test$sample)
dd_test <- data.frame(temp = as.factor(temp_test), hab = as.factor(hab_test))
ph <- predict(mm, XData = dd_test, XFormula = ~ temp + hab, 
              studyDesign = studyDesign_test, ranLevels = list(sample = rL_test))
ph_mat <- matrix(unlist(ph), nrow = length(sites_2010)*D)
hmsc_exp <- matrix(rowMeans(ph_mat), nrow = length(sites_2010), ncol = D)


# GLLVM (Generalized Linear Latent Variable Models)
vv_nb <- gllvm(y = y, family = "negative.binomial", num.lv = 1, seed = 12345)
vv_zinb <- gllvm(y = y, family = "ZINB", num.lv = 1, seed = 12345)
levels(dd_test$hab) <- levels(dd$hab)
gllvm_nb_exp <- predict(vv_nb, newX = dd_test, level = 0, type = "response")
gllvm_zinb_exp <- predict(vv_zinb, newX = dd_test, level = 0, type = "response")


# COMPARISON in terms of predictive accuracy on the TEST SET
mean((as.matrix(y_test) - mad_exp)^2)
mean((as.matrix(y_test) - gjam_exp)^2)
mean((as.matrix(y_test) - hmsc_exp)^2)
mean((as.matrix(y_test) - gllvm_nb_exp)^2)
mean((as.matrix(y_test) - gllvm_zinb_exp)^2)

mean(abs(as.matrix(y_test) - mad_exp))
mean(abs(as.matrix(y_test) - gjam_exp))
mean(abs(as.matrix(y_test) - hmsc_exp))
mean(abs(as.matrix(y_test) - gllvm_nb_exp))
mean(abs(as.matrix(y_test) - gllvm_zinb_exp))



# FIGURE: Correlation matrices
cor_names <- c("CC", "EM", "CR", "EJ")
cor_emp <- cor(y)
colnames(cor_emp) <- rownames(cor_emp) <- cor_names

cov_mad <- cov_4_FB(pp$p_n, D, supp)
cor_mad <- cor_4_FB(cov_mad)
colnames(cor_mad) <- rownames(cor_mad) <- cor_names

OmegaCor <- computeAssociations(mm)
cor_hmsc <- OmegaCor[[1]][[1]]
colnames(cor_hmsc) <- rownames(cor_hmsc) <- cor_names

cor_gjam <- gg$parameters$corMu
colnames(cor_gjam) <- rownames(cor_gjam) <- cor_names

cor_nb <- getResidualCor(vv_nb)
colnames(cor_nb) <- rownames(cor_nb) <- cor_names
cor_zinb <- getResidualCor(vv_zinb)
colnames(cor_zinb) <- rownames(cor_zinb) <- cor_names

c_emp <- ggcorrplot(cor_emp, colors = c("red", "white", "blue"), tl.srt = 0, tl.cex = 11, 
                    tl.col = "black", show.legend = FALSE, title = "Empirical")
c_mad <- ggcorrplot(cor_mad, colors = c("red", "white", "blue"), tl.srt = 0, tl.cex = 11, 
                    tl.col = "black", show.legend = FALSE, title = "MAD")
c_hmsc <- ggcorrplot(cor_hmsc, colors = c("red", "white", "blue"), tl.srt = 0, tl.cex = 11, 
                     tl.col = "black", show.legend = FALSE, title = "HMSC")
c_gjam <- ggcorrplot(cor_gjam, colors = c("red", "white", "blue"), tl.srt = 0, tl.cex = 11, 
                     tl.col = "black", show.legend = FALSE, title = "GJAM")
c_nb <- ggcorrplot(cor_nb, colors = c("red", "white", "blue"), tl.srt = 0, tl.cex = 11, 
                   tl.col = "black", show.legend = FALSE, title = "GLLVM - NB")
c_zinb <- ggcorrplot(cor_zinb, colors = c("red", "white", "blue"), tl.srt = 0, tl.cex = 11, 
                     tl.col = "black", show.legend = FALSE, title = "GLLVM - ZINB")
cplot <- gridExtra::grid.arrange(c_emp, c_mad, c_hmsc, c_gjam, c_nb, c_zinb, nrow = 2)
# ggsave(file = "cor_plot.pdf", cplot, height = 6, width = 9)





