
rm(list = ls())

library(pROC)

drg <- function(x, mean, sd, supp) {
  num <- pnorm(x+0.5, mean = mean, sd = sd) - pnorm(x-0.5, mean = mean, sd = sd)
  den <- 1 - pnorm(min(supp)-0.5, mean = mean, sd = sd)
  pmax(1e-250, (num / den))
}
k_star_bin <- function(y, y_n, delta) {
  if (y == y_n) {
    return(1-delta)
  } else {
    return(delta)
  }
}
mad_class2 <- function(y, X, p_0, hyper, support_x) {
  sigma_1 <- hyper[1]
  sigma_2 <- hyper[2]
  delta <- hyper[3]
  N <- length(y)
  p_n <- p_0
  for (n in 1:N) {
    x_obs <- X[n,]
    y_obs <- y[n]
    pn_obs <- p_n[x_obs[1]+1, x_obs[2]+1, y_obs+1]
    w_n <- (2 - 1/n)*((n+1)^(-1))
    ker <- array(dim = dim(p_0))
    for (y_val in 0:1) {
      for (x_2 in support_x) {
        for (x_1 in support_x) {
          log_k0_y_yn <- log(drg(x_1, x_obs[1], sigma_1, support_x)) + 
            log(drg(x_2, x_obs[2], sigma_2, support_x)) + log(k_star_bin(y_val, y_obs, delta))
          log_k0_yn_y <- log(drg(x_obs[1], x_1, sigma_1, support_x)) +
            log(drg(x_obs[2], x_2, sigma_2, support_x)) + log(k_star_bin(y_obs, y_val, delta))
          pn <- p_n[x_1+1, x_2+1, y_val+1]
          
          log_gamma <- pmin(0, log(pn) + log_k0_yn_y - log(pn_obs) - log_k0_y_yn)
          ker[x_1+1, x_2+1, y_val+1] <- exp(log_gamma + log_k0_y_yn)
        }
      }
    }
    p_n <- (1-w_n)*p_n + w_n*ker
    p_n[x_obs[1]+1, x_obs[2]+1, y_obs+1] <- p_n[x_obs[1]+1, x_obs[2]+1, y_obs+1] + w_n*(1-sum(ker))
  }
  p_n
}
mad_class2_preq <- function(y, X, p_0, hyper, support_x) {
  sigma_1 <- hyper[1]
  sigma_2 <- hyper[2]
  delta <- hyper[3]
  N <- length(y)
  p_n <- p_0
  out <- 0
  for (n in 1:N) {
    x_obs <- X[n,]
    y_obs <- y[n]
    pn_obs <- p_n[x_obs[1]+1, x_obs[2]+1, y_obs+1]
    out <- out + log(pn_obs)
    w_n <- (2 - 1/n)*((n+1)^(-1))
    ker <- array(dim = dim(p_0))
    for (y_val in 0:1) {
      for (x_2 in support_x) {
        for (x_1 in support_x) {
          log_k0_y_yn <- log(drg(x_1, x_obs[1], sigma_1, support_x)) + 
            log(drg(x_2, x_obs[2], sigma_2, support_x)) + log(k_star_bin(y_val, y_obs, delta))
          log_k0_yn_y <- log(drg(x_obs[1], x_1, sigma_1, support_x)) +
            log(drg(x_obs[2], x_2, sigma_2, support_x)) + log(k_star_bin(y_obs, y_val, delta))
          pn <- p_n[x_1+1, x_2+1, y_val+1]
          
          log_gamma <- pmin(0, log(pn) + log_k0_yn_y - log(pn_obs) - log_k0_y_yn)
          ker[x_1+1, x_2+1, y_val+1] <- exp(log_gamma + log_k0_y_yn)
        }
      }
    }
    p_n <- (1-w_n)*p_n + w_n*ker
    p_n[x_obs[1]+1, x_obs[2]+1, y_obs+1] <- p_n[x_obs[1]+1, x_obs[2]+1, y_obs+1] + w_n*(1-sum(ker))
  }
  out
}

copula_class2 <- function(y, X, p_0, rho, support_x) {
  rho_1 <- rho[1]
  rho_2 <- rho[2]
  rho_3 <- rho[3]
  N <- length(y)
  p_n <- p_0
  for (n in 1:N) {
    x_obs <- X[n,]
    y_obs <- y[n]
    w_n <- (2 - 1/n)*((n+1)^(-1))
    p_n_temp <- array(dim = dim(p_0))
    for (y_val in 0:1) {
      for (x_2 in support_x) {
        for (x_1 in support_x) {
          pn_old <- p_n[x_1+1, x_2+1, y_val+1]
          pn_x1 <- sum(p_n[x_1+1,,])
          pn_x2_x1 <- sum(p_n[x_1+1, x_2+1,]) / pn_x1
          pn_x <- sum(p_n[x_1+1, x_2+1,])
          pn_y_x <- p_n[x_1+1, x_2+1, y_val+1] / pn_x
          pn_x2y_x1 <- p_n[x_1+1, x_2+1, y_val+1] / pn_x1
          p_n_temp[x_1+1, x_2+1, y_val+1] <- pn_old*(1 - w_n + w_n*(1-rho_1)*(1-rho_2)*(1-rho_3))
          if (x_1 == x_obs[1]) p_n_temp[x_1+1, x_2+1, y_val+1] <- p_n_temp[x_1+1, x_2+1, y_val+1] + w_n*rho_1*(1-rho_2)*(1-rho_3)*pn_old/pn_x1
          if (x_2 == x_obs[2]) p_n_temp[x_1+1, x_2+1, y_val+1] <- p_n_temp[x_1+1, x_2+1, y_val+1] + w_n*(1-rho_1)*rho_2*(1-rho_3)*pn_old/pn_x2_x1
          if (y_val == y_obs) p_n_temp[x_1+1, x_2+1, y_val+1] <- p_n_temp[x_1+1, x_2+1, y_val+1] + w_n*(1-rho_1)*(1-rho_2)*rho_3*pn_x
          if (x_1 == x_obs[1] & x_2 == x_obs[2]) p_n_temp[x_1+1, x_2+1, y_val+1] <- p_n_temp[x_1+1, x_2+1, y_val+1] + w_n*rho_1*rho_2*(1-rho_3)*pn_y_x
          if (x_1 == x_obs[1] & y_val == y_obs) p_n_temp[x_1+1, x_2+1, y_val+1] <- p_n_temp[x_1+1, x_2+1, y_val+1] + w_n*rho_1*(1-rho_2)*rho_3*pn_old/(pn_x1*pn_y_x)
          if (x_2 == x_obs[2] & y_val == y_obs) p_n_temp[x_1+1, x_2+1, y_val+1] <- p_n_temp[x_1+1, x_2+1, y_val+1] + w_n*(1-rho_1)*rho_2*rho_3*pn_old/pn_x2y_x1
        }
      }
    }
    p_n_temp[x_obs[1]+1, x_obs[2]+1, y_obs+1] <- p_n_temp[x_obs[1]+1, x_obs[2]+1, y_obs+1] + w_n*rho_1*rho_2*rho_3
    p_n <- p_n_temp
  }
  p_n
} 
copula_class2_preq <- function(y, X, p_0, rho, support_x) {
  rho_1 <- rho[1]
  rho_2 <- rho[2]
  rho_3 <- rho[3]
  N <- length(y)
  p_n <- p_0
  out <- 0
  for (n in 1:N) {
    x_obs <- X[n,]
    y_obs <- y[n]
    out <- out + log(p_n[x_obs[1]+1, x_obs[2]+1, y_obs+1])
    w_n <- (2 - 1/n)*((n+1)^(-1))
    p_n_temp <- array(dim = dim(p_0))
    for (y_val in 0:1) {
      for (x_2 in support_x) {
        for (x_1 in support_x) {
          pn_old <- p_n[x_1+1, x_2+1, y_val+1]
          pn_x1 <- sum(p_n[x_1+1,,])
          pn_x2_x1 <- sum(p_n[x_1+1, x_2+1,]) / pn_x1
          pn_x <- sum(p_n[x_1+1, x_2+1,])
          pn_y_x <- p_n[x_1+1, x_2+1, y_val+1] / pn_x
          pn_x2y_x1 <- p_n[x_1+1, x_2+1, y_val+1] / pn_x1
          p_n_temp[x_1+1, x_2+1, y_val+1] <- pn_old*(1 - w_n + w_n*(1-rho_1)*(1-rho_2)*(1-rho_3))
          if (x_1 == x_obs[1]) p_n_temp[x_1+1, x_2+1, y_val+1] <- p_n_temp[x_1+1, x_2+1, y_val+1] + w_n*rho_1*(1-rho_2)*(1-rho_3)*pn_old/pn_x1
          if (x_2 == x_obs[2]) p_n_temp[x_1+1, x_2+1, y_val+1] <- p_n_temp[x_1+1, x_2+1, y_val+1] + w_n*(1-rho_1)*rho_2*(1-rho_3)*pn_old/pn_x2_x1
          if (y_val == y_obs) p_n_temp[x_1+1, x_2+1, y_val+1] <- p_n_temp[x_1+1, x_2+1, y_val+1] + w_n*(1-rho_1)*(1-rho_2)*rho_3*pn_x
          if (x_1 == x_obs[1] & x_2 == x_obs[2]) p_n_temp[x_1+1, x_2+1, y_val+1] <- p_n_temp[x_1+1, x_2+1, y_val+1] + w_n*rho_1*rho_2*(1-rho_3)*pn_y_x
          if (x_1 == x_obs[1] & y_val == y_obs) p_n_temp[x_1+1, x_2+1, y_val+1] <- p_n_temp[x_1+1, x_2+1, y_val+1] + w_n*rho_1*(1-rho_2)*rho_3*pn_old/(pn_x1*pn_y_x)
          if (x_2 == x_obs[2] & y_val == y_obs) p_n_temp[x_1+1, x_2+1, y_val+1] <- p_n_temp[x_1+1, x_2+1, y_val+1] + w_n*(1-rho_1)*rho_2*rho_3*pn_old/pn_x2y_x1
        }
      }
    }
    p_n_temp[x_obs[1]+1, x_obs[2]+1, y_obs+1] <- p_n_temp[x_obs[1]+1, x_obs[2]+1, y_obs+1] + w_n*rho_1*rho_2*rho_3
    p_n <- p_n_temp
  }
  out
} 
copula_class2_flip <- function(y, X, p_0, rho, support_x) {
  rho_1 <- rho[1]
  rho_2 <- rho[2]
  rho_3 <- rho[3]
  N <- length(y)
  p_n <- p_0
  for (n in 1:N) {
    x_obs <- X[n,]
    y_obs <- y[n]
    w_n <- (2 - 1/n)*((n+1)^(-1))
    p_n_temp <- array(dim = dim(p_0))
    for (y_val in 0:1) {
      for (x_2 in support_x) {
        for (x_1 in support_x) {
          pn_old <- p_n[x_1+1, x_2+1, y_val+1]
          pn_x2 <- sum(p_n[,x_2+1,])
          pn_x1_x2 <- sum(p_n[x_1+1, x_2+1,]) / pn_x2
          pn_x <- sum(p_n[x_1+1, x_2+1,])
          pn_y_x <- p_n[x_1+1, x_2+1, y_val+1] / pn_x
          pn_x1y_x2 <-  p_n[x_1+1, x_2+1, y_val+1] / pn_x2
          p_n_temp[x_1+1, x_2+1, y_val+1] <- pn_old*(1 - w_n + w_n*(1-rho_1)*(1-rho_2)*(1-rho_3))
          if (x_2 == x_obs[2]) p_n_temp[x_1+1, x_2+1, y_val+1] <- p_n_temp[x_1+1, x_2+1, y_val+1] + w_n*(1-rho_1)*rho_2*(1-rho_3)*pn_old/pn_x2
          if (x_1 == x_obs[1]) p_n_temp[x_1+1, x_2+1, y_val+1] <- p_n_temp[x_1+1, x_2+1, y_val+1] + w_n*rho_1*(1-rho_2)*(1-rho_3)*pn_old/pn_x1_x2
          if (y_val == y_obs) p_n_temp[x_1+1, x_2+1, y_val+1] <- p_n_temp[x_1+1, x_2+1, y_val+1] + w_n*(1-rho_1)*(1-rho_2)*rho_3*pn_x
          if (x_1 == x_obs[1] & x_2 == x_obs[2]) p_n_temp[x_1+1, x_2+1, y_val+1] <- p_n_temp[x_1+1, x_2+1, y_val+1] + w_n*rho_1*rho_2*(1-rho_3)*pn_y_x
          if (x_2 == x_obs[2] & y_val == y_obs) p_n_temp[x_1+1, x_2+1, y_val+1] <- p_n_temp[x_1+1, x_2+1, y_val+1] + w_n*(1-rho_1)*rho_2*rho_3*pn_old/(pn_x2*pn_y_x)
          if (x_1 == x_obs[1] & y_val == y_obs) p_n_temp[x_1+1, x_2+1, y_val+1] <- p_n_temp[x_1+1, x_2+1, y_val+1] + w_n*rho_1*(1-rho_2)*rho_3*pn_old/pn_x1y_x2
        }
      }
    }
    p_n_temp[x_obs[1]+1, x_obs[2]+1, y_obs+1] <- p_n_temp[x_obs[1]+1, x_obs[2]+1, y_obs+1] + w_n*rho_1*rho_2*rho_3
    p_n <- p_n_temp
  }
  p_n
}
copula_class2_preq_flip <- function(y, X, p_0, rho, support_x) {
  rho_1 <- rho[1]
  rho_2 <- rho[2]
  rho_3 <- rho[3]
  N <- length(y)
  p_n <- p_0
  out <- 0
  for (n in 1:N) {
    x_obs <- X[n,]
    y_obs <- y[n]
    out <- out + log(p_n[x_obs[1]+1, x_obs[2]+1, y_obs+1])
    w_n <- (2 - 1/n)*((n+1)^(-1))
    p_n_temp <- array(dim = dim(p_0))
    for (y_val in 0:1) {
      for (x_2 in support_x) {
        for (x_1 in support_x) {
          pn_old <- p_n[x_1+1, x_2+1, y_val+1]
          pn_x2 <- sum(p_n[,x_2+1,])
          pn_x1_x2 <- sum(p_n[x_1+1, x_2+1,]) / pn_x2
          pn_x <- sum(p_n[x_1+1, x_2+1,])
          pn_y_x <- p_n[x_1+1, x_2+1, y_val+1] / pn_x
          pn_x1y_x2 <-  p_n[x_1+1, x_2+1, y_val+1] / pn_x2
          p_n_temp[x_1+1, x_2+1, y_val+1] <- pn_old*(1 - w_n + w_n*(1-rho_1)*(1-rho_2)*(1-rho_3))
          if (x_2 == x_obs[2]) p_n_temp[x_1+1, x_2+1, y_val+1] <- p_n_temp[x_1+1, x_2+1, y_val+1] + w_n*(1-rho_1)*rho_2*(1-rho_3)*pn_old/pn_x2
          if (x_1 == x_obs[1]) p_n_temp[x_1+1, x_2+1, y_val+1] <- p_n_temp[x_1+1, x_2+1, y_val+1] + w_n*rho_1*(1-rho_2)*(1-rho_3)*pn_old/pn_x1_x2
          if (y_val == y_obs) p_n_temp[x_1+1, x_2+1, y_val+1] <- p_n_temp[x_1+1, x_2+1, y_val+1] + w_n*(1-rho_1)*(1-rho_2)*rho_3*pn_x
          if (x_1 == x_obs[1] & x_2 == x_obs[2]) p_n_temp[x_1+1, x_2+1, y_val+1] <- p_n_temp[x_1+1, x_2+1, y_val+1] + w_n*rho_1*rho_2*(1-rho_3)*pn_y_x
          if (x_2 == x_obs[2] & y_val == y_obs) p_n_temp[x_1+1, x_2+1, y_val+1] <- p_n_temp[x_1+1, x_2+1, y_val+1] + w_n*(1-rho_1)*rho_2*rho_3*pn_old/(pn_x2*pn_y_x)
          if (x_1 == x_obs[1] & y_val == y_obs) p_n_temp[x_1+1, x_2+1, y_val+1] <- p_n_temp[x_1+1, x_2+1, y_val+1] + w_n*rho_1*(1-rho_2)*rho_3*pn_old/pn_x1y_x2
        }
      }
    }
    p_n_temp[x_obs[1]+1, x_obs[2]+1, y_obs+1] <- p_n_temp[x_obs[1]+1, x_obs[2]+1, y_obs+1] + w_n*rho_1*rho_2*rho_3
    p_n <- p_n_temp
  }
  out
}

N <- 100
nSim <- 50
prob_w <- 0.7
L1 <- 3; L2 <- 12; L3 <- 9
B0 <- 6; B1 <- -2.1; B2 <- 0
supp <- 0:40

roc_A <- roc_B <- roc_mad <- list()
for (s in 1:nSim) {
  print(s)
  
  # Data simulation
  set.seed(s)
  xx <- matrix(nrow = N, ncol = 2)
  w <- rbinom(N, 1, prob = prob_w)
  for (i in 1:N) {
    if (w[i] == 1) xx[i,1] <- rpois(1, lambda = L1)
    if (w[i] == 0) xx[i,1] <- rpois(1, lambda = L2)
  }
  xx[,2] <- rpois(N, lambda = L3)
  my <- plogis(B0 + B1*xx[,1] + B2*xx[,2])
  yy <- rbinom(N, 1, my)
  
  # Out-of-sample 
  N_new <- 1e5
  xx_new <- matrix(nrow = N_new, ncol = 2)
  set.seed(2025)
  w_new <- rbinom(N_new, 1, prob = prob_w)
  for (i in 1:N_new) {
    if (w_new[i] == 1) xx_new[i,1] <- rpois(1, lambda = L1)
    if (w_new[i] == 0) xx_new[i,1] <- rpois(1, lambda = L2)
  }
  xx_new[,2] <- rpois(N_new, lambda = L3)
  my_new <- plogis(B0 + B1*xx_new[,1] + B2*xx_new[,2])
  yy_new <- rbinom(N_new, 1, my_new)
  
  p00 <- array(1, dim = c(length(supp), length(supp), 2)); p00 <- p00 / sum(p00)
  
  # Ordering 1
  opt_A <- nlminb(start = rep(0.5, 3), lower = rep(1e-10, 3), upper = rep(1-1e-10, 3),
                objective = function(r) -copula_class2_preq(yy, xx, p00, r, supp), control = list(rel.tol = 1e-4))
  ppc_A <- copula_class2(yy, xx, p00, opt_A$par, supp)
  
  # Ordering 2
  opt_B <- nlminb(start = rep(0.5, 3), lower = rep(1e-10, 3), upper = rep(1-1e-10, 3),
                  objective = function(r) -copula_class2_preq_flip(yy, xx, p00, r, supp), control = list(rel.tol = 1e-4))
  ppc_B <- copula_class2_flip(yy, xx, p00, opt_B$par, supp)
  
  # MAD
  opt_mad <- nlminb(start = c(2, 2, 0.5), lower = rep(1e-10, 3), upper = c(Inf, Inf, 1-1e-10),
                    objective = function(r) -mad_class2_preq(yy, xx, p00, r, supp), control = list(rel.tol = 1e-4))
  ppm <- mad_class2(yy, xx, p00, opt_mad$par, supp)
  
  # Predictions
  pred_mad <- pred_A <- pred_B <- rep(NA, N_new)
  for (i in 1:N_new) {
    pred_mad[i] <- ppm[xx_new[i,1]+1, xx_new[i,2]+1, 2] / sum(ppm[xx_new[i,1]+1, xx_new[i,2]+1,])
    pred_A[i] <- ppc_A[xx_new[i,1]+1, xx_new[i,2]+1, 2] / sum(ppc_A[xx_new[i,1]+1, xx_new[i,2]+1,])
    pred_B[i] <- ppc_B[xx_new[i,1]+1, xx_new[i,2]+1, 2] / sum(ppc_B[xx_new[i,1]+1, xx_new[i,2]+1,])
  }
  roc_mad[[s]] <- roc(yy_new, pred_mad)
  roc_A[[s]] <- roc(yy_new, pred_A)
  roc_B[[s]] <- roc(yy_new, pred_B)
  
  name_file <- paste0("SIM_copula_", N, ".RData")
  save(roc_A, roc_B, roc_mad, file = file(name_file, blocking = TRUE))
}

nSim <- 50

load("SIM_copula_50.RData")
auc_mad_50 <- auc_A_50 <- auc_B_50 <- rep(NA, length(roc_mad))
for (i in 1:length(roc_mad)) {
  auc_mad_50[i] <- auc(roc_mad[[i]])
  auc_A_50[i] <- auc(roc_A[[i]])
  auc_B_50[i] <- auc(roc_B[[i]])
}
(AB_50 <- sum(auc_A_50 > auc_B_50) / nSim)
sum(auc_B_50 < auc_mad_50)
boxplot(auc_A_50, auc_B_50, auc_mad_50)
cbind(auc_A_50, auc_B_50, auc_mad_50)

load("SIM_copula_50_05.RData")
auc_mad_50_05 <- auc_A_50_05 <- auc_B_50_05 <- rep(NA, length(roc_mad))
for (i in 1:length(roc_mad)) {
  auc_mad_50_05[i] <- auc(roc_mad[[i]])
  auc_A_50_05[i] <- auc(roc_A[[i]])
  auc_B_50_05[i] <- auc(roc_B[[i]])
}
(AB_50_05 <- sum(auc_A_50_05 > auc_B_50_05) / nSim)
sum(auc_B_50_05 < auc_mad_50_05)
boxplot(auc_A_50_05, auc_B_50_05, auc_mad_50_05)
cbind(auc_A_50_05, auc_B_50_05, auc_mad_50_05)

load("SIM_copula_100.RData")
auc_mad_100 <- auc_A_100 <- auc_B_100 <- rep(NA, length(roc_mad))
for (i in 1:length(roc_mad)) {
  auc_mad_100[i] <- auc(roc_mad[[i]])
  auc_A_100[i] <- auc(roc_A[[i]])
  auc_B_100[i] <- auc(roc_B[[i]])
}
(AB_100 <- sum(auc_A_100 > auc_B_100) / nSim)
sum(auc_B_100 < auc_mad_100)
boxplot(auc_A_100, auc_B_100, auc_mad_100)
cbind(auc_A_100, auc_B_100, auc_mad_100)

load("SIM_copula_100_05.RData")
auc_mad_100_05 <- auc_A_100_05 <- auc_B_100_05 <- rep(NA, length(roc_mad))
for (i in 1:length(roc_mad)) {
  auc_mad_100_05[i] <- auc(roc_mad[[i]])
  auc_A_100_05[i] <- auc(roc_A[[i]])
  auc_B_100_05[i] <- auc(roc_B[[i]])
}
(AB_100_05 <- sum(auc_A_100_05 > auc_B_100_05) / nSim)
sum(auc_B_100_05 < auc_mad_100_05)
boxplot(auc_A_100_05, auc_B_100_05, auc_mad_100_05)
cbind(auc_A_100_05, auc_B_100_05, auc_mad_100_05)

load("SIM_copula_150.RData")
auc_mad_150 <- auc_A_150 <- auc_B_150 <- rep(NA, length(roc_mad))
for (i in 1:length(roc_mad)) {
  auc_mad_150[i] <- auc(roc_mad[[i]])
  auc_A_150[i] <- auc(roc_A[[i]])
  auc_B_150[i] <- auc(roc_B[[i]])
}
(AB_150 <- sum(auc_A_150 > auc_B_150) / nSim)
sum(auc_B_150 < auc_mad_150)
boxplot(auc_A_150, auc_B_150, auc_mad_150)
cbind(auc_A_150, auc_B_150, auc_mad_150)

load("SIM_copula_150_05.RData")
auc_mad_150_05 <- auc_A_150_05 <- auc_B_150_05 <- rep(NA, length(roc_mad))
for (i in 1:length(roc_mad)) {
  auc_mad_150_05[i] <- auc(roc_mad[[i]])
  auc_A_150_05[i] <- auc(roc_A[[i]])
  auc_B_150_05[i] <- auc(roc_B[[i]])
}
(AB_150_05 <- sum(auc_A_150_05 > auc_B_150_05) / nSim)
sum(auc_B_150_05 < auc_mad_150_05)
boxplot(auc_A_150_05, auc_B_150_05, auc_mad_150_05)
cbind(auc_A_150_05, auc_B_150_05, auc_mad_150_05)


tab <- matrix(c(AB_50, AB_50_05, AB_100, AB_100_05, AB_150, AB_150_05), nrow = 2, ncol = 3)
colnames(tab) <- c("n=50", "n=100", "n=150")
rownames(tab) <- c("B=0", "B=0.5")
tab
  

# Figure S3 of the Supplement
library(ggplot2)
df <- data.frame(auc = c(auc_A_50, auc_B_50, auc_mad_50,
                         auc_A_50_05, auc_B_50_05, auc_mad_50_05,
                         auc_A_100, auc_B_100, auc_mad_100,
                         auc_A_100_05, auc_B_100_05, auc_mad_100_05,
                         auc_A_150, auc_B_150, auc_mad_150,
                         auc_A_150_05, auc_B_150_05, auc_mad_150_05),
                 nn = rep(c(50, 100, 150), each = nSim*6),
                 beta_val = rep(c(0, 0.5), each = nSim*3, 3),
                 type = rep(c("COP-A", "COP-B", "MAD"), each = nSim, 6))
df$nn <- as.factor(df$nn)
levels(df$nn) <- c("n = 50", "n = 100", "n = 150")
pl_cc <- ggplot(df) +
  facet_grid(rows = vars(beta_val), cols = vars(nn),
             labeller = label_bquote(beta[.(2)] == .(beta_val))) +
  geom_boxplot(aes(x = type, y = auc), fill = "gray90", outlier.size = .1) +
  ylim(c(0.72, 1)) + xlab("") + ylab("AUC") +
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
pl_cc
ggsave(pl_cc, filename = "pl_cc.pdf", width = 9, height = 6)

