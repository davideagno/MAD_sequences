
k0_int <- function(x, marg, support) {
  dx <- pmin(max(support), ceiling(x + marg))
  sx <- pmax(min(support), floor(x - marg))
  sx:dx
}
rel_freq <- function(x, support) {
  D <- length(support)
  out <- rep(0, D)
  N <- length(x)
  out[as.numeric(dimnames(table(x))$x) + 1] <- table(x) / N
  out
}
rho <- function(n, lambda, N_star) {
  lambda + (1 - lambda)*exp(-(1/N_star)*n)
}

pred_4_FB <- function(y, temp, hab, p_0, p_0_az, alpha, lambda, N_star, sigma, delta, support_y, support_temp, support_hab) {
  
  d <- 4
  require(crch)
  left <- min(support_y) - 0.5
  interval <- rep(0, d)
  for (k in 1:d) {
    extremes <- crch::qtnorm(c(0.025, 0.975), mean = 50, sd = sigma[k], left = left)
    interval[k] <- max(0, (extremes[2] - extremes[1]) / 2)
  }
  
  lsy <- length(support_y)
  lst <- length(support_temp)
  lsh <- length(support_hab)
  delta_temp <- delta[1]
  delta_hab <- delta[2]
  
  p_0[y[1,1]+1, y[1,2]+1, y[1,3]+1, y[1,4]+1, temp[1]+1, hab[1]] <- p_0_az
  out_z <- p_0_az
  preqLL <- 0
  
  for (n in 1:NROW(y)) {
    y_obs <- c(y[n,])
    temp_obs <- temp[n]
    hab_obs <- hab[n]
    idx_obs <- y_obs + 1
    
    # Recover p_{n-1}(y_n, x_n)
    p_n_obs <- p_0[idx_obs[1], idx_obs[2], idx_obs[3], idx_obs[4], temp_obs+1, hab_obs]
    if (p_n_obs == 0) p_n_obs <- out_z
    
    # Update prequential log-likelihood
    preqLL <- preqLL + log(p_n_obs) # - log(sum(p_0[,,,, temp_obs+1, hab_obs]))
    
    # Ball around y_obs
    ball <- list()
    for (k in 1:d) {
      ball[[k]] <- k0_int(y_obs[k], interval[k], support_y)
    }
    ball[[d+1]] <- support_temp
    ball[[d+2]] <- support_hab
    
    # Index conversion
    idx_tr <- c()
    # idx_tr <- ball2idx_5(ball, support_hab, support_temp, lsy)
    idx_tr <- mat2idx(as.matrix(expand.grid(ball)), d, lsy, lst, lsh)
    
    # Compute \gamma(y,x, y_n,x_n) * k_0(y,x | y_n,x_n)
    p_0[idx_tr][p_0[idx_tr] == 0] <- out_z
    ker <- mhk_4_FB(p_0[idx_tr], p_n_obs, y_obs, temp_obs, hab_obs, ball,
                    support_temp, support_hab, sigma, delta_temp, delta_hab)
    
    # Update p_{n-1}
    rho_n <- rho(n, lambda, N_star)
    w_n <- 1 / ((alpha + n)^rho_n)
    p_0 <- (1 - w_n) * p_0
    p_0[idx_tr] <- p_0[idx_tr] + w_n*ker
    
    # Add \sum_z {\gamma(z,y_n)*k_0(z | y_n)} in y_n
    p_0[idx_obs[1], idx_obs[2], idx_obs[3], idx_obs[4], temp_obs+1, hab_obs] <- p_0[idx_obs[1], idx_obs[2], idx_obs[3], idx_obs[4], temp_obs+1, hab_obs] + w_n*(1 - sum(ker))
    
    # Update other values of p_n
    out_z <- (1 - w_n) * out_z
  }
  
  return(list(p_n = p_0, p_az = out_z, preqLL = preqLL))
}
pred_4_FB_perm <- function(y, temp, hab, p_0, p_0_az, alpha, lambda, N_star, sigma, delta, support_y, support_temp, support_hab, nPerm, cores = 1, seed = 1) {
  require(parallel)
  N <- length(temp)
  idx_perm <- matrix(nrow = N, ncol = nPerm)
  for (p in 1:nPerm) {
    set.seed(seed + p)
    idx_perm[,p] <- sample(1:N, size = N, replace = FALSE)
  }
  
  mclapply_function <- function(i) {
    pred_4_FB(y[idx_perm[,i],], temp[idx_perm[,i]], hab[idx_perm[,i]], p_0, p_0_az, alpha, 
              lambda, N_star, sigma, delta, support_y, support_temp, support_hab)
  }
  out <- mclapply(1:nPerm, mclapply_function, mc.cores = cores)
  
  p_n <- out[[1]]$p_n
  p_az <- out[[1]]$p_az
  preqLL <- out[[1]]$preqLL
  for (i in 2:nPerm) {
    p_n <- p_n + out[[i]]$p_n
    p_az <- p_az + out[[i]]$p_az
    preqLL <- preqLL + out[[i]]$preqLL
  }
  p_n <- p_n / nPerm
  p_az <- p_az / nPerm
  preqLL <- preqLL / nPerm
  
  return(list(p_n = p_n, p_az = p_az, preqLL = preqLL,
              alpha = alpha, lambda = lambda, N_star = N_star, sigma = sigma, delta = delta, nPerm = nPerm,
              support_y = support_y, support_temp = support_temp, support_hab = support_hab,
              y = y, temp = temp, hab = hab, p_0 = p_0, p_0_az = p_0_az))
}

pred_res_4_FB_one <- function(p_n_vec, p_n_az, idx_nz, sample_size, N, alpha, lambda, N_star, sigma, delta, 
                              support_y, support_temp, support_hab, interval, seed = 1) {
  
  D <- 4
  lsy <- length(support_y)
  lst <- length(support_temp)
  lsh <- length(support_hab)
  delta_temp <- delta[1]
  delta_hab <- delta[2]
  set.seed(seed)
  
  for (n in 1:N) {
    
    # Sample new observation
    y_tr_obs <- sample(idx_nz, 1, prob = p_n_vec)
    y_tr_obs_idx <- which(idx_nz == y_tr_obs)
    
    # Recover p_{n-1}(y_n, x_n)
    p_n_obs <- p_n_vec[y_tr_obs_idx]
    
    # Recover y_n
    rec_val <- idx2vec(y_tr_obs, D, lsy, lst, lsh)
    y_obs <- rec_val[1:D]
    
    # Ball around y_obs
    ball <- list()
    for (k in 1:D) {
      ball[[k]] <- k0_int(y_obs[k], interval[k], support_y)
    }
    ball[[D+1]] <- support_temp
    ball[[D+2]] <- support_hab
    
    # Index conversion
    ball_tr <- c()
    ball_tr <- mat2idx(as.matrix(expand.grid(ball)), D, lsy, lst, lsh)
    idx_ball_nz <- which(ball_tr %in% idx_nz)
    
    # Compute \gamma(y,x, y_n,x_n) * k_0(y,x | y_n,x_n)
    ker <- mhk_4_ps(p_n_vec[idx_ball_nz], p_n_obs, y_tr_obs, ball_tr, sigma, delta_temp, delta_hab, lsy, lst, lsh)
    
    # Update p_{n-1}
    n_cum <- sample_size + n
    rho_n <- rho(n_cum, lambda, N_star)
    w_n <- 1 / ((alpha + n_cum)^rho_n)
    p_n_vec <- (1 - w_n) * p_n_vec
    p_n_vec[idx_ball_nz] <- p_n_vec[idx_ball_nz] + w_n*ker
    
    # Add \sum_z {\gamma(z,y_n)*k_0(z | y_n)} in y_n
    p_n_vec[y_tr_obs_idx] <- p_n_vec[y_tr_obs_idx] + w_n*(1 - sum(ker))
    
    # Update other values of p_n
    p_n_az <- (1 - w_n) * p_n_az
  }
  return(p_n_vec)
}

pred_res_4_FB <- function(mad_object, N, B, cores = 6, seed = 1) {
  require(parallel)
  alpha <- mad_object$alpha
  sigma <- mad_object$sigma
  delta <- mad_object$delta
  support_y <- mad_object$support_y
  support_temp <- mad_object$support_temp
  support_hab <- mad_object$support_hab
  sample_size <- NROW(mad_object$y)
  lambda <- mad_object$lambda
  N_star <- mad_object$N_star
  
  p_n <- mad_object$p_n
  nz_log <- as.numeric(is_nonzero(p_n))
  idx_nz <- c(1:length(nz_log))[which(nz_log == 1)]
  p_n_vec <- p_n[idx_nz]
  p_az <- mad_object$p_az
  
  D <- 4
  require(crch)
  left <- min(support_y) - 0.5
  interval <- rep(0, D)
  for (k in 1:D) {
    extremes <- crch::qtnorm(c(0.025, 0.975), mean = 50, sd = sigma[k], left = left)
    interval[k] <- max(0, (extremes[2] - extremes[1]) / 2)
  }
  
  mclapply_function <- function(b) {
    seed_b <- seed + b
    pred_res_4_FB_one(p_n_vec, p_az, idx_nz, sample_size, N, alpha, lambda, N_star, sigma, delta, support_y,
                      support_temp, support_hab, interval, seed_b)
  }
  out <- mclapply(1:B, mclapply_function, mc.cores = cores)
  out <- matrix(unlist(out), nrow = length(p_n_vec), ncol = B)
  out
}

pred_cond_mad <- function(pred, pred_az, D, temperature, habitat, lsy) {
  out <- matrix(0, nrow = lsy, ncol = 4)
  for (i in 1:(lsy)) {
    out[i,1] <- sum(pred[i,,,, temperature+1, habitat])
    out[i,2] <- sum(pred[,i,,, temperature+1, habitat])
    out[i,3] <- sum(pred[,,i,, temperature+1, habitat])
    out[i,4] <- sum(pred[,,,i, temperature+1, habitat])
  }
  out[out == 0] <- pred_az
  out <- apply(out, 2, function(x) x / sum(x))
  return(out)
}
pred_cond_mad_hab <- function(pred, pred_az, habitat, lsy) {
  out <- matrix(0, nrow = lsy, ncol = 4)
  for (i in 1:(lsy)) {
    out[i,1] <- sum(pred[i,,,, , habitat])
    out[i,2] <- sum(pred[,i,,, , habitat])
    out[i,3] <- sum(pred[,,i,, , habitat])
    out[i,4] <- sum(pred[,,,i, , habitat])
  }
  out[out == 0] <- pred_az
  out <- apply(out, 2, function(x) x / sum(x))
  return(out)
}
pred_cond_mad_temp <- function(pred, pred_az, temperature, lsy) {
  out <- matrix(0, nrow = lsy, ncol = 4)
  for (i in 1:(lsy)) {
    out[i,1] <- sum(pred[i,,,, temperature,])
    out[i,2] <- sum(pred[,i,,, temperature,])
    out[i,3] <- sum(pred[,,i,, temperature,])
    out[i,4] <- sum(pred[,,,i, temperature,])
  }
  out[out == 0] <- pred_az
  out <- apply(out, 2, function(x) x / sum(x))
  return(out)
}
exp_cond_FB <- function(pred_cond, support) {
  out <- rep(0, NCOL(pred_cond))
  for (i in 1:length(out)) {
    out[i] <- sum(supp * pred_cond[,i])
  }
  return(out)
}
conditional_prediction_mad <- function(pred, pred_az, D, support_y, support_temp, support_hab) {
  out <- matrix(NA, nrow = length(support_temp) * length(support_hab), ncol = D + 2)
  idx <- 0
  for (h in supp_hab) {
    for (t in support_temp) {
      idx <- idx + 1
      pc <- pred_cond_mad(pred, pred_az, D, t, h, length(support_y))
      out[idx, 1:D] <- exp_cond_FB(pc, support_y)
      out[idx, (D+1):(D+2)] <- c(t, h)
    }
  }
  return(out)
}
predict_mad <- function(mad_object, temp_new, hab_new) {
  D <- NCOL(mad_object$y)
  exp_pred <- conditional_prediction_mad(mad_object$p_n, mad_object$p_az, D, 
                                         mad_object$support_y, mad_object$support_temp, 
                                         mad_object$support_hab)
  out <- matrix(nrow = length(temp_new), ncol = NCOL(mad_object$y))
  for (i in 1:NROW(out)) {
    idx <- which(exp_pred[,D+1] == temp_test[i] & exp_pred[,D+2] == hab_test[i])
    out[i,] <- exp_pred[idx, 1:D]
  }
  out
}
cov_4_FB <- function(p_n, D, support) {
  # Marginal predictives
  p_n_marg <- matrix(NA, nrow = length(support), ncol = D)
  for (i in 1:length(supp)) {
    p_n_marg[i,1] <- sum(p_n[i,,,,,])
    p_n_marg[i,2] <- sum(p_n[,i,,,,])
    p_n_marg[i,3] <- sum(p_n[,,i,,,])
    p_n_marg[i,4] <- sum(p_n[,,,i,,])
  }
  p_n_marg <- apply(p_n_marg, 2, function(x) x / sum(x))
  mu_y <- apply(p_n_marg, 2, function(x) sum(x * support))
  
  # Variances
  vv <- rep(NA, D)
  for (i in 1:D) {
    vv[i] <- sum((support - mu_y[i])^2 * p_n_marg[,i])
  }
  
  # Expceted product
  grid_val <- as.matrix(expand.grid(support, support))
  p_n_joint_marg <- matrix(NA, nrow = NROW(grid_val), ncol = 6)
  for (i in 1:NROW(grid_val)) {
    p_n_joint_marg[i,1] <- sum(p_n[grid_val[i,1]+1, grid_val[i,2]+1,,,,])
    p_n_joint_marg[i,2] <- sum(p_n[grid_val[i,1]+1,, grid_val[i,2]+1,,,])
    p_n_joint_marg[i,3] <- sum(p_n[grid_val[i,1]+1,,, grid_val[i,2]+1,,])
    p_n_joint_marg[i,4] <- sum(p_n[, grid_val[i,1]+1, grid_val[i,2]+1,,,])
    p_n_joint_marg[i,5] <- sum(p_n[, grid_val[i,1]+1,, grid_val[i,2]+1,,])
    p_n_joint_marg[i,6] <- sum(p_n[,, grid_val[i,1]+1, grid_val[i,2]+1,,])
  }
  p_n_joint_marg <- apply(p_n_joint_marg, 2, function(x) x / sum(x))
  product_dim <- apply(p_n_joint_marg, 2, function(x) sum(grid_val[,1] * grid_val[,2] * x))
  
  # Covariance matrix
  out <- matrix(NA, D, D)
  out[1,2] <- out[2,1] <- product_dim[1] - mu_y[1]*mu_y[2]
  out[1,3] <- out[3,1] <- product_dim[2] - mu_y[1]*mu_y[3]
  out[1,4] <- out[4,1] <- product_dim[3] - mu_y[1]*mu_y[4]
  out[2,3] <- out[3,2] <- product_dim[4] - mu_y[2]*mu_y[3]
  out[2,4] <- out[4,2] <- product_dim[5] - mu_y[2]*mu_y[4]
  out[3,4] <- out[4,3] <- product_dim[6] - mu_y[3]*mu_y[4]
  diag(out) <- vv
  return(out)
}
cor_4_FB <- function(cov_mad) {
  D <- NROW(cov_mad)
  cr <- matrix(NA, D, D)
  for (i in 1:D) {
    for (j in 1:D) {
      cr[i,j] <- cov_mad[i,j] / sqrt(cov_mad[i,i]*cov_mad[j,j])
    }
  }
  return(cr)
}

preq_4_FB <- function(y, temp, hab, p_0, p_0_az, alpha, lambda, N_star, sigma, delta, support_y, support_temp, support_hab) {
  
  d <- 4
  require(crch)
  left <- min(support_y) - 0.5
  interval <- rep(0, d)
  for (k in 1:d) {
    extremes <- crch::qtnorm(c(0.025, 0.975), mean = 50, sd = sigma[k], left = left)
    interval[k] <- max(0, (extremes[2] - extremes[1]) / 2)
  }
  
  lsy <- length(support_y)
  lst <- length(support_temp)
  lsh <- length(support_hab)
  delta_temp <- delta[1]
  delta_hab <- delta[2]
  
  p_0[y[1,1]+1, y[1,2]+1, y[1,3]+1, y[1,4]+1, temp[1]+1, hab[1]] <- p_0_az
  out_z <- p_0_az
  preqLL <- 0
  
  for (n in 1:NROW(y)) {
    y_obs <- c(y[n,])
    temp_obs <- temp[n]
    hab_obs <- hab[n]
    idx_obs <- y_obs + 1
    
    # Recover p_{n-1}(y_n, x_n)
    p_n_obs <- p_0[idx_obs[1], idx_obs[2], idx_obs[3], idx_obs[4], temp_obs+1, hab_obs]
    if (p_n_obs == 0) p_n_obs <- out_z
    
    # Update prequential log-likelihood
    preqLL <- preqLL + log(p_n_obs) # - log(sum(p_0[,,,, temp_obs+1, hab_obs]))
    
    # Ball around y_obs
    ball <- list()
    for (k in 1:d) {
      ball[[k]] <- k0_int(y_obs[k], interval[k], support_y)
    }
    ball[[d+1]] <- support_temp
    ball[[d+2]] <- support_hab
    
    # Index conversion
    idx_tr <- c()
    # idx_tr <- ball2idx_5(ball, support_hab, support_temp, lsy)
    idx_tr <- mat2idx(as.matrix(expand.grid(ball)), d, lsy, lst, lsh)
    
    # Compute \gamma(y,x, y_n,x_n) * k_0(y,x | y_n,x_n)
    p_0[idx_tr][p_0[idx_tr] == 0] <- out_z
    ker <- mhk_4_FB(p_0[idx_tr], p_n_obs, y_obs, temp_obs, hab_obs, ball,
                    support_temp, support_hab, sigma, delta_temp, delta_hab)
    
    # Update p_{n-1}
    rho_n <- rho(n, lambda, N_star)
    w_n <- 1 / ((alpha + n)^rho_n)
    p_0 <- (1 - w_n) * p_0
    p_0[idx_tr] <- p_0[idx_tr] + w_n*ker
    
    # Add \sum_z {\gamma(z,y_n)*k_0(z | y_n)} in y_n
    p_0[idx_obs[1], idx_obs[2], idx_obs[3], idx_obs[4], temp_obs+1, hab_obs] <- p_0[idx_obs[1], idx_obs[2], idx_obs[3], idx_obs[4], temp_obs+1, hab_obs] + w_n*(1 - sum(ker))
    
    # Update other values of p_n
    out_z <- (1 - w_n) * out_z
  }
  
  preqLL
}
preq_4_FB_perm <- function(y, temp, hab, p_0, p_0_az, alpha, lambda, N_star, sigma, delta, support_y, support_temp, support_hab, nPerm, cores = 1, seed = 1) {
  require(parallel)
  N <- length(temp)
  idx_perm <- matrix(nrow = N, ncol = nPerm)
  for (p in 1:nPerm) {
    set.seed(seed + p)
    idx_perm[,p] <- sample(1:N, size = N, replace = FALSE)
  }
  
  mclapply_function <- function(i) {
    preq_4_FB(y[idx_perm[,i],], temp[idx_perm[,i]], hab[idx_perm[,i]], p_0, p_0_az, alpha, lambda, N_star,
              sigma, delta, support_y, support_temp, support_hab)
  }
  out <- mclapply(1:nPerm, mclapply_function, mc.cores = cores)
  mean(unlist(out))
}

preq_cond_4_FB <- function(y, temp, hab, p_0, p_0_az, alpha, lambda, N_star, sigma, delta, support_y, support_temp, support_hab) {
  
  d <- 4
  require(crch)
  left <- min(support_y) - 0.5
  interval <- rep(0, d)
  for (k in 1:d) {
    extremes <- crch::qtnorm(c(0.025, 0.975), mean = 50, sd = sigma[k], left = left)
    interval[k] <- max(0, (extremes[2] - extremes[1]) / 2)
  }
  
  lsy <- length(support_y)
  lst <- length(support_temp)
  lsh <- length(support_hab)
  delta_temp <- delta[1]
  delta_hab <- delta[2]
  
  p_0[y[1,1]+1, y[1,2]+1, y[1,3]+1, y[1,4]+1, temp[1]+1, hab[1]] <- p_0_az
  out_z <- p_0_az
  preqLL <- 0
  
  for (n in 1:NROW(y)) {
    y_obs <- c(y[n,])
    temp_obs <- temp[n]
    hab_obs <- hab[n]
    idx_obs <- y_obs + 1
    
    # Recover p_{n-1}(y_n, x_n)
    p_n_obs <- p_0[idx_obs[1], idx_obs[2], idx_obs[3], idx_obs[4], temp_obs+1, hab_obs]
    if (p_n_obs == 0) p_n_obs <- out_z
    
    # Update prequential log-likelihood
    preqLL <- preqLL + log(p_n_obs) - log(sum(p_0[,,,, temp_obs+1, hab_obs]))
    
    # Ball around y_obs
    ball <- list()
    for (k in 1:d) {
      ball[[k]] <- k0_int(y_obs[k], interval[k], support_y)
    }
    ball[[d+1]] <- support_temp
    ball[[d+2]] <- support_hab
    
    # Index conversion
    idx_tr <- c()
    # idx_tr <- ball2idx_5(ball, support_hab, support_temp, lsy)
    idx_tr <- mat2idx(as.matrix(expand.grid(ball)), d, lsy, lst, lsh)
    
    # Compute \gamma(y,x, y_n,x_n) * k_0(y,x | y_n,x_n)
    p_0[idx_tr][p_0[idx_tr] == 0] <- out_z
    ker <- mhk_4_FB(p_0[idx_tr], p_n_obs, y_obs, temp_obs, hab_obs, ball,
                    support_temp, support_hab, sigma, delta_temp, delta_hab)
    
    # Update p_{n-1}
    rho_n <- rho(n, lambda, N_star)
    w_n <- 1 / ((alpha + n)^rho_n)
    p_0 <- (1 - w_n) * p_0
    p_0[idx_tr] <- p_0[idx_tr] + w_n*ker
    
    # Add \sum_z {\gamma(z,y_n)*k_0(z | y_n)} in y_n
    p_0[idx_obs[1], idx_obs[2], idx_obs[3], idx_obs[4], temp_obs+1, hab_obs] <- p_0[idx_obs[1], idx_obs[2], idx_obs[3], idx_obs[4], temp_obs+1, hab_obs] + w_n*(1 - sum(ker))
    
    # Update other values of p_n
    out_z <- (1 - w_n) * out_z
  }
  
  preqLL
}
preq_cond_4_FB_perm <- function(y, temp, hab, p_0, p_0_az, alpha, lambda, N_star, sigma, delta, support_y, support_temp, support_hab, nPerm, cores = 1, seed = 1) {
  require(parallel)
  N <- length(temp)
  idx_perm <- matrix(nrow = N, ncol = nPerm)
  for (p in 1:nPerm) {
    set.seed(seed + p)
    idx_perm[,p] <- sample(1:N, size = N, replace = FALSE)
  }
  
  mclapply_function <- function(i) {
    preq_cond_4_FB(y[idx_perm[,i],], temp[idx_perm[,i]], hab[idx_perm[,i]], p_0, p_0_az, alpha, lambda, N_star,
                   sigma, delta, support_y, support_temp, support_hab)
  }
  out <- mclapply(1:nPerm, mclapply_function, mc.cores = cores)
  mean(unlist(out))
}
