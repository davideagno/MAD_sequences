
rel_freq_class <- function(y, X) {
  N <- length(y)
  D <- NCOL(X) + 1
  out <- array(0, dim = rep(2, D))
  for (n in 1:N) {
    y_obs <- y[n]
    x_obs <- c(X[n,])
    idx_obs <- sum(c(y_obs, x_obs) * 2^c(0:10)) + 1
    out[idx_obs] <- out[idx_obs] + 1
  }
  out / N
}
rho <- function(n, lambda, N_star) {
  lambda + (1 - lambda)*exp(-(1/N_star)*n)
}

pred_class_one <- function(y, X, alpha, lambda, N_star, P_0, delta_y, delta_x) {
  P_0 <- c(P_0)
  for (n in 1:NROW(y)) {
    y_obs <- y[n]
    x_obs <- c(X[n,])
    idx_obs <- sum(c(y_obs, x_obs) * 2^c(0:10)) + 1
    p_n_obs <- P_0[idx_obs]
    
    # Compute \gamma(y,x, y_n,x_n) * k_0(y,x | y_n,x_n)
    ker <- mhk_class_2(P_0, p_n_obs, y_obs, x_obs, delta_y, delta_x)
    
    # Update p_{n-1}
    rho_n <- rho(n, lambda, N_star)
    w_n <- 1 / ((alpha + n)^rho_n)
    P_0 <- (1 - w_n)*P_0 + w_n*ker
    
    # Add \sum_z {\gamma(z,y_n)*k_0(z | y_n)} in y_n
    P_0[idx_obs] <- P_0[idx_obs] + w_n*(1 - sum(ker))
  }
  return(P_0)
}
pred_class <- function(y, X, alpha, lambda, N_star, P_0, delta, nPerm, cores = 1, seed = 1) {
  require(parallel)
  N <- length(y)
  idx_perm <- matrix(nrow = N, ncol = nPerm)
  for (p in 1:nPerm) {
    set.seed(seed + p)
    idx_perm[,p] <- sample(1:N, size = N, replace = FALSE)
  }
  
  delta_y <- delta[1]
  delta_x <- delta[2]
  mclapply_function <- function(p) {
    pred_class_one(y[idx_perm[,p]], X[idx_perm[,p],], alpha, lambda, N_star, P_0, delta_y, delta_x)
  }
  out <- mclapply(1:nPerm, mclapply_function, mc.cores = cores)
  out <- matrix(unlist(out), ncol = nPerm)
  out <- c(rowMeans(out))
  out <- array(out, dim = dim(P_0))
  out
}

preq_cond_class_one <- function(y, X, alpha, lambda, N_star, P_0, delta_y, delta_x) {
  P_0 <- c(P_0)
  preq <- 0
  for (n in 1:NROW(y)) {
    y_obs <- y[n]
    x_obs <- c(X[n,])
    idx_obs <- sum(c(y_obs, x_obs) * 2^c(0:10)) + 1
    p_n_obs <- P_0[idx_obs]
    
    # Prequential log-likelihoog
    idx_obs_0 <- sum(c(0, x_obs) * 2^c(0:10)) + 1
    idx_obs_1 <- sum(c(1, x_obs) * 2^c(0:10)) + 1
    p_cond_obs <- P_0[idx_obs] / (P_0[idx_obs_0] + P_0[idx_obs_1])
    preq <- preq + log(p_cond_obs)
    
    # Compute \gamma(y,x, y_n,x_n) * k_0(y,x | y_n,x_n)
    ker <- mhk_class_2(P_0, p_n_obs, y_obs, x_obs, delta_y, delta_x)
    
    # Update p_{n-1}
    rho_n <- rho(n, lambda, N_star)
    w_n <- 1 / ((alpha + n)^rho_n)
    P_0 <- (1 - w_n)*P_0 + w_n*ker
    
    # Add \sum_z {\gamma(z,y_n)*k_0(z | y_n)} in y_n
    P_0[idx_obs] <- P_0[idx_obs] + w_n*(1 - sum(ker))
  }
  return(preq)
}
preq_cond_class <- function(y, X, alpha, lambda, N_star, P_0, delta, nPerm, cores = 1, seed = 1) {
  require(parallel)
  N <- length(y)
  idx_perm <- matrix(nrow = N, ncol = nPerm)
  for (p in 1:nPerm) {
    set.seed(seed + p)
    idx_perm[,p] <- sample(1:N, size = N, replace = FALSE)
  }
  
  delta_y <- delta[1]
  delta_x <- delta[2]
  mclapply_function <- function(p) {
    preq_cond_class_one(y[idx_perm[,p]], X[idx_perm[,p],], alpha, lambda, N_star, P_0, delta_y, delta_x)
  }
  out <- mclapply(1:nPerm, mclapply_function, mc.cores = cores)
  mean(unlist(out))
}

preq_class_one <- function(y, X, alpha, lambda, N_star, P_0, delta_y, delta_x) {
  P_0 <- c(P_0)
  preq <- 0
  for (n in 1:NROW(y)) {
    y_obs <- y[n]
    x_obs <- c(X[n,])
    idx_obs <- sum(c(y_obs, x_obs) * 2^c(0:10)) + 1
    p_n_obs <- P_0[idx_obs]
    
    # Prequential log-likelihoog
    preq <- preq + log(p_n_obs)
    
    # Compute \gamma(y,x, y_n,x_n) * k_0(y,x | y_n,x_n)
    ker <- mhk_class_2(P_0, p_n_obs, y_obs, x_obs, delta_y, delta_x)
    
    # Update p_{n-1}
    rho_n <- rho(n, lambda, N_star,)
    w_n <- 1 / ((alpha + n)^rho_n)
    P_0 <- (1 - w_n)*P_0 + w_n*ker
    
    # Add \sum_z {\gamma(z,y_n)*k_0(z | y_n)} in y_n
    P_0[idx_obs] <- P_0[idx_obs] + w_n*(1 - sum(ker))
  }
  return(preq)
}
preq_class <- function(y, X, alpha, lambda, N_star, P_0, delta, nPerm, cores = 1, seed = 1) {
  require(parallel)
  N <- length(y)
  idx_perm <- matrix(nrow = N, ncol = nPerm)
  for (p in 1:nPerm) {
    set.seed(seed + p)
    idx_perm[,p] <- sample(1:N, size = N, replace = FALSE)
  }
  
  delta_y <- delta[1]
  delta_x <- delta[2]
  mclapply_function <- function(p) {
    preq_class_one(y[idx_perm[,p]], X[idx_perm[,p],], alpha, lambda, N_star, P_0, delta_y, delta_x)
  }
  out <- mclapply(1:nPerm, mclapply_function, mc.cores = cores)
  mean(unlist(out))
}

pred_cond_class <- function(p_n, x_new) {
  idx_0_tr <- sum(c(0, x_new) * 2^c(0:10)) + 1
  idx_1_tr <- sum(c(1, x_new) * 2^c(0:10)) + 1
  prob_0 <- p_n[idx_0_tr]
  prob_1 <- p_n[idx_1_tr]
  prob_1 / (prob_0 + prob_1)
}
predict_mad_class <- function(p_n, X_new) {
  N <- NROW(X_new)
  out <- rep(NA, N)
  for (i in 1:N) {
    out[i] <- pred_cond_class(p_n, X_new[i,])
  }
  out
}

pred_res_class_one <- function(p_n_vec, dim_joint_supp, sample_size, alpha, lambda, N_star, delta_y, delta_x, N, seed = 1) {
  set.seed(seed)
  for(n in 1:N) {
    # Sample new observation
    y_tr_obs <- sample(x = 1:dim_joint_supp, size = 1, prob = p_n_vec)
    rec_obs <- idx2vec_class(y_tr_obs, 11)
    y_obs <- rec_obs[1]
    x_obs <- rec_obs[-1]
    
    # Recover p_{n-1}(y_n, x_n)
    p_n_obs <- p_n_vec[y_tr_obs]
    
    # Compute \gamma(y,x, y_n,x_n) * k_0(y,x | y_n,x_n)
    ker <- mhk_class_2(p_n_vec, p_n_obs, y_obs, x_obs, delta_y, delta_x)
    
    # Update p_{n-1}
    n_ps <- sample_size + n
    rho_n <- rho(n_ps, lambda, N_star)
    w_n <- 1 / ((alpha + n_ps)^rho_n)
    p_n_vec <- (1 - w_n)*p_n_vec + w_n*ker
    
    # Add \sum_z {\gamma(z,y_n)*k_0(z | y_n)} in y_n
    p_n_vec[y_tr_obs] <- p_n_vec[y_tr_obs] + w_n*(1 - sum(ker))
  }
  return(p_n_vec)
}
pred_res_class <- function(p_n, sample_size, alpha, lambda, N_star, delta, N, B, cores = 1, seed = 1) {
  require(parallel)
  p_n <- c(p_n)
  djs <- length(p_n)
  delta_y <- delta[1]
  delta_x <- delta[2]
  
  mclapply_function <- function(b) {
    seed_b <- seed + b
    pred_res_class_one(p_n, djs, sample_size, alpha, lambda, N_star, delta_y, delta_x, N, seed_b)
  }
  out <- mclapply(1:B, mclapply_function, mc.cores = cores)
  out <- matrix(unlist(out), nrow = djs, ncol = B)
  out
}