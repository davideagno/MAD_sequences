

k0_int <- function(x, marg, support) {
  dx <- pmin(max(support), ceiling(x + marg))
  sx <- pmax(min(support), floor(x - marg))
  sx:dx
}

pred_regr_one <- function(y, X, alpha, P_0, sigma, delta_x, support_y) {
  P_0 <- c(P_0)
  d <- NCOL(X) + 1
  lsy <- length(support_y)
  
  require(crch)
  left <- min(support_y) - 0.5
  extremes <- crch::qtnorm(c(0.005, 0.995), mean = 50, sd = sigma, left = left)
  interval <- max(0, (extremes[2] - extremes[1]) / 2)
  
  ball <- list()
  ball[[1]] <- 0
  for (k in 2:d) {
    ball[[k]] <- 0:1
  }
  
  for (n in 1:NROW(y)) {
    y_obs <- y[n]
    x_obs <- c(X[n,])
    idx_obs <- vec2idx_regr(c(y_obs, x_obs), d, lsy)
    p_n_obs <- P_0[idx_obs]
    
    # Ball around y_obs
    ball[[1]] <- k0_int(y_obs, interval, support_y)
    
    # Index conversion
    idx_tr <- c()
    idx_tr <- mat2idx_regr(as.matrix(expand.grid(ball)), d, lsy)
    
    # Compute \gamma(y,x, y_n,x_n) * k_0(y,x | y_n,x_n)
    ker <- mhk_regr(P_0[idx_tr], ball[[1]], p_n_obs, y_obs, x_obs, sigma, delta_x)
    
    # Update p_{n-1}
    w_n <- 1 / (alpha + n)
    P_0 <- (1 - w_n)*P_0
    P_0[idx_tr] <- P_0[idx_tr] + w_n*ker
    
    # Add \sum_z {\gamma(z,y_n)*k_0(z | y_n)} in y_n
    P_0[idx_obs] <- P_0[idx_obs] + w_n*(1 - sum(ker))
  }
  return(P_0)
}
pred_regr <- function(y, X, alpha, P_0, hyper, support_y, nPerm, cores = 1, seed = 1) {
  require(parallel)
  N <- length(y)
  idx_perm <- matrix(nrow = N, ncol = nPerm)
  for (p in 1:nPerm) {
    set.seed(seed + p)
    idx_perm[,p] <- sample(1:N, size = N, replace = FALSE)
  }
  
  sigma <- hyper[1]
  delta_x <- hyper[2]
  mclapply_function <- function(p) {
    pred_regr_one(y[idx_perm[,p]], X[idx_perm[,p],], alpha, P_0, sigma, delta_x, support_y)
  }
  out <- mclapply(1:nPerm, mclapply_function, mc.cores = cores)
  out <- matrix(unlist(out), ncol = nPerm)
  out <- c(rowMeans(out))
  out <- array(out, dim = dim(P_0))
  out
}

preq_regr_one <- function(y, X, alpha, P_0, sigma, delta_x, support_y) {
  P_0 <- c(P_0)
  d <- NCOL(X) + 1
  lsy <- length(support_y)
  
  require(crch)
  left <- min(support_y) - 0.5
  extremes <- crch::qtnorm(c(0.005, 0.995), mean = 50, sd = sigma, left = left)
  interval <- max(0, (extremes[2] - extremes[1]) / 2)
  
  ball <- list()
  ball[[1]] <- 0
  for (k in 2:d) {
    ball[[k]] <- 0:1
  }
  
  preq <- 0
  for (n in 1:NROW(y)) {
    y_obs <- y[n]
    x_obs <- c(X[n,])
    idx_obs <- vec2idx_regr(c(y_obs, x_obs), d, lsy)
    p_n_obs <- P_0[idx_obs]
    
    # Prequential log-likelihoog
    idx_x_obs <- vec2idx_regr(c(0, x_obs), d, lsy)
    p_x_obs <- sum(P_0[idx_x_obs:(idx_x_obs+lsy-1)])
    preq <- preq + log(p_n_obs)
    
    # Ball around y_obs
    ball[[1]] <- k0_int(y_obs, interval, support_y)
    
    # Index conversion
    idx_tr <- c()
    idx_tr <- mat2idx_regr(as.matrix(expand.grid(ball)), d, lsy)
    
    # Compute \gamma(y,x, y_n,x_n) * k_0(y,x | y_n,x_n)
    ker <- mhk_regr(P_0[idx_tr], ball[[1]], p_n_obs, y_obs, x_obs, sigma, delta_x)
    
    # Update p_{n-1}
    w_n <- 1 / (alpha + n)
    P_0 <- (1 - w_n)*P_0
    P_0[idx_tr] <- P_0[idx_tr] + w_n*ker
    
    # Add \sum_z {\gamma(z,y_n)*k_0(z | y_n)} in y_n
    P_0[idx_obs] <- P_0[idx_obs] + w_n*(1 - sum(ker))
  }
  return(preq)
}
preq_regr <- function(y, X, alpha, P_0, hyper, support_y, nPerm, cores = 1, seed = 1) {
  require(parallel)
  N <- length(y)
  idx_perm <- matrix(nrow = N, ncol = nPerm)
  for (p in 1:nPerm) {
    set.seed(seed + p)
    idx_perm[,p] <- sample(1:N, size = N, replace = FALSE)
  }
  
  sigma <- hyper[1]
  delta_x <- hyper[2]
  mclapply_function <- function(p) {
    preq_regr_one(y[idx_perm[,p]], X[idx_perm[,p],], alpha, P_0, sigma, delta_x, support_y)
  }
  out <- mclapply(1:nPerm, mclapply_function, mc.cores = cores)
  mean(unlist(out))
}

preq_cond_regr_one <- function(y, X, alpha, P_0, sigma, delta_x, support_y) {
  P_0 <- c(P_0)
  d <- NCOL(X) + 1
  lsy <- length(support_y)
  
  require(crch)
  left <- min(support_y) - 0.5
  extremes <- crch::qtnorm(c(0.005, 0.995), mean = 50, sd = sigma, left = left)
  interval <- max(0, (extremes[2] - extremes[1]) / 2)
  
  ball <- list()
  ball[[1]] <- 0
  for (k in 2:d) {
    ball[[k]] <- 0:1
  }
  
  preq <- 0
  for (n in 1:NROW(y)) {
    y_obs <- y[n]
    x_obs <- c(X[n,])
    idx_obs <- vec2idx_regr(c(y_obs, x_obs), d, lsy)
    p_n_obs <- P_0[idx_obs]
    
    # Prequential log-likelihoog
    idx_x_obs <- vec2idx_regr(c(0, x_obs), d, lsy)
    p_x_obs <- sum(P_0[idx_x_obs:(idx_x_obs+lsy-1)])
    preq <- preq + log(p_n_obs) - log(p_x_obs)
    
    # Ball around y_obs
    ball[[1]] <- k0_int(y_obs, interval, support_y)
    
    # Index conversion
    idx_tr <- c()
    idx_tr <- mat2idx_regr(as.matrix(expand.grid(ball)), d, lsy)
    
    # Compute \gamma(y,x, y_n,x_n) * k_0(y,x | y_n,x_n)
    ker <- mhk_regr(P_0[idx_tr], ball[[1]], p_n_obs, y_obs, x_obs, sigma, delta_x)
    
    # Update p_{n-1}
    w_n <- 1 / (alpha + n)
    P_0 <- (1 - w_n)*P_0
    P_0[idx_tr] <- P_0[idx_tr] + w_n*ker
    
    # Add \sum_z {\gamma(z,y_n)*k_0(z | y_n)} in y_n
    P_0[idx_obs] <- P_0[idx_obs] + w_n*(1 - sum(ker))
  }
  return(preq)
}
preq_cond_regr <- function(y, X, alpha, P_0, hyper, support_y, nPerm, cores = 1, seed = 1) {
  require(parallel)
  N <- length(y)
  idx_perm <- matrix(nrow = N, ncol = nPerm)
  for (p in 1:nPerm) {
    set.seed(seed + p)
    idx_perm[,p] <- sample(1:N, size = N, replace = FALSE)
  }
  
  sigma <- hyper[1]
  delta_x <- hyper[2]
  mclapply_function <- function(p) {
    preq_cond_regr_one(y[idx_perm[,p]], X[idx_perm[,p],], alpha, P_0, sigma, delta_x, support_y)
  }
  out <- mclapply(1:nPerm, mclapply_function, mc.cores = cores)
  mean(unlist(out))
}

preq_cond_05_regr_one <- function(y, X, alpha, P_0, sigma, delta_x, support_y) {
  P_0 <- c(P_0)
  d <- NCOL(X) + 1
  lsy <- length(support_y)
  
  require(crch)
  left <- min(support_y) - 0.5
  extremes <- crch::qtnorm(c(0.005, 0.995), mean = 50, sd = sigma, left = left)
  interval <- max(0, (extremes[2] - extremes[1]) / 2)
  
  ball <- list()
  ball[[1]] <- 0
  for (k in 2:d) {
    ball[[k]] <- 0:1
  }
  
  preq <- 0
  nTreshold <- 0.5*NROW(y)
  
  for (n in 1:NROW(y)) {
    y_obs <- y[n]
    x_obs <- c(X[n,])
    idx_obs <- vec2idx_regr(c(y_obs, x_obs), d, lsy)
    p_n_obs <- P_0[idx_obs]
    
    # Prequential log-likelihoog
    idx_x_obs <- vec2idx_regr(c(0, x_obs), d, lsy)
    p_x_obs <- sum(P_0[idx_x_obs:(idx_x_obs+lsy-1)])
    if (n >= nTreshold) preq <- preq + log(p_n_obs) - log(p_x_obs)
    
    # Ball around y_obs
    ball[[1]] <- k0_int(y_obs, interval, support_y)
    
    # Index conversion
    idx_tr <- c()
    idx_tr <- mat2idx_regr(as.matrix(expand.grid(ball)), d, lsy)
    
    # Compute \gamma(y,x, y_n,x_n) * k_0(y,x | y_n,x_n)
    ker <- mhk_regr(P_0[idx_tr], ball[[1]], p_n_obs, y_obs, x_obs, sigma, delta_x)
    
    # Update p_{n-1}
    w_n <- 1 / (alpha + n)
    P_0 <- (1 - w_n)*P_0
    P_0[idx_tr] <- P_0[idx_tr] + w_n*ker
    
    # Add \sum_z {\gamma(z,y_n)*k_0(z | y_n)} in y_n
    P_0[idx_obs] <- P_0[idx_obs] + w_n*(1 - sum(ker))
  }
  return(preq)
}
preq_cond_05_regr <- function(y, X, alpha, P_0, hyper, support_y, nPerm, cores = 1, seed = 1) {
  require(parallel)
  N <- length(y)
  idx_perm <- matrix(nrow = N, ncol = nPerm)
  for (p in 1:nPerm) {
    set.seed(seed + p)
    idx_perm[,p] <- sample(1:N, size = N, replace = FALSE)
  }
  
  sigma <- hyper[1]
  delta_x <- hyper[2]
  mclapply_function <- function(p) {
    preq_cond_05_regr_one(y[idx_perm[,p]], X[idx_perm[,p],], alpha, P_0, sigma, delta_x, support_y)
  }
  out <- mclapply(1:nPerm, mclapply_function, mc.cores = cores)
  mean(unlist(out))
}

pred_cond_regr <- function(p_n, x_new, d = 11, lsy = 101) {
  out <- rep(NA, lsy)
  for (i in 1:lsy) {
    idx_y_obs <- vec2idx_regr(c(i-1, x_new), d, lsy)
    out[i] <- p_n[idx_y_obs]
  }
  out / sum(out)
} 
predict_mad_regr <- function(p_n, X_new, d = 11, support_y = 0:100) {
  N <- NROW(X_new)
  lsy <- length(support_y)
  out <- rep(NA, N)
  for (i in 1:N) {
    p_cond <- pred_cond_regr(p_n, X_new[i,], d, lsy)
    out[i] <- sum(support_y * p_cond)
  }
  out
}

pred_res_regr_one <- function(p_n_vec, dim_joint_supp, sample_size, alpha, sigma, delta_x, support_y, ball, interval, N, seed = 1) {
  lsy <- length(support_y)
  set.seed(seed)
  for (n in 1:N) {
    # Sample new observation
    idx_obs <- sample(x = 1:dim_joint_supp, size = 1, prob = p_n_vec)
    rec_obs <- idx2vec_regr(idx_obs, 11, lsy)
    y_obs <- rec_obs[1]
    x_obs <- rec_obs[-1]
    
    # Recover p_{n-1}(y_n, x_n)
    p_n_obs <- p_n_vec[idx_obs]
    
    # Ball around y_obs
    ball[[1]] <- k0_int(y_obs, interval, support_y)
    
    # Index conversion
    idx_tr <- c()
    idx_tr <- mat2idx_regr(as.matrix(expand.grid(ball)), 11, lsy)
    
    # Compute \gamma(y,x, y_n,x_n) * k_0(y,x | y_n,x_n)
    ker <- mhk_regr(p_n_vec[idx_tr], ball[[1]], p_n_obs, y_obs, x_obs, sigma, delta_x)
    
    # Update p_{n-1}
    w_n <- 1 / (alpha + sample_size + n)
    p_n_vec <- (1 - w_n)*p_n_vec
    p_n_vec[idx_tr] <- p_n_vec[idx_tr] + w_n*ker
    
    # Add \sum_z {\gamma(z,y_n)*k_0(z | y_n)} in y_n
    p_n_vec[idx_obs] <- p_n_vec[idx_obs] + w_n*(1 - sum(ker))
  }
  return(p_n_vec)
}
pred_res_regr <- function(p_n, sample_size, alpha, sigma, delta_x, support_y, N, B, cores = 1, seed = 1, dp = FALSE) {
  require(parallel)
  p_n <- c(p_n)
  djs <- length(p_n)
  
  if (dp == TRUE) {
    interval <- 0
  } else {
    require(crch)
    left <- min(support_y) - 0.5
    extremes <- crch::qtnorm(c(0.005, 0.995), mean = 50, sd = sigma, left = left)
    interval <- max(0, (extremes[2] - extremes[1]) / 2)  
  }
  
  ball <- list()
  ball[[1]] <- 0
  for (k in 2:11) {
    ball[[k]] <- 0:1
  }
  
  mclapply_function <- function(b) {
    seed_b <- seed + b
    pred_res_regr_one(p_n, djs, sample_size, alpha, sigma, delta_x, support_y, ball, interval, N, seed_b)
  }
  out <- mclapply(1:B, mclapply_function, mc.cores = cores)
  out <- matrix(unlist(out), nrow = djs, ncol = B)
  out
}

cond_exp_val_ps <- function(ps, x_new, support_y = 0:100) {
  lsy <- length(support_y)
  p_cond <- matrix(nrow = lsy, ncol = NCOL(ps))
  for (i in 1:NCOL(p_cond)) {
    p_cond[,i] <- pred_cond_regr(ps[,i], x_new, d = 11, lsy = lsy)
  }
  out <- rep(NA, NCOL(ps)) 
  for (j in 1:NCOL(ps)) {
    out[j] <- sum(support_y * p_cond[,j])
  }
  out
}
rel_freq_regr <- function(y, X, support_y) {
  N <- length(y)
  out <- array(0, dim = c(length(support_y), rep(2, NCOL(X))))
  for (n in 1:N) {
    idx_n <- vec2idx_regr(c(y[n], X[n,]), d = NCOL(X)+1, lsy = length(support_y))
    out[idx_n] <- out[idx_n] + 1
  }
  out / N
}