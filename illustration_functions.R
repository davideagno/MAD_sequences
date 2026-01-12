emp_hpd_lower <- function(x, conf_level = 0.95) robustBLME::hpd(x, conf_level)[1]
emp_hpd_upper <- function(x, conf_level = 0.95) robustBLME::hpd(x, conf_level)[2]
pred_mad <- function(y, p_0, alpha, lambda, N_star, sd_rg, support, nPerm, seed = 1, cores = 5) {
  idx_perm <- y_poss_perm <- matrix(nrow = length(y), ncol = nPerm) 
  set.seed(seed)
  for (p in 1:nPerm) {
    idx_perm[,p] <- sample(1:length(y), size = length(y), replace = FALSE)
    y_poss_perm[,p] <- y[idx_perm[,p]]
  }
  
  out <- mclapply(1:nPerm, mclapply_function_cpp, mc.cores = cores, y = y_poss_perm, 
                  p_0 = p_0, sd_rg = sd_rg, support = support, lambda = lambda, 
                  N_star = N_star, alpha = alpha)
  out <- rowMeans(matrix(unlist(out), nrow = length(support), ncol = nPerm))
  out
}
preq_ll <- function(sd_rg, y, p_0, alpha, lambda, N_star, support, nPerm, seed = 1, cores = 5) {
  idx_perm <- y_poss_perm <- matrix(nrow = length(y), ncol = nPerm) 
  set.seed(seed)
  for (p in 1:nPerm) {
    idx_perm[,p] <- sample(1:length(y), size = length(y), replace = FALSE)
    y_poss_perm[,p] <- y[idx_perm[,p]]
  }
  
  out <- mclapply(1:nPerm, mclapply_function_pq_cpp, mc.cores = cores, y = y_poss_perm,
                  p_0 = p_0, sd_rg = sd_rg, support = support, lambda = lambda, 
                  N_star = N_star, alpha = alpha)
  out <- mean(unlist(out))
  out
}
pred_res <- function(y, p_n, alpha, lambda, N_star, support, N, B, sd_rg, seed = 1, cores = 5) {
  require(parallel)
  m <- length(y)
  
  set.seed(seed)
  out <- mclapply(1:B, mclapply_function_ps_cpp, mc.cores = cores,
                  p_n = p_n, m = m, N = N, sd_rg = sd_rg, lambda = lambda,
                  N_star = N_star, alpha = alpha, support = support)
  out <- matrix(unlist(out), nrow = length(support), ncol = B)
  out
}
mclapply_function_ps_conv_val <- function(b, p_n, m, N, sd_rg, lambda, N_star, alpha, support, val, seed) {
  set.seed(seed*b)
  mclapply_function_ps_conv_val_cpp(val, b, p_n, m, N, sd_rg, lambda, N_star, alpha, support)
}
pred_res_conv_val <- function(val, y, p_n, alpha, lambda, N_star, support, N, B, sd_rg, seed = 1, cores = 5) {
  require(parallel)
  m <- length(y)
  
  out <- mclapply(1:B, mclapply_function_ps_conv_val, mc.cores = cores,
                  p_n = p_n, m = m, N = N, sd_rg = sd_rg, lambda = lambda,
                  N_star = N_star, alpha = alpha, support = support, val = val, seed = seed)
  out <- matrix(unlist(out), nrow = N - m, ncol = B)
  out
}

mclapply_function_ps_conv <- function(b, p_n, m, N, sd_rg, lambda, N_star, alpha, support, seed) {
  set.seed(seed*b)
  mclapply_ps_conv_cpp(p_n, m, N, sd_rg, lambda, N_star, alpha, support)
}
pred_res_conv <- function(y, p_n, alpha, lambda, N_star, support, N, B, sd_rg, seed = 1, cores = 5) {
  require(parallel)
  m <- length(y)
  
  out <- mclapply(1:B, mclapply_function_ps_conv, mc.cores = cores,
                  p_n = p_n, m = m, N = N, sd_rg = sd_rg, lambda = lambda,
                  N_star = N_star, alpha = alpha, support = support, seed = seed)
  out <- matrix(unlist(out), nrow = N - m, ncol = B)
  out
}

pred_mad_cons <- function(y, sample_size, p_0, alpha, lambda, N_star, sd_rg, support, nPerm, seed = 1, cores = 5) {
  idx_perm <- y_poss_perm <- matrix(nrow = length(y), ncol = nPerm) 
  set.seed(seed)
  for (p in 1:nPerm) {
    idx_perm[,p] <- sample(1:length(y), size = length(y), replace = FALSE)
    y_poss_perm[,p] <- y[idx_perm[,p]]
  }
  
  out <- mclapply(1:nPerm, mclapply_function_cons_cpp, mc.cores = cores, y = y_poss_perm, 
                  p_0 = p_0, sd_rg = sd_rg, support = support, lambda = lambda, 
                  N_star = N_star, alpha = alpha, sample_size = sample_size)
  out <- rowMeans(matrix(unlist(out), nrow = length(support), ncol = nPerm))
  out
}
