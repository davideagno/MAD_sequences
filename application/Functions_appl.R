
mclapply_function_ps <- function(b, p_n_vec, p_n_az, idx_nz, sample_size, N, 
                                 alpha, lambda, N_star, sigma, delta, support_y,
                                 support_temp, support_hab, interval, seed = 1) {
  set.seed(seed + b)
  pred_res_4_FB_one(p_n_vec, p_n_az, idx_nz, sample_size, N, alpha, lambda,
                    N_star, sigma, delta, support_y, support_temp,
                    support_hab, interval)
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
  
  out <- mclapply(1:B, mclapply_function_ps, mc.cores = cores,
                  p_n_vec = p_n_vec, p_n_az = p_az, idx_nz = idx_nz, sample_size = sample_size,
                  N = N, alpha = alpha, lambda = lambda, N_star = N_star,
                  sigma = sigma, delta = delta, support_y = support_y,
                  support_temp, support_hab, interval = interval, seed = seed)
  out <- matrix(unlist(out), nrow = length(p_n_vec), ncol = B)
  out
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