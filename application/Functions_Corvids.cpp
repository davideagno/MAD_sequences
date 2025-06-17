#include <RcppDist.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppDist)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double log_drg(double x, double mean, double sd, double left = -0.5) {
  double upp = std::numeric_limits<float>::infinity();
  return log(p_truncnorm(x+0.5, mean, sd, left, upp) - p_truncnorm(x-0.5, mean, sd, left, upp));
}

// [[Rcpp::export]]
double p_min(double x) {
  NumericVector out = NumericVector::create(0,x);
  return min(out);
}

// [[Rcpp::export]]
double log_k0_cat(double x, int size, double x_obs, double delta) {
  double out = 0;
  if (x == x_obs) out = 1.0-delta;
  else out = delta / (size - 1);
  return(log(out));
}

//[[Rcpp::export]]
double idx_p(IntegerVector y, int d, int lsy, int lsh, int lst) {
  NumericVector y_idx(d);
  for (int i = 0; i < d; i++) {
    y_idx[i] = y[i] * pow(lsy, i);
  }
  double sum_y_idx = sum(y_idx);
  double t_idx = y[d] * pow(lsy, d);
  double h_idx = (y[d+1] - 1.0) * lst * pow(lsy, d);
  return sum_y_idx + t_idx + h_idx + 1;
}

//[[Rcpp::export]]
NumericVector mat2idx(IntegerMatrix mat, int d, int lsy, int lst, int lsh) {
  int N = mat.nrow();
  NumericVector out(N);
  for (int i = 0; i < N; i++) {
    out(i) = idx_p(mat(i,_), d, lsy, lsh, lst);
  }
  return out;
}


// [[Rcpp::export]]
NumericVector mhk_5_FB(NumericVector p_n, double p_n_obs,
                       NumericVector y, double temp, double hab, 
                       List ball, NumericVector support_temp, NumericVector support_hab, 
                       NumericVector sigma, double delta_temp, double delta_hab, double left = -0.5) {
  
  NumericVector ball_0 = ball[0];
  NumericVector ball_1 = ball[1];
  NumericVector ball_2 = ball[2];
  NumericVector ball_3 = ball[3];
  NumericVector ball_4 = ball[4];
  int lst = support_temp.length();
  int lsh = support_hab.length();
  
  int idx_p = 0;
  NumericVector out(p_n.length());
  
  double log_k0_hab_n, log_k0_n_hab;
  double log_k0_temp_n, log_k0_n_temp;
  double log_k0_y_n_4, log_k0_n_y_4;
  double log_k0_y_n_3, log_k0_n_y_3;
  double log_k0_y_n_2, log_k0_n_y_2;
  double log_k0_y_n_1, log_k0_n_y_1;
  double log_k0_y_n_0, log_k0_n_y_0;
  double log_k0_yx_n, log_k0_n_yx, log_p_n_xy, log_gamma, log_gamma_2;
  
  for (int hab_val = min(support_hab); hab_val <= max(support_hab); hab_val++) {
    log_k0_hab_n = log_k0_cat(hab_val, lsh, hab, delta_hab);
    log_k0_n_hab = log_k0_cat(hab, lsh, hab_val, delta_hab);
    for (int temp_val = min(support_temp); temp_val <= max(support_temp); temp_val++) {
      log_k0_temp_n = log_k0_cat(temp_val, lst, temp, delta_temp);
      log_k0_n_temp = log_k0_cat(temp, lst, temp_val, delta_temp);
      for (int y_val_4 = min(ball_4); y_val_4 <= max(ball_4); y_val_4++) {
        log_k0_y_n_4 = log_drg(y_val_4, y[4], sigma[4], left);
        log_k0_n_y_4 = log_drg(y[4], y_val_4, sigma[4], left);
        for (int y_val_3 = min(ball_3); y_val_3 <= max(ball_3); y_val_3++) {
          log_k0_y_n_3 = log_drg(y_val_3, y[3], sigma[3], left);
          log_k0_n_y_3 = log_drg(y[3], y_val_3, sigma[3], left);
          for (int y_val_2 = min(ball_2); y_val_2 <= max(ball_2); y_val_2++) {
            log_k0_y_n_2 = log_drg(y_val_2, y[2], sigma[2], left);
            log_k0_n_y_2 = log_drg(y[2], y_val_2, sigma[2], left);
            for (int y_val_1 = min(ball_1); y_val_1 <= max(ball_1); y_val_1++) {
              log_k0_y_n_1 = log_drg(y_val_1, y[1], sigma[1], left);
              log_k0_n_y_1 = log_drg(y[1], y_val_1, sigma[1], left);
              for (int y_val_0 = min(ball_0); y_val_0 <= max(ball_0); y_val_0++) {
                log_k0_y_n_0 = log_drg(y_val_0, y[0], sigma[0], left);
                log_k0_n_y_0 = log_drg(y[0], y_val_0, sigma[0], left);
                
                // Compute k_0(y,x | y_n,x_n) and k_0(y_n,x_n | y,x)
                log_k0_yx_n = log_k0_y_n_0 + log_k0_y_n_1 + log_k0_y_n_2 + log_k0_y_n_3 +
                  log_k0_y_n_4 + log_k0_temp_n + log_k0_hab_n;
                log_k0_n_yx = log_k0_n_y_0 + log_k0_n_y_1 + log_k0_n_y_2 + log_k0_n_y_3 +
                  log_k0_n_y_4 + log_k0_n_temp + log_k0_n_hab;
                
                // Recover p_{n-1}(y,x)
                log_p_n_xy = log(p_n[idx_p]);
                
                // Compute \gamma_{n-1}(y,x, y_n,x_n)
                log_gamma = log_p_n_xy + log_k0_n_yx - log(p_n_obs) - log_k0_yx_n;
                log_gamma_2 = p_min(log_gamma);
                
                // Compute \gamma_{n-1}(y,x, y_n,x_n)*k_0(y,x | y_n,x_n)
                out[idx_p] = exp(log_gamma_2 + log_k0_yx_n);
                
                // Update the index
                idx_p++;
              } } } } } } }
  
  return(out);
}


// [[Rcpp::export]]
NumericVector idx2vec(int idx, int D, int lsy, int lst, int lsh) {
  NumericVector out(D+2);
  int num = idx - 1;
  for (int i = 0; i < D; i++) {
    out[i] = num % lsy;
    num = num / lsy;
  }
  out[D] = num % lst;
  num = num / lst;
  out[D+1] = (num % lsh) + 1;
  return(out);
}

// [[Rcpp::export]]
NumericVector mhk_4_ps(NumericVector p_n_vec, double p_n_obs, int idx_tr_obs,
                       IntegerVector idx_tr, NumericVector sigma, double delta_temp,
                       double delta_hab, int lsy, int lst, int lsh, double left = -0.5) {
  
  int D = 4;
  int n = p_n_vec.length();
  NumericVector out(n);
  
  NumericVector val;
  double temp_val, hab_val;
  double y_0_val, y_1_val, y_2_val, y_3_val;
  
  NumericVector val_obs = idx2vec(idx_tr_obs, D, lsy, lst, lsh);
  double y_0 = val_obs[0];
  double y_1 = val_obs[1];
  double y_2 = val_obs[2];
  double y_3 = val_obs[3];
  double temp = val_obs[4];
  double hab = val_obs[5];
  
  double log_k0_hab_n, log_k0_n_hab;
  double log_k0_temp_n, log_k0_n_temp;
  double log_k0_y_n_3, log_k0_n_y_3;
  double log_k0_y_n_2, log_k0_n_y_2;
  double log_k0_y_n_1, log_k0_n_y_1;
  double log_k0_y_n_0, log_k0_n_y_0;
  double log_k0_yx_n, log_k0_n_yx, log_p_n_xy, log_gamma, log_gamma_2;
  
  for (int i = 0; i < n; i++) {
    val = idx2vec(idx_tr[i], D, lsy, lst, lsh);
    y_0_val = val[0];
    y_1_val = val[1];
    y_2_val = val[2];
    y_3_val = val[3];
    temp_val = val[4];
    hab_val = val[5];
    
    log_k0_y_n_0 = log_drg(y_0_val, y_0, sigma[0], left);
    log_k0_n_y_0 = log_drg(y_0, y_0_val, sigma[0], left);
    log_k0_y_n_1 = log_drg(y_1_val, y_1, sigma[1], left);
    log_k0_n_y_1 = log_drg(y_1, y_1_val, sigma[1], left);
    log_k0_y_n_2 = log_drg(y_2_val, y_2, sigma[2], left);
    log_k0_n_y_2 = log_drg(y_2, y_2_val, sigma[2], left);
    log_k0_y_n_3 = log_drg(y_3_val, y_3, sigma[3], left);
    log_k0_n_y_3 = log_drg(y_3, y_3_val, sigma[3], left);
    log_k0_temp_n = log_k0_cat(temp_val, lst, temp, delta_temp);
    log_k0_n_temp = log_k0_cat(temp, lst, temp_val, delta_temp);
    log_k0_hab_n = log_k0_cat(hab_val, lsh, hab, delta_hab);
    log_k0_n_hab = log_k0_cat(hab, lsh, hab_val, delta_hab);
    
    // Compute k_0(y,x | y_n,x_n) and k_0(y_n,x_n | y,x)
    log_k0_yx_n = log_k0_y_n_0 + log_k0_y_n_1 + log_k0_y_n_2 + log_k0_y_n_3 +
      log_k0_temp_n + log_k0_hab_n;
    log_k0_n_yx = log_k0_n_y_0 + log_k0_n_y_1 + log_k0_n_y_2 + log_k0_n_y_3 +
      log_k0_n_temp + log_k0_n_hab;
    
    // Recover p_{n-1}(y,x)
    log_p_n_xy = log(p_n_vec[i]);
    
    // Compute \gamma_{n-1}(y,x, y_n,x_n)
    log_gamma = log_p_n_xy + log_k0_n_yx - log(p_n_obs) - log_k0_yx_n;
    log_gamma_2 = p_min(log_gamma);
    
    // Compute \gamma_{n-1}(y,x, y_n,x_n)*k_0(y,x | y_n,x_n)
    out[i] = exp(log_gamma_2 + log_k0_yx_n);
  }
  
  return(out);
}


// [[Rcpp::export]]
NumericVector mhk_4_FB(NumericVector p_n, double p_n_obs,
                       NumericVector y, double temp, double hab, 
                       List ball, NumericVector support_temp, NumericVector support_hab, 
                       NumericVector sigma, double delta_temp, double delta_hab, double left = -0.5) {
  
  NumericVector ball_0 = ball[0];
  NumericVector ball_1 = ball[1];
  NumericVector ball_2 = ball[2];
  NumericVector ball_3 = ball[3];
  int lst = support_temp.length();
  int lsh = support_hab.length();
  
  int idx_p = 0;
  NumericVector out(p_n.length());
  
  double log_k0_hab_n, log_k0_n_hab;
  double log_k0_temp_n, log_k0_n_temp;
  double log_k0_y_n_3, log_k0_n_y_3;
  double log_k0_y_n_2, log_k0_n_y_2;
  double log_k0_y_n_1, log_k0_n_y_1;
  double log_k0_y_n_0, log_k0_n_y_0;
  double log_k0_yx_n, log_k0_n_yx, log_p_n_xy, log_gamma, log_gamma_2;
  
  for (int hab_val = min(support_hab); hab_val <= max(support_hab); hab_val++) {
    log_k0_hab_n = log_k0_cat(hab_val, lsh, hab, delta_hab);
    log_k0_n_hab = log_k0_cat(hab, lsh, hab_val, delta_hab);
    for (int temp_val = min(support_temp); temp_val <= max(support_temp); temp_val++) {
      log_k0_temp_n = log_k0_cat(temp_val, lst, temp, delta_temp);
      log_k0_n_temp = log_k0_cat(temp, lst, temp_val, delta_temp);
      for (int y_val_3 = min(ball_3); y_val_3 <= max(ball_3); y_val_3++) {
        log_k0_y_n_3 = log_drg(y_val_3, y[3], sigma[3], left);
        log_k0_n_y_3 = log_drg(y[3], y_val_3, sigma[3], left);
        for (int y_val_2 = min(ball_2); y_val_2 <= max(ball_2); y_val_2++) {
          log_k0_y_n_2 = log_drg(y_val_2, y[2], sigma[2], left);
          log_k0_n_y_2 = log_drg(y[2], y_val_2, sigma[2], left);
          for (int y_val_1 = min(ball_1); y_val_1 <= max(ball_1); y_val_1++) {
            log_k0_y_n_1 = log_drg(y_val_1, y[1], sigma[1], left);
            log_k0_n_y_1 = log_drg(y[1], y_val_1, sigma[1], left);
            for (int y_val_0 = min(ball_0); y_val_0 <= max(ball_0); y_val_0++) {
              log_k0_y_n_0 = log_drg(y_val_0, y[0], sigma[0], left);
              log_k0_n_y_0 = log_drg(y[0], y_val_0, sigma[0], left);
              
              // Compute k_0(y,x | y_n,x_n) and k_0(y_n,x_n | y,x)
              log_k0_yx_n = log_k0_y_n_0 + log_k0_y_n_1 + log_k0_y_n_2 + log_k0_y_n_3 +
                log_k0_temp_n + log_k0_hab_n;
              log_k0_n_yx = log_k0_n_y_0 + log_k0_n_y_1 + log_k0_n_y_2 + log_k0_n_y_3 +
                log_k0_n_temp + log_k0_n_hab;
              
              // Recover p_{n-1}(y,x)
              log_p_n_xy = log(p_n[idx_p]);
              
              // Compute \gamma_{n-1}(y,x, y_n,x_n)
              log_gamma = log_p_n_xy + log_k0_n_yx - log(p_n_obs) - log_k0_yx_n;
              log_gamma_2 = p_min(log_gamma);
              
              // Compute \gamma_{n-1}(y,x, y_n,x_n)*k_0(y,x | y_n,x_n)
              out[idx_p] = exp(log_gamma_2 + log_k0_yx_n);
              
              // Update the index
              idx_p++;
            } } } } } }
  
  return(out);
}

