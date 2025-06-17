#include <RcppDist.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppDist)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double log_k0_cat(double x, int size, double x_obs, double delta) {
  double out = 0;
  if (x == x_obs) out = 1.0-delta;
  else out = delta / (size - 1);
  return(log(out));
}

// [[Rcpp::export]]
double log_k0_drg(double x, double mean, double sd, double left = -0.5) {
  double upp = std::numeric_limits<float>::infinity();
  return log(p_truncnorm(x+0.5, mean, sd, left, upp) - p_truncnorm(x-0.5, mean, sd, left, upp));
}

// [[Rcpp::export]]
double p_min(double x) {
  NumericVector out = NumericVector::create(0,x);
  return min(out);
}

// [[Rcpp::export]]
NumericVector idx2vec_regr(int idx, int D, int lsy) {
  NumericVector out(D);
  int num = idx - 1;
  out[0] = num % lsy;
  num = num / lsy;
  for (int i = 1; i < D-1; i++) {
    out[i] = num % 2;
    num = num / 2;
  }
  out[D-1] = num % 2;
  return(out);
}


// [[Rcpp::export]]
NumericVector mhk_regr(NumericVector p_n, NumericVector ball_y, double p_n_obs,
                       double y, NumericVector x, double sigma, double delta_x, double left = -0.5) {
  
  double log_k0_y_n, log_k0_n_y;
  double log_k0_x_n_0, log_k0_n_x_0, log_k0_x_n_1, log_k0_n_x_1;
  double log_k0_x_n_2, log_k0_n_x_2, log_k0_x_n_3, log_k0_n_x_3;
  double log_k0_x_n_4, log_k0_n_x_4, log_k0_x_n_5, log_k0_n_x_5;
  double log_k0_x_n_6, log_k0_n_x_6, log_k0_x_n_7, log_k0_n_x_7;
  double log_k0_x_n_8, log_k0_n_x_8, log_k0_x_n_9, log_k0_n_x_9;
  double log_k0_yx_n, log_k0_n_yx, log_p_n_xy, log_gamma, log_gamma_2;
  
  int idx_p = 0;
  NumericVector out(p_n.length());
  
  for (int x9 = 0; x9 <= 1; x9++) {
    log_k0_x_n_9 = log_k0_cat(x9, 2, x[9], delta_x);
    log_k0_n_x_9 = log_k0_cat(x[9], 2, x9, delta_x);
    for (int x8 = 0; x8 <= 1; x8++) {
      log_k0_x_n_8 = log_k0_cat(x8, 2, x[8], delta_x);
      log_k0_n_x_8 = log_k0_cat(x[8], 2, x8, delta_x);
      for (int x7 = 0; x7 <= 1; x7++) {
        log_k0_x_n_7 = log_k0_cat(x7, 2, x[7], delta_x);
        log_k0_n_x_7 = log_k0_cat(x[7], 2, x7, delta_x);
        for (int x6 = 0; x6 <= 1; x6++) {
          log_k0_x_n_6 = log_k0_cat(x6, 2, x[6], delta_x);
          log_k0_n_x_6 = log_k0_cat(x[6], 2, x6, delta_x);
          for (int x5 = 0; x5 <= 1; x5++) {
            log_k0_x_n_5 = log_k0_cat(x5, 2, x[5], delta_x);
            log_k0_n_x_5 = log_k0_cat(x[5], 2, x5, delta_x);
            for (int x4 = 0; x4 <= 1; x4++) {
              log_k0_x_n_4 = log_k0_cat(x4, 2, x[4], delta_x);
              log_k0_n_x_4 = log_k0_cat(x[4], 2, x4, delta_x);
              for (int x3 = 0; x3 <= 1; x3++) {
                log_k0_x_n_3 = log_k0_cat(x3, 2, x[3], delta_x);
                log_k0_n_x_3 = log_k0_cat(x[3], 2, x3, delta_x);
                for (int x2 = 0; x2 <= 1; x2++) {
                  log_k0_x_n_2 = log_k0_cat(x2, 2, x[2], delta_x);
                  log_k0_n_x_2 = log_k0_cat(x[2], 2, x2, delta_x);
                  for (int x1 = 0; x1 <= 1; x1++) {
                    log_k0_x_n_1 = log_k0_cat(x1, 2, x[1], delta_x);
                    log_k0_n_x_1 = log_k0_cat(x[1], 2, x1, delta_x);
                    for (int x0 = 0; x0 <= 1; x0++) {
                      log_k0_x_n_0 = log_k0_cat(x0, 2, x[0], delta_x);
                      log_k0_n_x_0 = log_k0_cat(x[0], 2, x0, delta_x);
                      
                      for (int y_val = min(ball_y); y_val <= max(ball_y); y_val++) {
                        log_k0_y_n = log_k0_drg(y_val, y, sigma, left);
                        log_k0_n_y = log_k0_drg(y, y_val, sigma, left);
                        
                        // Compute k_0(y,x | y_n,x_n) and k_0(y_n,x_n | y,x)
                        log_k0_yx_n = log_k0_y_n + log_k0_x_n_0 + log_k0_x_n_1 + log_k0_x_n_2 +
                          log_k0_x_n_3 + log_k0_x_n_4 + log_k0_x_n_5 + log_k0_x_n_6 + 
                          log_k0_x_n_7 + log_k0_x_n_8 + log_k0_x_n_9;
                        log_k0_n_yx = log_k0_n_y + log_k0_n_x_0 + log_k0_n_x_1 + log_k0_n_x_2 +
                          log_k0_n_x_3 + log_k0_n_x_4 + log_k0_n_x_5 + log_k0_n_x_6 + 
                          log_k0_n_x_7 + log_k0_n_x_8 + log_k0_n_x_9;
                        
                        // Recover p_{n-1}(y,x)
                        log_p_n_xy = log(p_n[idx_p]);
                        
                        // Compute \gamma_{n-1}(y,x, y_n,x_n)
                        log_gamma = log_p_n_xy + log_k0_n_yx - log(p_n_obs) - log_k0_yx_n;
                        log_gamma_2 = p_min(log_gamma);
                        
                        // Compute \gamma_{n-1}(y,x, y_n,x_n)*k_0(y,x | y_n,x_n)
                        out[idx_p] = exp(log_gamma_2 + log_k0_yx_n);
                        
                        // Update the index
                        idx_p++;
                      } } } } } } } } } } }
  return(out);
}

//[[Rcpp::export]]
double vec2idx_regr(IntegerVector y, int d, int lsy) {
  NumericVector idx(d);
  idx[0] = y[0];
  for(int i = 1; i < d; i++) {
    idx[i] = y[i] * lsy * pow(2, i-1);
  }
  return sum(idx) + 1;
}

//[[Rcpp::export]]
NumericVector mat2idx_regr(IntegerMatrix mat, int d, int lsy) {
  int N = mat.nrow();
  NumericVector out(N);
  for (int i = 0; i < N; i++) {
    out(i) = vec2idx_regr(mat(i,_), d, lsy);
  }
  return out;
}


