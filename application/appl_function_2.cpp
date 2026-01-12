#include <RcppDist.h>
#include <RcppArmadillo.h>
#include <unordered_map>
#include <unordered_set>

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
  return std::min(0.0, x);
}

// [[Rcpp::export]]
double log_k0_cat(double x, int size, double x_obs, double delta) {
  double out = (x == x_obs) ? (1.0 - delta) : (delta / (size - 1));
  return log(out);
}

//[[Rcpp::export]]
double idx_p(arma::ivec y, int d, int lsy, int lsh, int lst) {
  double sum_y_idx = 0;
  int multiplier = 1;
  for (int i = 0; i < d; i++) {
    sum_y_idx += y(i) * multiplier;
    multiplier *= lsy;
  }
  double t_idx = y(d) * multiplier;
  double h_idx = (y(d+1) - 1.0) * lst * multiplier;
  return sum_y_idx + t_idx + h_idx + 1;
}

//[[Rcpp::export]]
vec mat2idx(arma::imat mat, int d, int lsy, int lst, int lsh) {
  int N = mat.n_rows;
  vec out(N);
  for (int i = 0; i < N; i++) {
    out(i) = idx_p(mat.row(i).t(), d, lsy, lsh, lst);
  }
  return out;
}

// [[Rcpp::export]]
vec idx2vec(int idx, int D, int lsy, int lst, int lsh) {
  vec out(D+2);
  int num = idx - 1;
  for (int i = 0; i < D; i++) {
    out(i) = num % lsy;
    num = num / lsy;
  }
  out(D) = num % lst;
  num = num / lst;
  out(D+1) = (num % lsh) + 1;
  return out;
}

// [[Rcpp::export]]
arma::ivec k0_int(double x, double marg, arma::ivec support) {
  int min_sup = support.min();
  int max_sup = support.max();
  
  int dx = std::min(max_sup, (int)std::ceil(x + marg));
  int sx = std::max(min_sup, (int)std::floor(x - marg));
  
  int len = dx - sx + 1;
  arma::ivec out(len);
  for (int i = 0; i < len; i++) {
    out(i) = sx + i;
  }
  return out;
}

// Optimized expand_grid that directly computes indices
void expand_grid_recursive(const std::vector<arma::ivec>& ball, imat& result, 
                           ivec& current, int depth, int& row_idx) {
  if (depth == ball.size()) {
    result.row(row_idx++) = current.t();
    return;
  }
  
  for (uword i = 0; i < ball[depth].n_elem; i++) {
    current(depth) = ball[depth](i);
    expand_grid_recursive(ball, result, current, depth + 1, row_idx);
  }
}

imat expand_grid_ball(const std::vector<arma::ivec>& ball) {
  int D = ball.size();
  
  // Calculate total number of combinations
  int n_rows = 1;
  for (int i = 0; i < D; i++) {
    n_rows *= ball[i].n_elem;
  }
  
  imat result(n_rows, D);
  ivec current(D);
  int row_idx = 0;
  expand_grid_recursive(ball, result, current, 0, row_idx);
  
  return result;
}

// Optimized sampling function using cumulative probabilities
int weighted_sample(const vec& probs) {
  double u = R::runif(0, 1);
  double cumsum = 0.0;
  for (uword i = 0; i < probs.n_elem; i++) {
    cumsum += probs(i);
    if (u <= cumsum) return i;
  }
  return probs.n_elem - 1;
}

// [[Rcpp::export]]
vec mhk_4_ps(vec p_n_vec, double p_n_obs, int idx_tr_obs,
             arma::ivec idx_tr, vec sigma, double delta_temp,
             double delta_hab, int lsy, int lst, int lsh, double left = -0.5) {
  
  int D = 4;
  int n = p_n_vec.n_elem;
  vec out(n);
  
  vec val_obs = idx2vec(idx_tr_obs, D, lsy, lst, lsh);
  double y_0 = val_obs(0);
  double y_1 = val_obs(1);
  double y_2 = val_obs(2);
  double y_3 = val_obs(3);
  double temp = val_obs(4);
  double hab = val_obs(5);
  
  for (int i = 0; i < n; i++) {
    vec val = idx2vec(idx_tr(i), D, lsy, lst, lsh);
    double y_0_val = val(0);
    double y_1_val = val(1);
    double y_2_val = val(2);
    double y_3_val = val(3);
    double temp_val = val(4);
    double hab_val = val(5);
    
    double log_k0_y_n_0 = log_drg(y_0_val, y_0, sigma(0), left);
    double log_k0_n_y_0 = log_drg(y_0, y_0_val, sigma(0), left);
    double log_k0_y_n_1 = log_drg(y_1_val, y_1, sigma(1), left);
    double log_k0_n_y_1 = log_drg(y_1, y_1_val, sigma(1), left);
    double log_k0_y_n_2 = log_drg(y_2_val, y_2, sigma(2), left);
    double log_k0_n_y_2 = log_drg(y_2, y_2_val, sigma(2), left);
    double log_k0_y_n_3 = log_drg(y_3_val, y_3, sigma(3), left);
    double log_k0_n_y_3 = log_drg(y_3, y_3_val, sigma(3), left);
    double log_k0_temp_n = log_k0_cat(temp_val, lst, temp, delta_temp);
    double log_k0_n_temp = log_k0_cat(temp, lst, temp_val, delta_temp);
    double log_k0_hab_n = log_k0_cat(hab_val, lsh, hab, delta_hab);
    double log_k0_n_hab = log_k0_cat(hab, lsh, hab_val, delta_hab);
    
    double log_k0_yx_n = log_k0_y_n_0 + log_k0_y_n_1 + log_k0_y_n_2 + log_k0_y_n_3 +
      log_k0_temp_n + log_k0_hab_n;
    double log_k0_n_yx = log_k0_n_y_0 + log_k0_n_y_1 + log_k0_n_y_2 + log_k0_n_y_3 +
      log_k0_n_temp + log_k0_n_hab;
    
    double log_p_n_xy = log(p_n_vec(i));
    double log_gamma = log_p_n_xy + log_k0_n_yx - log(p_n_obs) - log_k0_yx_n;
    double log_gamma_2 = p_min(log_gamma);
    
    out(i) = exp(log_gamma_2 + log_k0_yx_n);
  }
  
  return out;
}

// [[Rcpp::export]]
vec pred_res_4_FB_one_conv(vec p_n_vec, double p_n_az, arma::ivec idx_nz, 
                           int sample_size, int N, double alpha, double lambda, 
                           double N_star, vec sigma, vec delta, 
                           arma::ivec support_y, arma::ivec support_temp, 
                           arma::ivec support_hab, vec interval) {
  
  int D = 4;
  int lsy = support_y.n_elem;
  int lst = support_temp.n_elem;
  int lsh = support_hab.n_elem;
  double delta_temp = delta(0);
  double delta_hab = delta(1);
  
  vec out(N);
  vec p_n_vec_start = p_n_vec;
  
  // Create hash map for fast lookup
  std::unordered_map<int, int> idx_nz_map;
  for (uword i = 0; i < idx_nz.n_elem; i++) {
    idx_nz_map[idx_nz(i)] = i;
  }
  
  for (int n = 0; n < N; n++) {
    
    if ((n + 1) % 100 == 0) {
      Rcpp::Rcout << "Iteration: " << n + 1 << std::endl;
      Rcpp::checkUserInterrupt();
    }
    
    // Sample new observation
    int y_tr_obs_idx = weighted_sample(p_n_vec);
    int y_tr_obs = idx_nz(y_tr_obs_idx);
    double p_n_obs = p_n_vec(y_tr_obs_idx);
    
    // Recover y_n
    vec rec_val = idx2vec(y_tr_obs, D, lsy, lst, lsh);
    
    // Ball around y_obs
    std::vector<arma::ivec> ball(D + 2);
    for (int k = 0; k < D; k++) {
      ball[k] = k0_int(rec_val(k), interval(k), support_y);
    }
    ball[D] = support_temp;
    ball[D+1] = support_hab;
    
    // Create Cartesian product
    imat ball_mat = expand_grid_ball(ball);
    vec ball_tr = mat2idx(ball_mat, D, lsy, lst, lsh);
    
    // Find indices in idx_nz using hash map
    std::vector<int> idx_ball_nz_vec;
    std::vector<int> ball_tr_matched;
    idx_ball_nz_vec.reserve(ball_tr.n_elem);
    ball_tr_matched.reserve(ball_tr.n_elem);
    
    for (uword i = 0; i < ball_tr.n_elem; i++) {
      int val = (int)ball_tr(i);
      auto it = idx_nz_map.find(val);
      if (it != idx_nz_map.end()) {
        idx_ball_nz_vec.push_back(it->second);
        ball_tr_matched.push_back(val);
      }
    }
    
    if (idx_ball_nz_vec.empty()) {
      out(n) = mean(abs(p_n_vec - p_n_vec_start));
      continue;
    }
    
    // Convert to Armadillo vectors
    uvec idx_ball_nz = conv_to<uvec>::from(idx_ball_nz_vec);
    arma::ivec ball_tr_int = conv_to<arma::ivec>::from(ball_tr_matched);
    
    // Compute gamma * k_0
    vec ker = mhk_4_ps(p_n_vec.elem(idx_ball_nz), p_n_obs, y_tr_obs, 
                       ball_tr_int, sigma, delta_temp, delta_hab, lsy, lst, lsh);
    
    // Update p_{n-1}
    int n_cum = sample_size + n + 1;
    double rho_n = lambda + (1.0 - lambda) * std::exp(-(1.0/N_star) * n_cum);
    double w_n = 1.0 / std::pow(alpha + n_cum, rho_n);
    
    p_n_vec *= (1.0 - w_n);
    p_n_vec.elem(idx_ball_nz) += w_n * ker;
    
    // Add sum in y_n
    double ker_sum = sum(ker);
    p_n_vec(y_tr_obs_idx) += w_n * (1.0 - ker_sum);
    
    // Update other values
    p_n_az *= (1.0 - w_n);
    
    // Calculate output
    out(n) = sum(abs(p_n_vec - p_n_vec_start));
  }
  
  return out;
}


// [[Rcpp::export]]
vec pred_res_4_FB_one(vec p_n_vec, double p_n_az, arma::ivec idx_nz, 
                      int sample_size, int N, double alpha, double lambda, 
                      double N_star, vec sigma, vec delta, 
                      arma::ivec support_y, arma::ivec support_temp, 
                      arma::ivec support_hab, vec interval) {
  
  int D = 4;
  int lsy = support_y.n_elem;
  int lst = support_temp.n_elem;
  int lsh = support_hab.n_elem;
  double delta_temp = delta(0);
  double delta_hab = delta(1);
  
  // Create hash map for fast lookup
  std::unordered_map<int, int> idx_nz_map;
  for (uword i = 0; i < idx_nz.n_elem; i++) {
    idx_nz_map[idx_nz(i)] = i;
  }
  
  for (int n = 0; n < N; n++) {
    
    if ((n + 1) % 100 == 0) {
      Rcpp::Rcout << "Iteration: " << n + 1 << std::endl;
      Rcpp::checkUserInterrupt();
    }
    
    // Sample new observation
    int y_tr_obs_idx = weighted_sample(p_n_vec);
    int y_tr_obs = idx_nz(y_tr_obs_idx);
    double p_n_obs = p_n_vec(y_tr_obs_idx);
    
    // Recover y_n
    vec rec_val = idx2vec(y_tr_obs, D, lsy, lst, lsh);
    
    // Ball around y_obs
    std::vector<arma::ivec> ball(D + 2);
    for (int k = 0; k < D; k++) {
      ball[k] = k0_int(rec_val(k), interval(k), support_y);
    }
    ball[D] = support_temp;
    ball[D+1] = support_hab;
    
    // Create Cartesian product
    imat ball_mat = expand_grid_ball(ball);
    vec ball_tr = mat2idx(ball_mat, D, lsy, lst, lsh);
    
    // Find indices in idx_nz using hash map
    std::vector<int> idx_ball_nz_vec;
    std::vector<int> ball_tr_matched;
    idx_ball_nz_vec.reserve(ball_tr.n_elem);
    ball_tr_matched.reserve(ball_tr.n_elem);
    
    for (uword i = 0; i < ball_tr.n_elem; i++) {
      int val = (int)ball_tr(i);
      auto it = idx_nz_map.find(val);
      if (it != idx_nz_map.end()) {
        idx_ball_nz_vec.push_back(it->second);
        ball_tr_matched.push_back(val);
      }
    }
    
    // Convert to Armadillo vectors
    uvec idx_ball_nz = conv_to<uvec>::from(idx_ball_nz_vec);
    arma::ivec ball_tr_int = conv_to<arma::ivec>::from(ball_tr_matched);
    
    // Compute gamma * k_0
    vec ker = mhk_4_ps(p_n_vec.elem(idx_ball_nz), p_n_obs, y_tr_obs, 
                       ball_tr_int, sigma, delta_temp, delta_hab, lsy, lst, lsh);
    
    // Update p_{n-1}
    int n_cum = sample_size + n + 1;
    double rho_n = lambda + (1.0 - lambda) * std::exp(-(1.0/N_star) * n_cum);
    double w_n = 1.0 / std::pow(alpha + n_cum, rho_n);
    
    p_n_vec *= (1.0 - w_n);
    p_n_vec.elem(idx_ball_nz) += w_n * ker;
    
    // Add sum in y_n
    double ker_sum = sum(ker);
    p_n_vec(y_tr_obs_idx) += w_n * (1.0 - ker_sum);
    
    // Update other values
    p_n_az *= (1.0 - w_n);
  }
  
  return p_n_vec;
}