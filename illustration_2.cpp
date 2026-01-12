#include <random>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// Helper function to compute the standard normal CDF
double pnorm_cpp(double x, double mean, double sd) {
  return 0.5 * (1.0 + std::erf((x - mean) / (sd * std::sqrt(2.0))));
}

// Density of a rounded Gaussian distribution
// [[Rcpp::export]]
vec drg_cpp(ivec x, double mean, double sd) {
  int D = x.n_elem;
  vec out(D);
  
  for (int i = 0; i < D; i++) {
    out(i) = pnorm_cpp(x(i) + 0.5, mean, sd) - pnorm_cpp(x(i) - 0.5, mean, sd);
  }
  
  // Normalize
  double sum = accu(out);
  out = out / sum;
  
  // Floor small values
  out.elem(find(out <= 1e-250)).fill(1e-250);
  
  return out;
}

// Metropolis-Hastings kernel with rounded Gaussian
// [[Rcpp::export]]
vec mhk_rg_cpp(vec p_n, int y, double sd_rg, ivec support) {
  int n = p_n.n_elem;
  vec out(n);
  
  // Find index of y in support
  uvec idx_match = find(support == y);
  if (idx_match.n_elem == 0) {
    Rcpp::stop("Value y not found in support");
  }
  int idx_y = idx_match(0);
  
  // Compute acceptance probabilities
  vec prob_acc = clamp(p_n / p_n(idx_y), 0.0, 1.0);
  
  // Get rounded Gaussian kernel centered at y
  vec k0 = drg_cpp(support, static_cast<double>(y), sd_rg);
  
  // Compute transition probabilities
  out = prob_acc % k0;  // element-wise multiplication
  
  // Add rejection probability to diagonal
  double rej_prob = 1.0 - accu(out);
  out(idx_y) += rej_prob;
  
  return out;
}

// Stochastic approximation with Metropolis-Hastings kernel
// [[Rcpp::export]]
vec mclapply_function_cpp(int i, imat y, vec p_0, double sd_rg, 
                      ivec support, double lambda, 
                      double N_star, double alpha) {
  
  // Extract i-th column of y (R uses 1-based indexing)
  ivec y_perm = y.col(i - 1);
  
  // Initialize p_n with p_0
  vec p_n = p_0;
  
  // Iterate through y_perm
  int n_len = y_perm.n_elem;
  for (int n = 1; n <= n_len; n++) {
    // Compute MH kernel
    vec ker = mhk_rg_cpp(p_n, y_perm(n - 1), sd_rg, support);
    
    // Compute rho_n
    double rho_n = lambda + (1.0 - lambda) * std::exp(-(1.0 / N_star) * n);
    
    // Compute weight w_n
    double w_n = 1.0 / std::pow(alpha + n, rho_n);
    
    // Update p_n
    p_n = (1.0 - w_n) * p_n + w_n * ker;
  }
  
  return p_n;
}

// Helper function to sample from a discrete distribution
// [[Rcpp::export]]
int sample_discrete(ivec support, vec probs) {
  // Generate random number in [0,1)
  double u = R::runif(0.0, 1.0);
  
  // Compute cumulative probabilities
  double cumsum = 0.0;
  for (int i = 0; i < probs.n_elem; i++) {
    cumsum += probs(i);
    if (u < cumsum) {
      return support(i);
    }
  }
  
  // Return last element if rounding errors occur
  return support(probs.n_elem - 1);
}


// Stochastic approximation with predictive sampling
// [[Rcpp::export]]
vec mclapply_function_ps_cpp(int b, vec p_n, int m, int N, double sd_rg, 
                             double lambda, double N_star, double alpha, 
                             ivec support) {
  
  // Initialize prob with p_n
  vec prob = p_n;
  
  // Iterate from (m+1) to N
  for (int n = m + 1; n <= N; n++) {
    // Sample from support using current probabilities
    int y_past = sample_discrete(support, prob);
    
    // Compute MH kernel
    vec ker = mhk_rg_cpp(p_n, y_past, sd_rg, support);
    
    // Compute rho_n
    double rho_n = lambda + (1.0 - lambda) * std::exp(-(1.0 / N_star) * n);
    
    // Compute weight w_n
    double w_n = 1.0 / std::pow(alpha + n, rho_n);
    
    // Update prob
    prob = (1.0 - w_n) * prob + w_n * ker;
  }
  
  return prob;
}


// Stochastic approximation with pseudo-log-likelihood computation
// [[Rcpp::export]]
double mclapply_function_pq_cpp(int i, imat y, vec p_0, double sd_rg, 
                                ivec support, double lambda, 
                                double N_star, double alpha) {
  
  // Initialize pseudo-log-likelihood
  double pll = 0.0;
  
  // Extract i-th column of y (R uses 1-based indexing)
  ivec y_perm = y.col(i - 1);
  
  // Initialize p_n with p_0
  vec p_n = p_0;
  
  // Iterate through y_perm
  int n_len = y_perm.n_elem;
  for (int n = 1; n <= n_len; n++) {
    // Find index where support equals y_perm[n]
    uvec idx_match = find(support == y_perm(n - 1));
    if (idx_match.n_elem == 0) {
      Rcpp::stop("Value y_perm[n] not found in support");
    }
    int idx_y = idx_match(0);
    
    // Update pseudo-log-likelihood
    pll += std::log(p_n(idx_y));
    
    // Compute MH kernel
    vec ker = mhk_rg_cpp(p_n, y_perm(n - 1), sd_rg, support);
    
    // Compute rho_n
    double rho_n = lambda + (1.0 - lambda) * std::exp(-(1.0 / N_star) * n);
    
    // Compute weight w_n
    double w_n = 1.0 / std::pow(alpha + n, rho_n);
    
    // Update p_n
    p_n = (1.0 - w_n) * p_n + w_n * ker;
  }
  
  return pll;
}



// Stochastic approximation with predictive sampling and convergence tracking
// [[Rcpp::export]]
vec mclapply_function_ps_conv_val_cpp(int val, int b, vec p_n, int m, int N, double sd_rg, 
                                    double lambda, double N_star, double alpha, 
                                    ivec support) {
  
  // Initialize output vector
  vec out(N - m);
  out.fill(datum::nan);
  
  // Initialize prob with p_n
  vec prob = p_n;
  
  std::mt19937 generator(b);
  
  // Iterate from (m+1) to N
  int idx = 0;
  for (int n = m + 1; n <= N; n++) {
    
    // Sample from support using current probabilities
    int y_past = sample_discrete(support, prob);
    
    // Compute MH kernel
    vec ker = mhk_rg_cpp(p_n, y_past, sd_rg, support);
    
    // Compute rho_n
    double rho_n = lambda + (1.0 - lambda) * std::exp(-(1.0 / N_star) * n);
    
    // Compute weight w_n
    double w_n = 1.0 / std::pow(alpha + n, rho_n);
    
    // Update prob
    prob = (1.0 - w_n) * prob + w_n * ker;
    
    // Store sum of first 26 elements (0-indexed: 0 to 25)
    out(idx) = accu(prob.subvec(0, (val - 1)));
    idx++;
  }
  
  return out;
}


// Stochastic approximation with convergence tracking via mean absolute difference
// [[Rcpp::export]]
vec mclapply_ps_conv_cpp(vec p_n, int m, int N, double sd_rg, 
                         double lambda, double N_star, double alpha, 
                         ivec support) {
  
  // Initialize output vector
  vec out(N - m);
  out.fill(datum::nan);
  
  // Initialize prob with p_n
  vec prob = p_n;
  
  // Iterate from (m+1) to N
  for (int n = m + 1; n <= N; n++) {
    // Sample from support using current probabilities
    int y_past = sample_discrete(support, prob);
    
    // Compute MH kernel (fixed typo: y_pas -> y_past)
    vec ker = mhk_rg_cpp(p_n, y_past, sd_rg, support);
    
    // Compute rho_n
    double rho_n = lambda + (1.0 - lambda) * std::exp(-(1.0 / N_star) * n);
    
    // Compute weight w_n
    double w_n = 1.0 / std::pow(alpha + n, rho_n);
    
    // Update prob
    prob = (1.0 - w_n) * prob + w_n * ker;
    
    // Store mean absolute difference between real vectors prob and p_n
    out(n - m - 1) = sum(abs(prob - p_n));
    // out(n - m - 1) = sum(square(prob - p_n));
  }
  
  return out;
}


// Stochastic approximation with Metropolis-Hastings kernel
// [[Rcpp::export]]
vec mclapply_function_cons_cpp(int i, imat y, int sample_size, vec p_0, double sd_rg, 
                               ivec support, double lambda, 
                               double N_star, double alpha) {
  
  // Extract i-th column of y (R uses 1-based indexing)
  ivec y_perm = y.col(i - 1);
  
  // Initialize p_n with p_0
  vec p_n = p_0;
  
  // Iterate through y_perm
  int n_len = y_perm.n_elem;
  for (int n = 1; n <= n_len; n++) {
    // Compute MH kernel
    vec ker = mhk_rg_cpp(p_n, y_perm(n - 1), sd_rg, support);
    
    // Compute rho_n
    double rho_n = lambda + (1.0 - lambda) * std::exp(-(1.0 / N_star) * n);
    
    // Compute weight w_n
    double w_n = 1.0 / std::pow(alpha + sample_size + n, rho_n);
    
    // Update p_n
    p_n = (1.0 - w_n) * p_n + w_n * ker;
  }
  
  return p_n;
}