# Supporting information for "Nonparametric predictive inference for discrete data via Metropolis-adjusted Dirichlet sequences"

Code and data for the analysis contained in [[1]](#1).

- #### Regression
  - xx
  - xx

- #### Classification
  - xx
  - xx

- #### Application
  - [application.R](application.R) contains the analysis of Section 5. [predictive_resampling.R](predictive_resampling.R) contains the code for predictive resampling.
  - [Functions_Corvids.R](Functions_Corvids.R) and [Functions_Corvids.cpp](Functions_Corvids.cpp) contains the R and Cpp functions for MAD seqeunces employed in the application, respectively.
  - [post_processing_plots.R](post_processing_plots.R) and [post_processing_plots.R](post_processing_plots.R) contain the code for obtaining the related figures.
  - [corvus_2009.csv](corvus_2009.csv) contains the train test; [corvus_2010.csv](corvus_2010.csv) contains the test set; [sites_loc_09.csv](sites_loc_09.csv) contains the geographical information for each location.

- [Illustrative_Example.R](Illustrative_Example.R) replicate the illustrative example of Section 2.5.

## Demo
Demo which replicates the MAD sequence with adaptive weigths for the illustrative example of Section 2.5.

Ensure you have R (≥ 4.0) installed along with the mentioned packages.
The source files contained in XXXXXXX.
```r
library(parallel)
library(rmp)
library(dplyr)
library(tidyr)
library(ggplot2)
library(latex2exp)

# Load C++ functions
Rcpp::sourceCpp("illustration_2.cpp")

# Load R functions
source("functions_mad_github.R")
```

We simulate data from a mixture of four Poisson distributions
```r
# Sample size
nn <- 500
set.seed(1)

# Generate from mixture: 0.3*Pois(1) + 0.2*Pois(18) + 0.3*Pois(40) + 0.2*Pois(62)
y <- sample(c(1, 18, 40, 62), size = nn, replace = TRUE, prob = c(3, 2, 3, 2))
for (i in 1:nn) {
  y[i] <- rpois(1, lambda = y[i])
}

# Support for discrete distribution
supp <- 0:100

# True probability mass function
p_true <- 0.3*dpois(supp, lambda = 1) + 
          0.2*dpois(supp, lambda = 18) + 
          0.3*dpois(supp, lambda = 40) + 
          0.2*dpois(supp, lambda = 62)
```

Model Configuration
```r
# Base measure
p_0 <- rep(1/length(supp), length(supp))

# Concentration parameter
alpha <- 1

# Adaptive weighting parameters
N_star <- 500
lambda <- 0.75

# Number of permutations
PP <- 10

# Number of forward samples for predictive resampling
N <- 1e4

# Sample size of the martingale posterior sample
B <- 1000
```

The bandwidth parameter is optimized by maximizing the prequential log-likelihood
```r
opt <- nlminb(start = 1, lower = 1e-10, upper = Inf,
              objective = function(x) -preq_ll(x, y, p_0, alpha, lambda, N_star, supp, PP, cores = 1))

# Optimal bandwidth
sd <- opt$par
```

Posterior mean and predictive resampling
```r
# Posterior mean equals the predictive distribution at time n
pred <- pred_mad(y, p_0, alpha, lambda, N_star, sd, supp, PP, cores = 5)

# Predictive resampling for full posterior inference
ps <- pred_res(y, pred, alpha, lambda, N_star, supp, N, B, sd, cores = 8)

# Compute 95% credible intervals
ci_low <- apply(ps, 1, emp_hpd_lower)
ci_upp <- apply(ps, 1, emp_hpd_upper)
```

Visualization
```r
# Prepare data for plotting
df <- data.frame(supp = supp, pm = pred, ci_low = ci_low, ci_upp = ci_upp, title = "MAD sequence")

# Create plot
ggplot(df) + 
  facet_wrap(~ title) + 
  geom_segment(aes(x = supp, xend = supp, y = ci_low, yend = ci_upp, col = "z_ci"), 
               alpha = 0.15, size = 2) +
  geom_segment(aes(x = supp, xend = supp, y = 0, yend = pm, col = "pm"), 
               alpha = 1, size = .4) + 
  geom_point(aes(x = supp, y = p_true), col = "forestgreen", size = .5) +
  scale_color_manual(
    name = "", 
    values = c("pm" = "red", "z_ci" = "blue3"), 
    labels = c("pm" = TeX("$P_n$"), "z_ci" = "Posterior 0.95 CI")
  ) +
  ylim(c(0, 0.185)) + 
  xlim(c(0, 85)) +
  xlab("Support") + 
  ylab("Probability") +
  theme_light() +
  theme(
    legend.text.align = 0,
    strip.text = element_text(size = 12, colour = "black"),
    strip.background = element_rect(fill = "gray82"),
    panel.grid.major = element_line(size = 0.3, colour = "gray93"),
    panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.title = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5),
    legend.position = "top"
  )
```


## References
<a id="1">[1]</a> 
Agnoletto, D., Rigon, T., Dunson, D. B. (2025).
Nonparametric predictive inference for discrete data via Metropolis-adjusted Dirichlet sequences.
[ArXiv](https://arxiv.org/abs/2507.08629).
