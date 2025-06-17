
rm(list = ls())

load("ps_a.RData")

var_cc <- var_em <- var_cr <- var_ej <- rep(NA, NCOL(ps))
cov_cc_em <- cov_cc_cr <- cov_cc_ej <- cov_em_cr <- cov_em_ej <- cov_cr_ej <- rep(NA, NCOL(ps))
grid_supp <- as.matrix(expand.grid(0:50, 0:50))
for (s in 1:NCOL(ps)) {
  print(s)
  post <- SVT_SparseArray(dim = dim(p_n))
  post[idx_nz] <- ps[,s]
  
  pm_cc <- pm_em <- pm_cr <- pm_ej <- rep(NA, 51)
  for (j in 1:51) {
    pm_cc[j] <- sum(post[j,,,,,])
    pm_em[j] <- sum(post[,j,,,,])
    pm_cr[j] <- sum(post[,,j,,,])
    pm_ej[j] <- sum(post[,,,j,,])
  }
  pm_cc <- pm_cc / sum(pm_cc)
  pm_em <- pm_em / sum(pm_em)
  pm_cr <- pm_cr / sum(pm_cr)
  pm_ej <- pm_ej / sum(pm_ej)
  E_cc <- sum(0:50 * pm_cc)
  E_em <- sum(0:50 * pm_em)
  E_cr <- sum(0:50 * pm_cr)
  E_ej <- sum(0:50 * pm_ej)
  var_cc[s] <- sum(((0:50 - E_cc)^2) * pm_cc)
  var_em[s] <- sum(((0:50 - E_em)^2) * pm_em)
  var_cr[s] <- sum(((0:50 - E_cr)^2) * pm_cr)
  var_ej[s] <- sum(((0:50 - E_ej)^2) * pm_ej)
  
  pj_cc_em <- pj_cc_cr <- pj_cc_ej <- pj_em_cr <- pj_em_ej <- pj_cr_ej <- rep(NA, NROW(grid_supp))
  for (i in 1:NROW(grid_supp)) {
    pj_cc_em[i] <- sum(post[grid_supp[i,1]+1, grid_supp[i,2]+1,,,,])
    pj_cc_cr[i] <- sum(post[grid_supp[i,1]+1,, grid_supp[i,2]+1,,,])
    pj_cc_ej[i] <- sum(post[grid_supp[i,1]+1,,, grid_supp[i,2]+1,,])
    pj_em_cr[i] <- sum(post[, grid_supp[i,1]+1, grid_supp[i,2]+1,,,])
    pj_em_ej[i] <- sum(post[, grid_supp[i,1]+1,, grid_supp[i,2]+1,,])
    pj_cr_ej[i] <- sum(post[,, grid_supp[i,1]+1, grid_supp[i,2]+1,,])
  }
  pj_cc_em <- pj_cc_em / sum(pj_cc_em)
  pj_cc_cr <- pj_cc_cr / sum(pj_cc_cr)
  pj_cc_ej <- pj_cc_ej / sum(pj_cc_ej)
  pj_em_cr <- pj_em_cr / sum(pj_em_cr)
  pj_em_ej <- pj_em_ej / sum(pj_em_ej)
  pj_cr_ej <- pj_cr_ej / sum(pj_cr_ej)
  j_prod_cc_em <- sum(grid_supp[,1] * grid_supp[,2] * pj_cc_em)
  j_prod_cc_cr <- sum(grid_supp[,1] * grid_supp[,2] * pj_cc_cr)
  j_prod_cc_ej <- sum(grid_supp[,1] * grid_supp[,2] * pj_cc_ej)
  j_prod_em_cr <- sum(grid_supp[,1] * grid_supp[,2] * pj_em_cr)
  j_prod_em_ej <- sum(grid_supp[,1] * grid_supp[,2] * pj_em_ej)
  j_prod_cr_ej <- sum(grid_supp[,1] * grid_supp[,2] * pj_cr_ej)
  cov_cc_em[s] <- j_prod_cc_em - E_cc*E_em
  cov_cc_cr[s] <- j_prod_cc_cr - E_cc*E_cr
  cov_cc_ej[s] <- j_prod_cc_ej - E_cc*E_ej
  cov_em_cr[s] <- j_prod_em_cr - E_em*E_cr
  cov_em_ej[s] <- j_prod_em_ej - E_em*E_ej
  cov_cr_ej[s] <- j_prod_cr_ej - E_cr*E_ej
}
save(var_cc, var_em, var_cr, var_ej, file = file("ps_variance.RData", blocking = TRUE))
save(cov_cc_em, cov_cc_cr, cov_cc_ej, cov_em_cr, cov_em_ej, cov_cr_ej, file = file("ps_covariance.RData", blocking = TRUE))

cor_cc_em <- cor_cc_cr <- cor_cc_ej <- cor_em_cr <- cor_em_ej <- cor_cr_ej <- rep(NA, NCOL(ps))
for (i in 1:NCOL(ps)) {
  cor_cc_em[i] <- cov_cc_em[i] / sqrt(var_cc[i] * var_em[i])
  cor_cc_cr[i] <- cov_cc_cr[i] / sqrt(var_cc[i] * var_cr[i])
  cor_cc_ej[i] <- cov_cc_ej[i] / sqrt(var_cc[i] * var_ej[i])
  cor_em_cr[i] <- cov_em_cr[i] / sqrt(var_em[i] * var_cr[i])
  cor_em_ej[i] <- cov_em_ej[i] / sqrt(var_em[i] * var_ej[i])
  cor_cr_ej[i] <- cov_cr_ej[i] / sqrt(var_cr[i] * var_ej[i])
}
save(cor_cc_em, cor_cc_cr, cor_cc_ej, cor_em_cr, cor_em_ej, cor_cr_ej, file = file("ps_correlation.RData", blocking = TRUE))



