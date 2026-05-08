data {
  // dimensions
  int<lower=1> nobs;
  int<lower=1> p;
  int<lower=1> ni;
  int<lower=1> d;

  // indices
  array[nobs] int<lower=1, upper=ni> person;
  array[nobs] int<lower=1, upper=p>  itm;
  array[nobs] int<lower=1, upper=d>  time;

  // time metric
  vector[nobs] age_c;
  vector[nobs] age2_c;

  // responses
  array[nobs] int y;
  array[p]    int<lower=0, upper=1> is_binary;

  // ordinal structure
  array[p] int<lower=2> k_item;
  int<lower=1> k_max;

  // predictors
  int<lower=0> nfpreds;
  int<lower=0> ntvpreds;
  matrix[ni,   nfpreds] xf_person;
  matrix[nobs, nfpreds] xf;
  matrix[nobs, ntvpreds] xtv;

  // DIF pattern (col 1 = time-varying, col 2 = invariant)
  matrix[p, 2] ldf;
  int<lower=0> mtv;
  int<lower=0> mf;

  // prior scales
  real<lower=0> sigma_l;
  real<lower=0> sigma_nu;
  real<lower=0> sigma_cor;
  real<lower=0> sigma_f;
  real<lower=0> sigma_di;
}

parameters {
  // item parameters
  vector<lower=0>[p] lp;
  vector[p] np;

  // DIF parameters
  vector[mf]  l_diff;
  vector[mtv] l_diftv;
  vector[mf]  n_diff;
  vector[mtv] n_diftv;

  // fixed growth means (intercept mean = 0 for identification)
  real mu_slp;
  real mu_quad;   // fixed quadratic; individuals share same curvature

  // random intercept + slope SDs (2-factor; no random quadratic)
  real<lower=0> phi_int;
  real<lower=0> phi_slp;
  real<lower=0> eti_sd;

  // covariate impact on growth factor means / SDs (2 growth factors)
  matrix[2, nfpreds] b_mu;
  matrix[2, nfpreds] b_phi;

  // Cholesky of 2×2 correlation between intercept and slope
  cholesky_factor_corr[2] L_Omega;

  // non-centered random effects (2 × ni)
  matrix[2, ni] fac_dist;
  matrix[d, ni] fac_eti_raw;

  // ordinal thresholds
  array[p] ordered[k_max - 1] tau;
}

transformed parameters {
  vector[p] ldiff  = rep_vector(0, p);
  vector[p] ldiftv = rep_vector(0, p);
  vector[p] ndiff  = rep_vector(0, p);
  vector[p] ndiftv = rep_vector(0, p);

  vector<lower=0>[2] phi_eta;
  phi_eta[1] = phi_int;
  phi_eta[2] = phi_slp;

  matrix[d, ni] fac_eti = fac_eti_raw * eti_sd;

  {
    int tmp;

    tmp = 0;
    for (i in 1:p) {
      if (ldf[i, 2] == 1) { tmp += 1; ldiff[i] = l_diff[tmp]; }
    }
    tmp = 0;
    for (i in 1:p) {
      if (ldf[i, 1] == 1) { tmp += 1; ldiftv[i] = l_diftv[tmp]; }
    }
    tmp = 0;
    for (i in 1:p) {
      if (ldf[i, 2] == 1) { tmp += 1; ndiff[i] = n_diff[tmp]; }
    }
    tmp = 0;
    for (i in 1:p) {
      if (ldf[i, 1] == 1) { tmp += 1; ndiftv[i] = n_diftv[tmp]; }
    }
  }
}

model {
  // local declarations first
  matrix[2, ni] fac_gr;

  // priors
  lp      ~ normal(0, sigma_l);
  np      ~ normal(0, sigma_nu);
  eti_sd  ~ normal(0, sigma_f);

  l_diff  ~ normal(0, sigma_di);
  l_diftv ~ normal(0, sigma_di);
  n_diff  ~ normal(0, sigma_di);
  n_diftv ~ normal(0, sigma_di);

  to_vector(b_mu)  ~ normal(0, sigma_f);
  to_vector(b_phi) ~ normal(0, sigma_f);

  mu_slp  ~ normal(0, sigma_f);
  mu_quad ~ normal(0, sigma_f);

  phi_int ~ normal(0, sigma_f);
  phi_slp ~ normal(0, sigma_f);

  L_Omega ~ lkj_corr_cholesky(sigma_cor);
  to_vector(fac_dist)    ~ normal(0, 1);
  to_vector(fac_eti_raw) ~ normal(0, 1);
  for (it in 1:p) {
    tau[it] ~ normal(0, 1.5);
  }

  // growth factors: random intercept + slope, fixed quadratic mean
  for (k in 1:ni) {
    vector[2] mu_eta;
    vector[2] sd_eta;
    mu_eta[1] = 0;
    mu_eta[2] = mu_slp;

    sd_eta = phi_eta .* exp(b_phi * (xf_person[k, ]'));
    fac_gr[, k] = mu_eta
                + b_mu * (xf_person[k, ]')
                + diag_pre_multiply(sd_eta, L_Omega) * fac_dist[, k];
  }

  // likelihood
  for (j in 1:nobs) {
    int ti = time[j];
    int it = itm[j];
    int pe = person[j];

    // individual trajectory: random intercept + slope, shared quadratic
    real eta_j = fac_gr[1, pe]
               + fac_gr[2, pe] * age_c[j]
               + mu_quad        * age2_c[j]
               + fac_eti[ti, pe];

    real nu  = np[it]
             + ndiff[it]  * (nfpreds > 0 ? xf[j, 1] : 0.0)
             + ndiftv[it] * age_c[j];
    real lam = lp[it] * exp(ldiff[it]  * (nfpreds > 0 ? xf[j, 1] : 0.0)
                           + ldiftv[it] * age_c[j]);

    if (is_binary[it] == 1) {
      y[j] ~ bernoulli_logit(nu + lam * eta_j);
    } else {
      y[j] ~ ordered_logistic(nu + lam * eta_j, head(tau[it], k_item[it] - 1));
    }
  }
}

generated quantities {
  // recover correlation between random intercept and slope
  matrix[2, 2] Omega = multiply_lower_tri_self_transpose(L_Omega);
}
