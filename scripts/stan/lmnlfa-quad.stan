data {
  // dimensions
  int<lower=1> nobs;              // total person-occasion-item observations
  int<lower=1> p;                 // number of items
  int<lower=1> ni;                // number of individuals
  int<lower=1> d;                 // number of time points (waves)

  // indices
  array[nobs] int<lower=1, upper=ni> person;
  array[nobs] int<lower=1, upper=p>  itm;
  array[nobs] int<lower=1, upper=d>  time;

  // time metric (centered age and centered age^2, computed in R)
  vector[nobs] age_c;
  vector[nobs] age2_c;

  // responses
  array[nobs] int y;                       // holds both 0/1 and 1..k responses
  array[p]    int<lower=0, upper=1> is_binary; // 1 if item is 0/1, else 0

  // for ordinal items:
  array[p] int<lower=2> k_item;  // number of categories per item (set 2 for binary, harmless)
  int<lower=1> k_max;            // max categories across items (e.g. 4)

  // predictors
  int<lower=0> nfpreds;           // # time-invariant predictors
  int<lower=0> ntvpreds;          // # time-varying predictors
  matrix[ni,   nfpreds] xf_person;  // person-level covariates
  matrix[nobs, nfpreds] xf;         // time-invariant covariates (baseline age, etc.)
  matrix[nobs, ntvpreds] xtv;       // time-varying covariates

  // DIF pattern (col 1 = time-varying DIF by age, col 2 = invariant DIF)
  matrix[p, 2] ldf;
  int<lower=0> mtv;   // count of items flagged for time-varying DIF
  int<lower=0> mf;    // count of items flagged for invariant DIF

  // prior scales
  real<lower=0> sigma_l;
  real<lower=0> sigma_nu;
  real<lower=0> sigma_cor;
  real<lower=0> sigma_f;
  real<lower=0> sigma_di;
}

parameters {
  // baseline loadings/intercepts
  vector<lower=0>[p] lp;
  vector[p] np;

  // DIF params (mapped to item-length vectors via ldf)
  vector[mf]  l_diff;    // loading DIF invariant
  vector[mtv] l_diftv;   // loading DIF time-varying
  vector[mf]  n_diff;    // intercept DIF invariant
  vector[mtv] n_diftv;   // intercept DIF time-varying

  // growth factor means (intercept mean fixed to 0 for identification)
  real mu_slp;
  real mu_quad;

  // growth factor SDs
  real<lower=0> phi_int;
  real<lower=0> phi_slp;
  real<lower=0> phi_quad;
  real<lower=0> eti_sd;

  // impact of time-invariant predictors on growth factor means/SDs
  matrix[3, nfpreds] b_mu;
  matrix[3, nfpreds] b_phi;

  // Cholesky factor of 3×3 correlation matrix among growth factors
  cholesky_factor_corr[3] L_Omega;

  // non-centered random effects
  matrix[3, ni] fac_dist;     // std-normal draws for growth random effects
  matrix[d, ni] fac_eti_raw;  // time-specific residuals (before scaling)

  // ordinal thresholds (k_max-1 per item; ignored for binary items)
  array[p] ordered[k_max - 1] tau;
}

transformed parameters {
  // expand DIF vectors to item length (0 for non-flagged items)
  vector[p] ldiff  = rep_vector(0, p);
  vector[p] ldiftv = rep_vector(0, p);
  vector[p] ndiff  = rep_vector(0, p);
  vector[p] ndiftv = rep_vector(0, p);

  vector<lower=0>[3] phi_eta;
  phi_eta[1] = phi_int;
  phi_eta[2] = phi_slp;
  phi_eta[3] = phi_quad;

  matrix[d, ni] fac_eti = fac_eti_raw * eti_sd;

  {
    int tmp;

    tmp = 0;
    for (i in 1:p) {
      if (ldf[i, 2] == 1) {
        tmp += 1;
        ldiff[i] = l_diff[tmp];
      }
    }

    tmp = 0;
    for (i in 1:p) {
      if (ldf[i, 1] == 1) {
        tmp += 1;
        ldiftv[i] = l_diftv[tmp];
      }
    }

    tmp = 0;
    for (i in 1:p) {
      if (ldf[i, 2] == 1) {
        tmp += 1;
        ndiff[i] = n_diff[tmp];
      }
    }

    tmp = 0;
    for (i in 1:p) {
      if (ldf[i, 1] == 1) {
        tmp += 1;
        ndiftv[i] = n_diftv[tmp];
      }
    }
  }
}

model {
  // --- local declarations must come first ---
  matrix[3, ni] fac_gr;

  // --- priors ---
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

  phi_int  ~ normal(0, sigma_f);
  phi_slp  ~ normal(0, sigma_f);
  phi_quad ~ normal(0, sigma_f);

  L_Omega ~ lkj_corr_cholesky(sigma_cor);
  to_vector(fac_dist)    ~ normal(0, 1);
  to_vector(fac_eti_raw) ~ normal(0, 1);
  for (it in 1:p) {
    tau[it] ~ normal(0, 1.5);
  }

  // --- growth factors per person (non-centered parameterization) ---
  for (k in 1:ni) {
    vector[3] mu_eta;
    vector[3] sd_eta;
    mu_eta[1] = 0;
    mu_eta[2] = mu_slp;
    mu_eta[3] = mu_quad;

    sd_eta = phi_eta .* exp(b_phi * (xf_person[k, ]'));
    fac_gr[, k] = mu_eta
                + b_mu * (xf_person[k, ]')
                + diag_pre_multiply(sd_eta, L_Omega) * fac_dist[, k];
  }

  // --- likelihood ---
  for (j in 1:nobs) {
    int ti = time[j];
    int it = itm[j];
    int pe = person[j];

    real eta_j = fac_gr[1, pe]
               + fac_gr[2, pe] * age_c[j]
               + fac_gr[3, pe] * age2_c[j]
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
  // posterior correlation matrix among growth factors
  matrix[3, 3] Omega = multiply_lower_tri_self_transpose(L_Omega);
}
