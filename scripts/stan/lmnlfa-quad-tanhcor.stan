// lmnlfa-quad-tanhcor.stan
// Identical to lmnlfa-quad.stan EXCEPT for how the intercept/slope
// correlation is parameterized.
//
// lmnlfa-quad.stan uses cholesky_factor_corr[2] + lkj_corr_cholesky, which
// can propose points where the Cholesky factor's diagonal hits exactly 0
// (rho -> +-1), triggering frequent "Random variable is 0, but must be
// positive" rejections and (per the growth-only smoke test) severe
// ill-conditioning even with no DIF/impact in the model.
//
// Here the single correlation is instead an unconstrained real, mapped
// through tanh() to (-1, 1) -- it can get arbitrarily close to +-1 but can
// never reach it exactly, so the Cholesky factor of the resulting 2x2
// correlation matrix never degenerates. This mirrors the approach already
// prototyped in lmnlfa-linear.stan (z_cor / tanh), applied here to the
// growth-only isolation model (no DIF, no covariate impact) with the
// current quad-model's data/parameter layout otherwise unchanged, so
// results are directly comparable to lmnlfa-quad.stan and the existing
// R helpers (which key off "Omega[1,2]") work unmodified.

data {
  int<lower=1> nobs;
  int<lower=2> p;  // >=2 so item 1 can serve as the scale-identifying marker
  int<lower=1> ni;
  int<lower=1> d;

  array[nobs] int<lower=1, upper=ni> person;
  array[nobs] int<lower=1, upper=p>  itm;
  array[nobs] int<lower=1, upper=d>  time;

  vector[nobs] age_c;
  vector[nobs] age2_c;

  array[nobs] int y;
  array[p]    int<lower=0, upper=1> is_binary;

  array[p] int<lower=2> k_item;
  int<lower=1> k_max;

  int<lower=0> nfpreds;
  int<lower=0> ntvpreds;
  matrix[ni,   nfpreds] xf_person;
  matrix[nobs, nfpreds] xf;
  matrix[nobs, ntvpreds] xtv;

  matrix[p, 2] ldf;
  int<lower=0> mtv;
  int<lower=0> mf;

  real<lower=0> sigma_l;
  real<lower=0> sigma_nu;
  real<lower=0> sigma_cor;
  real<lower=0> sigma_f;
  real<lower=0> sigma_di;
}

parameters {
  // item 1's loading is fixed to 1 (see transformed parameters) as a marker
  // item for scale identification -- lp and the growth-factor SDs
  // (phi_int, phi_slp) are otherwise only jointly identified through their
  // product, which produces a strong non-identified ridge (confirmed in the
  // static model: cor(lp, phi0) ~ -0.92). With two SDs riding the same
  // ridge against a shared loading vector, plus their correlation, this is
  // the most likely root cause of the extreme threshold values and 100%
  // max-treedepth saturation seen in the growth-only smoke tests.
  vector<lower=0>[p - 1] lp_free;
  vector[p] np;

  vector[mf]  l_diff;
  vector[mtv] l_diftv;
  vector[mf]  n_diff;
  vector[mtv] n_diftv;

  real mu_slp;

  real<lower=0> phi_int;
  real<lower=0> phi_slp;

  matrix[2, nfpreds] b_mu;
  matrix[2, nfpreds] b_phi;

  // unconstrained intercept/slope correlation (replaces cholesky_factor_corr)
  real z_cor;

  matrix[2, ni] fac_dist;

  // occasion-specific residual ("wobble"): lets a person's true trajectory
  // deviate from the perfectly linear-in-age growth curve at each wave,
  // rather than forcing every observation to sit exactly on that line.
  // Present in the source paper's model (fac_eti) and in an earlier draft
  // of this model (lmnlfa-linear.stan's fac_eti_raw) but missing here until
  // now -- `ti` was already being computed in the likelihood loop below but
  // never used, which is the tell that this was dropped rather than never
  // planned.
  matrix[d, ni] fac_eti_raw;
  real<lower=0> eti_sd;

  array[p] ordered[k_max - 1] tau;
}

transformed parameters {
  vector[p] lp;
  lp[1] = 1;
  lp[2:p] = lp_free;

  vector[p] ldiff  = rep_vector(0, p);
  vector[p] ldiftv = rep_vector(0, p);
  vector[p] ndiff  = rep_vector(0, p);
  vector[p] ndiftv = rep_vector(0, p);

  vector<lower=0>[2] phi_eta;
  phi_eta[1] = phi_int;
  phi_eta[2] = phi_slp;

  // rho is strictly in (-1, 1): the Cholesky factor below never degenerates
  real<lower=-1, upper=1> rho = tanh(z_cor);
  matrix[2, 2] L_Omega;
  L_Omega[1, 1] = 1;
  L_Omega[2, 1] = rho;
  L_Omega[1, 2] = 0;
  L_Omega[2, 2] = sqrt(1 - square(rho));

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
  matrix[2, ni] fac_gr;

  // lambda ~ N+(1, 1.5) matched to the source paper (lp_free's <lower=0>
  // constraint makes this N+(1, sigma_l)).
  //
  // nu is NOT matched to the paper's N(-2, 1.5): that assumes a rare,
  // hard-to-endorse binary item. These PDS items are ordinal (1-4) and
  // empirically skew the opposite way -- pooled across waves, all 8 items
  // are ABOVE their lowest category 80-95% of the time (implied logits
  // +1.4 to +2.9, both reporters). A -2 center would actively point the
  // prior in the wrong direction, so keep nu centered at 0.
  lp_free ~ normal(1, sigma_l);
  np      ~ normal(0, sigma_nu);

  l_diff  ~ normal(0, sigma_di);
  l_diftv ~ normal(0, sigma_di);
  n_diff  ~ normal(0, sigma_di);
  n_diftv ~ normal(0, sigma_di);

  to_vector(b_mu)  ~ normal(0, sigma_f);
  to_vector(b_phi) ~ normal(0, sigma_f);

  mu_slp  ~ normal(0, sigma_f);

  phi_int ~ normal(0, sigma_f);
  phi_slp ~ normal(0, sigma_f);

  z_cor ~ normal(0, sigma_cor);
  to_vector(fac_dist)    ~ normal(0, 1);

  eti_sd ~ normal(0, sigma_f);
  to_vector(fac_eti_raw) ~ normal(0, 1);

  for (it in 1:p) {
    tau[it] ~ normal(0, 1.5);
  }

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

  for (j in 1:nobs) {
    int ti = time[j];
    int it = itm[j];
    int pe = person[j];

    real eta_j = fac_gr[1, pe]
               + fac_gr[2, pe] * age_c[j]
               + fac_eti_raw[ti, pe] * eti_sd;

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
  // same name/shape as lmnlfa-quad.stan so existing R helpers work unchanged
  matrix[2, 2] Omega = multiply_lower_tri_self_transpose(L_Omega);
}
