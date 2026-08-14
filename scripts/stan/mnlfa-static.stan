// mnlfa-static.stan
// Stage 1+2 debug model: single-wave, cross-sectional MNLFA.
// Mean + variance impact of one predictor (age_c) on the latent factor;
// no DIF, no growth structure. Used to confirm the measurement model and
// impact structure converge before adding longitudinal complexity in
// lmnlfa-quad.stan.

data {
  int<lower=1> nobs;
  int<lower=2> p;  // >=2 so item 1 can serve as the scale-identifying marker
  int<lower=1> ni;

  array[nobs] int<lower=1, upper=ni> person;
  array[nobs] int<lower=1, upper=p>  itm;

  array[nobs] int y;
  array[p] int<lower=2> k_item;
  int<lower=1> k_max;

  vector[ni] age_c;  // person-level centered age (impact predictor)

  real<lower=0> sigma_l;
  real<lower=0> sigma_nu;
  real<lower=0> sigma_f;
}

parameters {
  // item 1's loading is fixed to 1 (see transformed parameters) as a marker
  // item for scale identification -- lp and phi0 are otherwise only jointly
  // identified through their product, which produces a strong non-identified
  // ridge (see mnlfa_static.R fit diagnostics: cor(lp, phi0) ~ -0.92).
  vector<lower=0>[p - 1] lp_free;
  vector[p] np;
  array[p] ordered[k_max - 1] tau;

  real b_mu;            // mean impact of age_c on the factor
  real b_phi;           // variance impact of age_c on the factor (log scale)
  real<lower=0> phi0;   // factor SD at age_c = 0

  vector[ni] eta_raw;   // non-centered latent factor
}

transformed parameters {
  vector[p] lp;
  lp[1] = 1;
  lp[2:p] = lp_free;

  vector[ni] mu_eta = b_mu * age_c;
  vector<lower=0>[ni] sd_eta = phi0 * exp(b_phi * age_c);
  vector[ni] eta = mu_eta + sd_eta .* eta_raw;
}

model {
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
  np ~ normal(0, sigma_nu);
  for (it in 1:p) {
    tau[it] ~ normal(0, 1.5);
  }

  b_mu  ~ normal(0, sigma_f);
  b_phi ~ normal(0, sigma_f);
  phi0  ~ normal(0, sigma_f);

  eta_raw ~ normal(0, 1);

  for (j in 1:nobs) {
    int it = itm[j];
    int pe = person[j];
    y[j] ~ ordered_logistic(np[it] + lp[it] * eta[pe], head(tau[it], k_item[it] - 1));
  }
}
