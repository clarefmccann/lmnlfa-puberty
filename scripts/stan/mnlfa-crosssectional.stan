// mnlfa-crosssectional.stan
// Cross-sectional MNLFA: one randomly selected wave per person, informant
// (parent/youth) as a within-person covariate rather than doubling the item
// set. Generalizes mnlfa-static.stan (which already validated cleanly: 0
// divergences, 0 max-treedepth) from a single impact covariate to K impact
// covariates (mean + variance) and a separate set of M DIF covariates on
// every item's loading and intercept.
//
// Design notes:
//   - Impact covariates (ximp, person-level: age, race contrasts, WHtR)
//     affect the population distribution of the latent factor.
//   - DIF covariates (xdif, observation-level: impact covariates + informant)
//     affect how items are scored. Informant is DIF-only -- a person's true
//     latent trait doesn't have an "informant effect", only the observed
//     item responses can.
//   - No item-by-covariate screening (ldf) here: every item is tested for
//     DIF on every covariate in xdif, matching the source paper's saturated
//     empirical-model approach rather than the staged screening procedure.
//   - Same identification fixes validated on the static/growth-only models:
//     item 1's loading fixed to 1 (marker item, avoids the loading/factor-SD
//     ridge), lambda ~ N+(1, sigma_l), nu centered at 0 (not the paper's -2,
//     since these ordinal items are not rare/skewed the way the paper's
//     binary items are).

data {
  int<lower=1> nobs;
  int<lower=2> p;   // number of items (4)
  int<lower=1> ni;  // number of unique people

  array[nobs] int<lower=1, upper=ni> person;
  array[nobs] int<lower=1, upper=p>  itm;

  array[nobs] int y;
  array[p] int<lower=2> k_item;
  int<lower=1> k_max;

  int<lower=1> kimp;              // number of impact covariates
  int<lower=1> kdif;              // number of DIF covariates
  matrix[ni, kimp]   ximp;        // person-level impact covariates
  matrix[nobs, kdif] xdif;        // observation-level DIF covariates

  real<lower=0> sigma_l;
  real<lower=0> sigma_nu;
  real<lower=0> sigma_f;
  real<lower=0> sigma_di;
}

parameters {
  vector<lower=0>[p - 1] lp_free;  // item 1 fixed to 1 (marker item)
  vector[p] np;
  array[p] ordered[k_max - 1] tau;

  vector[kimp] b_mu;    // mean impact
  vector[kimp] b_phi;   // variance impact (log scale)
  real<lower=0> phi0;   // factor SD at ximp = 0

  matrix[p, kdif] l_dif;  // loading DIF, item x covariate
  matrix[p, kdif] n_dif;  // intercept DIF, item x covariate

  vector[ni] eta_raw;   // non-centered latent factor
}

transformed parameters {
  vector[p] lp;
  lp[1] = 1;
  lp[2:p] = lp_free;
}

model {
  lp_free ~ normal(1, sigma_l);
  np      ~ normal(0, sigma_nu);
  for (it in 1:p) {
    tau[it] ~ normal(0, 1.5);
  }

  b_mu  ~ normal(0, sigma_f);
  b_phi ~ normal(0, sigma_f);
  phi0  ~ normal(0, sigma_f);

  to_vector(l_dif) ~ normal(0, sigma_di);
  to_vector(n_dif) ~ normal(0, sigma_di);

  eta_raw ~ normal(0, 1);

  vector[ni] mu_eta = ximp * b_mu;
  vector[ni] sd_eta = phi0 * exp(ximp * b_phi);  // exp() already ensures > 0
  vector[ni] eta = mu_eta + sd_eta .* eta_raw;

  for (j in 1:nobs) {
    int it = itm[j];
    int pe = person[j];
    real nu  = np[it] + xdif[j] * n_dif[it]';
    real lam = lp[it] * exp(xdif[j] * l_dif[it]');
    y[j] ~ ordered_logistic(nu + lam * eta[pe], head(tau[it], k_item[it] - 1));
  }
}
