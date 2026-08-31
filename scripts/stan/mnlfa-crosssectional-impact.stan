// mnlfa-crosssectional-impact.stan
// Stage A of the staged cross-sectional MNLFA procedure: impact ONLY
// (age, race, WHtR on the factor mean/variance), no DIF at all. This
// isolates clean impact estimates before DIF is introduced, avoiding the
// impact-vs-DIF confound diagnosed in the fully-saturated fit (Rhat
// 1.37-1.54, ESS 7-9/4000 despite 0 divergences/max-treedepth -- a
// between-chain multimodality problem, not a geometry problem, caused by
// age's effect being only weakly separable between "real impact" and
// "DIF" when both are freely estimated together with no screening).
//
// Structurally identical to mnlfa-static.stan, generalized from 1 impact
// covariate to kimp. Same identification fixes: item 1's loading fixed to
// 1 (marker item), lambda ~ N+(1, sigma_l), nu centered at 0.

data {
  int<lower=1> nobs;
  int<lower=2> p;
  int<lower=1> ni;

  array[nobs] int<lower=1, upper=ni> person;
  array[nobs] int<lower=1, upper=p>  itm;

  array[nobs] int y;
  array[p] int<lower=2> k_item;
  int<lower=1> k_max;

  int<lower=1> kimp;
  matrix[ni, kimp] ximp;

  real<lower=0> sigma_l;
  real<lower=0> sigma_nu;
  real<lower=0> sigma_f;
}

parameters {
  vector<lower=0>[p - 1] lp_free;
  vector[p] np;
  array[p] ordered[k_max - 1] tau;

  vector[kimp] b_mu;
  vector[kimp] b_phi;
  real<lower=0> phi0;

  vector[ni] eta_raw;
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

  eta_raw ~ normal(0, 1);

  vector[ni] mu_eta = ximp * b_mu;
  vector[ni] sd_eta = phi0 * exp(ximp * b_phi);
  vector[ni] eta = mu_eta + sd_eta .* eta_raw;

  for (j in 1:nobs) {
    int it = itm[j];
    int pe = person[j];
    y[j] ~ ordered_logistic(np[it] + lp[it] * eta[pe], head(tau[it], k_item[it] - 1));
  }
}
