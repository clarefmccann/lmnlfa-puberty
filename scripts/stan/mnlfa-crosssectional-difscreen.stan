// mnlfa-crosssectional-difscreen.stan
// Stage B of the staged cross-sectional MNLFA procedure: DIF screening with
// impact FIXED at Stage A's posterior means (b_mu_fixed, b_phi_fixed are
// data, not parameters). This breaks the impact-vs-DIF confound directly:
// DIF terms can no longer "steal" variance that impact would otherwise
// explain, because impact isn't free to move anymore. phi0 stays free (it's
// a baseline scale nuisance parameter, not in direct competition with DIF
// the way b_mu/b_phi are).
//
// Every item x DIF-covariate combination is freely estimated here (no
// pre-screening within this stage) -- the output of this stage is used to
// decide which item x covariate DIF terms are worth carrying into the
// final combined fit (Stage C), analogous to the source paper's item-by-
// item DIF screening step, but done jointly across items here for
// tractability rather than one item at a time.

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
  int<lower=1> kdif;
  matrix[ni, kimp]   ximp;
  matrix[nobs, kdif] xdif;

  vector[kimp] b_mu_fixed;
  vector[kimp] b_phi_fixed;

  real<lower=0> sigma_l;
  real<lower=0> sigma_nu;
  real<lower=0> sigma_f;
  real<lower=0> sigma_di;
}

parameters {
  vector<lower=0>[p - 1] lp_free;
  vector[p] np;
  array[p] ordered[k_max - 1] tau;

  real<lower=0> phi0;

  matrix[p, kdif] l_dif;
  matrix[p, kdif] n_dif;

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

  phi0 ~ normal(0, sigma_f);

  to_vector(l_dif) ~ normal(0, sigma_di);
  to_vector(n_dif) ~ normal(0, sigma_di);

  eta_raw ~ normal(0, 1);

  vector[ni] mu_eta = ximp * b_mu_fixed;
  vector[ni] sd_eta = phi0 * exp(ximp * b_phi_fixed);
  vector[ni] eta = mu_eta + sd_eta .* eta_raw;

  for (j in 1:nobs) {
    int it = itm[j];
    int pe = person[j];
    real nu  = np[it] + xdif[j] * n_dif[it]';
    real lam = lp[it] * exp(xdif[j] * l_dif[it]');
    y[j] ~ ordered_logistic(nu + lam * eta[pe], head(tau[it], k_item[it] - 1));
  }
}
