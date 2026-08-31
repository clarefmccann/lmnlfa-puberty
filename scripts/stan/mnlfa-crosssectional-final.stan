// mnlfa-crosssectional-final.stan
// Stage C of the staged cross-sectional MNLFA procedure: impact free again
// (b_mu, b_phi estimated, not fixed), combined with ONLY the item x
// covariate DIF terms retained from Stage B's screening -- not the full
// saturated set. l_pattern/n_pattern (data, 0/1) mark which cells get a
// free parameter; everything else is fixed at 0. This mirrors the sparse-
// DIF mechanism already used in lmnlfa-quad.stan (there: item x
// {time-varying, invariant}; here: item x covariate), generalized to a
// full 2D grid since DIF covariates aren't split into just two categories.
//
// If Stage B finds no significant DIF anywhere, ml = mn = 0 and this
// degenerates back to the impact-only model (Stage A) -- Stan handles
// zero-length parameter vectors fine, same as lmnlfa-quad.stan already
// does when mf = mtv = 0.

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

  matrix[p, kdif] l_pattern;  // 0/1: which item x covariate loading-DIF terms are free
  matrix[p, kdif] n_pattern;  // 0/1: same for intercept-DIF
  int<lower=0> ml;            // sum(l_pattern)
  int<lower=0> mn;            // sum(n_pattern)

  real<lower=0> sigma_l;
  real<lower=0> sigma_nu;
  real<lower=0> sigma_f;
  real<lower=0> sigma_di;
}

parameters {
  vector<lower=0>[p - 1] lp_free;
  vector[p] np;
  array[p] ordered[k_max - 1] tau;

  vector[kimp] b_mu;
  vector[kimp] b_phi;
  real<lower=0> phi0;

  vector[ml] l_dif_free;
  vector[mn] n_dif_free;

  vector[ni] eta_raw;
}

transformed parameters {
  vector[p] lp;
  lp[1] = 1;
  lp[2:p] = lp_free;

  matrix[p, kdif] l_dif = rep_matrix(0, p, kdif);
  matrix[p, kdif] n_dif = rep_matrix(0, p, kdif);
  {
    int tmp = 0;
    for (i in 1:p) {
      for (k in 1:kdif) {
        if (l_pattern[i, k] == 1) { tmp += 1; l_dif[i, k] = l_dif_free[tmp]; }
      }
    }
    tmp = 0;
    for (i in 1:p) {
      for (k in 1:kdif) {
        if (n_pattern[i, k] == 1) { tmp += 1; n_dif[i, k] = n_dif_free[tmp]; }
      }
    }
  }
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

  l_dif_free ~ normal(0, sigma_di);
  n_dif_free ~ normal(0, sigma_di);

  eta_raw ~ normal(0, 1);

  vector[ni] mu_eta = ximp * b_mu;
  vector[ni] sd_eta = phi0 * exp(ximp * b_phi);
  vector[ni] eta = mu_eta + sd_eta .* eta_raw;

  for (j in 1:nobs) {
    int it = itm[j];
    int pe = person[j];
    real nu  = np[it] + xdif[j] * n_dif[it]';
    real lam = lp[it] * exp(xdif[j] * l_dif[it]');
    y[j] ~ ordered_logistic(nu + lam * eta[pe], head(tau[it], k_item[it] - 1));
  }
}
