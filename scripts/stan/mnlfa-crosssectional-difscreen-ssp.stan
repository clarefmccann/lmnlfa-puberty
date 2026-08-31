// mnlfa-crosssectional-difscreen-ssp.stan
// Replaces the plain-normal-prior + CI-threshold DIF screening in
// mnlfa-crosssectional-difscreen.stan with genuine spike-and-slab
// regularization (SSP), matching Chen & Bauer (2024, Psych Methods)
// Equations 14-17 and Appendix A1.
//
// For each item x DIF-covariate cell (i, k):
//   - a Beta(0.5, 0.5) inclusion parameter r[i,k] (U-shaped: pushes toward
//     0 or 1, a continuous relaxation of "is this cell in the model")
//   - a Bayesian-lasso (Laplace) "slab" magnitude on each of loading-DIF
//     and intercept-DIF, l_dif_star[i,k] and n_dif_star[i,k]
//   - actual DIF = slab magnitude x inclusion: l_dif = l_dif_star .* r,
//     n_dif = n_dif_star .* r -- one shared r per cell governs both,
//     matching the paper's "DIF effects on the intercept and the loading
//     from each covariate were assigned one shared inclusion parameter."
// The lasso penalties (phi_l for loadings, phi_n for intercepts) get a
// Gamma(4000, 200) hyperprior, matching the paper's empirical example --
// at its mean this implies a fairly tight default DIF prior (~0.05 SD),
// so only real effects escape the shrinkage.
//
// Impact (b_mu, b_phi) stays fixed at Stage A's estimates (data, not
// parameters), same as the plain screening model -- unrelated to the
// SSP-vs-plain-prior distinction, this piece was already correct.
//
// Selection: after fitting, threshold each cell's posterior mean
// inclusion parameter at 0.8 (the paper's reported threshold from prior
// simulation work) to decide which DIF terms carry into the final model.

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
}

parameters {
  vector<lower=0>[p - 1] lp_free;
  vector[p] np;
  array[p] ordered[k_max - 1] tau;

  real<lower=0> phi0;

  matrix<lower=0, upper=1>[p, kdif] r;   // shared inclusion parameter per cell
  matrix[p, kdif] l_dif_star;            // latent loading-DIF magnitude (lasso slab)
  matrix[p, kdif] n_dif_star;            // latent intercept-DIF magnitude (lasso slab)
  real<lower=0> phi_l;                   // lasso penalty, loadings
  real<lower=0> phi_n;                   // lasso penalty, intercepts

  vector[ni] eta_raw;
}

transformed parameters {
  vector[p] lp;
  lp[1] = 1;
  lp[2:p] = lp_free;

  matrix[p, kdif] l_dif = l_dif_star .* r;
  matrix[p, kdif] n_dif = n_dif_star .* r;
}

model {
  real u = 1;  // fixed SD scale, matching the paper's u = 1

  lp_free ~ normal(1, sigma_l);
  np      ~ normal(0, sigma_nu);
  for (it in 1:p) {
    tau[it] ~ normal(0, 1.5);
  }

  phi0 ~ normal(0, sigma_f);

  phi_l ~ gamma(4000, 200);
  phi_n ~ gamma(4000, 200);

  to_vector(r) ~ beta(0.5, 0.5);
  to_vector(l_dif_star) ~ double_exponential(0, u / phi_l);
  to_vector(n_dif_star) ~ double_exponential(0, u / phi_n);

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
