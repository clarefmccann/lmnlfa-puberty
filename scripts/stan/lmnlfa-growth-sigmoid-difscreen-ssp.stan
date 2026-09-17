// lmnlfa-growth-sigmoid-difscreen-ssp.stan
// Alternative to lmnlfa-growth-sigmoid-difscreen.stan's Stage B: same growth
// curve, same impact-fixed setup, but DIF terms are selected via Bayesian
// spike-and-slab regularization (Chen & Bauer, 2024, "Modeling Construct
// Change Over Time Amidst Potential Changes in Construct Measurement:
// A Longitudinal Moderated Factor Analysis Approach") instead of freely
// estimating every DIF term with a generic diffuse prior and screening
// them post hoc with Benjamini-Hochberg FDR.
//
// SPIKE-AND-SLAB CONSTRUCTION (mirrors Chen & Bauer Eqs. 14-15):
//   l_dif[i,k] = l_star[i,k] * r_incl[i,k]
//   n_dif[i,k] = n_star[i,k] * r_incl[i,k]
//   l_star[i,k], n_star[i,k] ~ double_exponential(0, ssp_scale / phi[k])
//   r_incl[i,k] ~ beta(0.5, 0.5)
//   phi_l[k], phi_n[k] ~ gamma(phi_shape, phi_rate)
//
// r_incl is the SHARED inclusion parameter for item i x covariate k,
// governing BOTH the loading and intercept DIF terms for that cell (Chen &
// Bauer: "DIF effects on the intercept and the loading from each covariate
// were assigned one shared inclusion parameter to improve DIF evaluation
// accuracy on the item level"). A separate Laplace/lasso penalty (phi) is
// estimated per covariate and per parameter type (loading vs intercept),
// each with its own gamma hyperprior, matching the paper's per-covariate
// penalty structure.
//
// After fitting, DIF is "selected" by thresholding the posterior mean of
// r_incl (Chen & Bauer used 0.7-0.8) rather than BH-FDR + magnitude floor.
// No anchor items are declared in advance -- unimportant DIF effects are
// pulled toward zero by the lasso prior and their inclusion parameter stays
// near its 0.5 prior mean; real effects push r_incl toward 1.
//
// NOTE: ssp_scale, phi_shape, phi_rate are NEW data inputs (not present in
// the BH-FDR difscreen.stan) -- see the R driver for how they're set. They
// are NOT copied from Chen & Bauer's own hyperprior values, which were
// calibrated for their binary items and different covariate scaling; ours
// are a reasonable starting point for this model's effect-size scale and
// may need tuning -- watch r_incl's posterior distribution (it should be
// clearly bimodal-ish, not uniformly stuck at 0.5) as a diagnostic.

data {
  int<lower=1> nobs;
  int<lower=2> p;
  int<lower=1> ni;
  int<lower=1> d;

  array[nobs] int<lower=1, upper=ni> person;
  array[nobs] int<lower=1, upper=p>  itm;
  array[nobs] int<lower=1, upper=d>  time;

  vector[nobs] age_c;

  array[nobs] int y;
  array[p]    int<lower=0, upper=1> is_binary;

  array[p] int<lower=2> k_item;
  int<lower=1> k_max;

  int<lower=1> kimp;
  int<lower=1> kdif;
  matrix[ni, kimp]   ximp;
  matrix[nobs, kdif] xdif;  // age_c, race_c1-3, whtr_c, informant_c

  vector[kimp] b_mu_logk_fixed;
  vector[kimp] b_mu_alpha_fixed;

  real<lower=0> sigma_l;
  real<lower=0> sigma_nu;
  real<lower=0> sigma_cor;
  real<lower=0> sigma_f;
  real<lower=0> sigma_k;

  // spike-and-slab hyperparameters (NOT used by the BH-FDR difscreen model)
  real<lower=0> ssp_scale;    // "u" in Chen & Bauer Eq. 14 -- known SD-like scale
  real<lower=0> phi_shape;    // gamma hyperprior shape for lasso penalties
  real<lower=0> phi_rate;     // gamma hyperprior rate for lasso penalties
}

parameters {
  vector<lower=0>[p] lp;
  vector[p] np;

  real mu_logk;
  real<lower=0> phi_logk;

  real mu_alpha;
  real<lower=0> phi_alpha;

  real z_cor;

  matrix[2, ni] fac_dist;

  matrix[d, ni] fac_eti_raw;
  real<lower=0> eti_sd;

  array[p] ordered[k_max - 1] tau;

  // spike-and-slab DIF parameters
  matrix[p, kdif] l_star;
  matrix[p, kdif] n_star;
  matrix<lower=0, upper=1>[p, kdif] r_incl;
  vector<lower=0>[kdif] phi_l;
  vector<lower=0>[kdif] phi_n;
}

transformed parameters {
  vector<lower=0>[2] phi_gr;
  phi_gr[1] = phi_logk;
  phi_gr[2] = phi_alpha;

  real<lower=-1, upper=1> rho = tanh(z_cor);
  matrix[2, 2] L_Omega;
  L_Omega[1, 1] = 1;
  L_Omega[2, 1] = rho;
  L_Omega[1, 2] = 0;
  L_Omega[2, 2] = sqrt(1 - square(rho));

  matrix[2, ni] fac_gr;
  for (k in 1:ni) {
    vector[2] mu_eta;
    mu_eta[1] = mu_logk + ximp[k] * b_mu_logk_fixed;
    mu_eta[2] = mu_alpha + ximp[k] * b_mu_alpha_fixed;
    fac_gr[, k] = mu_eta + diag_pre_multiply(phi_gr, L_Omega) * fac_dist[, k];
  }

  // spike-and-slab reconstruction: shared inclusion parameter per (item,
  // covariate) cell gates BOTH the loading and intercept DIF term
  matrix[p, kdif] l_dif = l_star .* r_incl;
  matrix[p, kdif] n_dif = n_star .* r_incl;
}

model {
  lp ~ normal(1, sigma_l);
  np ~ normal(0, sigma_nu);

  mu_logk  ~ normal(0, sigma_k);
  phi_logk ~ normal(0, sigma_k);

  mu_alpha  ~ normal(0, sigma_f);
  phi_alpha ~ normal(0, sigma_f);

  z_cor ~ normal(0, sigma_cor);
  to_vector(fac_dist) ~ normal(0, 1);

  eti_sd ~ normal(0, sigma_f);
  to_vector(fac_eti_raw) ~ normal(0, 1);

  for (it in 1:p) {
    tau[it] ~ normal(0, 1.5);
  }

  // spike-and-slab priors on DIF: Laplace ("lasso") magnitude x Beta(.5,.5)
  // inclusion, one lasso penalty per covariate (shared across items) for
  // loadings and intercepts separately
  phi_l ~ gamma(phi_shape, phi_rate);
  phi_n ~ gamma(phi_shape, phi_rate);
  for (k in 1:kdif) {
    to_vector(l_star[, k]) ~ double_exponential(0, ssp_scale / phi_l[k]);
    to_vector(n_star[, k]) ~ double_exponential(0, ssp_scale / phi_n[k]);
  }
  to_vector(r_incl) ~ beta(0.5, 0.5);

  for (j in 1:nobs) {
    int ti = time[j];
    int it = itm[j];
    int pe = person[j];

    real k_i = exp(fac_gr[1, pe]);
    real alpha_i = fac_gr[2, pe];
    real sig = inv_logit(k_i * (age_c[j] - alpha_i));
    real eta_true = 1 + 4 * sig + fac_eti_raw[ti, pe] * eti_sd;
    real eta_j = eta_true - 3;

    real nu  = np[it] + xdif[j] * n_dif[it]';
    real lam = lp[it] * exp(xdif[j] * l_dif[it]');

    if (is_binary[it] == 1) {
      y[j] ~ bernoulli_logit(nu + lam * eta_j);
    } else {
      y[j] ~ ordered_logistic(nu + lam * eta_j, head(tau[it], k_item[it] - 1));
    }
  }
}
