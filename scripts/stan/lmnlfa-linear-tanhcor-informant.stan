// lmnlfa-linear-tanhcor-informant.stan
// Named "linear" (not "quad") deliberately: despite the "quad" in
// lmnlfa-quad-tanhcor.stan's name, that model (and every variant since,
// including this one) has never actually implemented a quadratic growth
// term -- age2_c was declared but never referenced in any equation. This
// file drops age2_c entirely rather than carry forward unused, misleading
// vestigial structure.
//
// Next increment beyond the validated growth-only model
// (lmnlfa-quad-tanhcor.stan): restructures from 8 separate items
// (parent/youth treated as unrelated items) to 4 shared items with
// informant (parent vs youth) as an explicit DIF covariate -- matching
// the cross-sectional analysis's design (mnlfa-crosssectional*.stan),
// which found this DIF to be real and pervasive in both sexes.
//
// A person's true pubertal trajectory doesn't depend on who's reporting
// it, only the item responses can -- so informant enters only as a DIF
// covariate (loading + intercept shift), never as growth-factor impact,
// same principle as the cross-sectional model.
//
// Growth structure (intercept/slope/correlation via tanh, marker-item
// identification, occasion-specific "wobble" residual) is unchanged from
// the validated growth-only model. The old nfpreds/xf/ldiff/ldiftv
// machinery (dormant in growth-only, since nfpreds=0 there) is dropped
// here rather than carried forward unused -- growth-factor impact
// (race, WHtR) and age-varying DIF are deferred to their own later stage,
// matching the file-per-stage pattern used throughout this project.
//
// fac_gr (the per-person growth factor) is computed in transformed
// parameters rather than as a model{}-local, so it's saved automatically
// for post-hoc factor-score extraction without duplicating the
// computation in a separate generated quantities block.

data {
  int<lower=1> nobs;
  int<lower=2> p;  // >=2 so item 1 can serve as the scale-identifying marker
  int<lower=1> ni;
  int<lower=1> d;

  array[nobs] int<lower=1, upper=ni> person;
  array[nobs] int<lower=1, upper=p>  itm;
  array[nobs] int<lower=1, upper=d>  time;

  vector[nobs] age_c;
  vector[nobs] informant_c;  // +1 parent, -1 youth (effect-coded)

  array[nobs] int y;
  array[p]    int<lower=0, upper=1> is_binary;

  array[p] int<lower=2> k_item;
  int<lower=1> k_max;

  real<lower=0> sigma_l;
  real<lower=0> sigma_nu;
  real<lower=0> sigma_cor;
  real<lower=0> sigma_f;
  real<lower=0> sigma_di;
}

parameters {
  vector<lower=0>[p - 1] lp_free;
  vector[p] np;

  real mu_slp;
  real<lower=0> phi_int;
  real<lower=0> phi_slp;

  real z_cor;

  matrix[2, ni] fac_dist;

  matrix[d, ni] fac_eti_raw;
  real<lower=0> eti_sd;

  array[p] ordered[k_max - 1] tau;

  vector[p] l_dif_informant;
  vector[p] n_dif_informant;
}

transformed parameters {
  vector[p] lp;
  lp[1] = 1;
  lp[2:p] = lp_free;

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

  // per-person growth factor (intercept, slope) -- computed here rather
  // than as a model{}-local so it's saved for post-hoc factor scores
  matrix[2, ni] fac_gr;
  for (k in 1:ni) {
    vector[2] mu_eta;
    mu_eta[1] = 0;
    mu_eta[2] = mu_slp;
    fac_gr[, k] = mu_eta + diag_pre_multiply(phi_eta, L_Omega) * fac_dist[, k];
  }
}

model {
  lp_free ~ normal(1, sigma_l);
  np      ~ normal(0, sigma_nu);

  l_dif_informant ~ normal(0, sigma_di);
  n_dif_informant ~ normal(0, sigma_di);

  mu_slp  ~ normal(0, sigma_f);
  phi_int ~ normal(0, sigma_f);
  phi_slp ~ normal(0, sigma_f);

  z_cor ~ normal(0, sigma_cor);
  to_vector(fac_dist) ~ normal(0, 1);

  eti_sd ~ normal(0, sigma_f);
  to_vector(fac_eti_raw) ~ normal(0, 1);

  for (it in 1:p) {
    tau[it] ~ normal(0, 1.5);
  }

  for (j in 1:nobs) {
    int ti = time[j];
    int it = itm[j];
    int pe = person[j];

    real eta_j = fac_gr[1, pe]
               + fac_gr[2, pe] * age_c[j]
               + fac_eti_raw[ti, pe] * eti_sd;

    real nu  = np[it] + n_dif_informant[it] * informant_c[j];
    real lam = lp[it] * exp(l_dif_informant[it] * informant_c[j]);

    if (is_binary[it] == 1) {
      y[j] ~ bernoulli_logit(nu + lam * eta_j);
    } else {
      y[j] ~ ordered_logistic(nu + lam * eta_j, head(tau[it], k_item[it] - 1));
    }
  }
}

generated quantities {
  matrix[2, 2] Omega = multiply_lower_tri_self_transpose(L_Omega);
}
