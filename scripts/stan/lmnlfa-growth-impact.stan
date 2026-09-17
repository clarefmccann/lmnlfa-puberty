// lmnlfa-growth-impact.stan
// Stage A of the staged longitudinal MNLFA procedure: growth-factor impact
// ONLY (race, WHtR on the intercept and slope means), no DIF at all -- not
// even informant, matching mnlfa-crosssectional-impact.stan's precedent.
// This isolates clean impact estimates before DIF is introduced, avoiding
// the impact-vs-DIF confound diagnosed in the cross-sectional model's
// fully-saturated fit (Rhat 1.37-1.54, ESS 7-9/4000 -- age's effect on
// item responses was only weakly separable between "real impact" and
// "DIF" when both were estimated freely together with no screening).
//
// Age is NOT an impact covariate here (unlike the cross-sectional model):
// age is already this model's growth axis (mu_slp *is* the age effect),
// so "age impact on the growth factor" isn't a separable concept the way
// it is cross-sectionally. Age-varying DIF -- item parameters shifting
// with age independent of the latent trajectory -- is the correct
// longitudinal analog, and is tested in Stage B/C instead.
//
// Growth structure (intercept/slope/correlation via tanh, marker-item
// identification, occasion-specific "wobble" residual) is otherwise
// unchanged from the validated growth-only/informant-only models
// (lmnlfa-quad-tanhcor.stan / lmnlfa-linear-tanhcor-informant.stan).

data {
  int<lower=1> nobs;
  int<lower=2> p;  // >=2 so item 1 can serve as the scale-identifying marker
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
  matrix[ni, kimp] ximp;  // growth-factor impact covariates: race_c1-3, whtr_c

  real<lower=0> sigma_l;
  real<lower=0> sigma_nu;
  real<lower=0> sigma_cor;
  real<lower=0> sigma_f;
}

parameters {
  vector<lower=0>[p - 1] lp_free;
  vector[p] np;

  real mu_slp;
  vector[kimp] b_mu_int;  // race/WHtR impact on the intercept mean
  vector[kimp] b_mu_slp;  // race/WHtR impact on the slope mean
  real<lower=0> phi_int;
  real<lower=0> phi_slp;

  real z_cor;

  matrix[2, ni] fac_dist;

  matrix[d, ni] fac_eti_raw;
  real<lower=0> eti_sd;

  array[p] ordered[k_max - 1] tau;
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

  // per-person growth factor (intercept, slope), now with race/WHtR impact
  // on both means
  matrix[2, ni] fac_gr;
  for (k in 1:ni) {
    vector[2] mu_eta;
    mu_eta[1] = ximp[k] * b_mu_int;
    mu_eta[2] = mu_slp + ximp[k] * b_mu_slp;
    fac_gr[, k] = mu_eta + diag_pre_multiply(phi_eta, L_Omega) * fac_dist[, k];
  }
}

model {
  lp_free ~ normal(1, sigma_l);
  np      ~ normal(0, sigma_nu);

  mu_slp   ~ normal(0, sigma_f);
  b_mu_int ~ normal(0, sigma_f);
  b_mu_slp ~ normal(0, sigma_f);
  phi_int  ~ normal(0, sigma_f);
  phi_slp  ~ normal(0, sigma_f);

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

    if (is_binary[it] == 1) {
      y[j] ~ bernoulli_logit(np[it] + lp[it] * eta_j);
    } else {
      y[j] ~ ordered_logistic(np[it] + lp[it] * eta_j, head(tau[it], k_item[it] - 1));
    }
  }
}

generated quantities {
  matrix[2, 2] Omega = multiply_lower_tri_self_transpose(L_Omega);
}
