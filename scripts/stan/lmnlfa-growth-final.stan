// lmnlfa-growth-final.stan
// Stage C of the staged longitudinal MNLFA procedure: growth-factor impact
// free again (b_mu_int, b_mu_slp estimated), combined with ONLY the item x
// covariate DIF terms retained from Stage B's screening -- not the full
// saturated set. l_pattern/n_pattern (data, 0/1) mark which cells get a
// free parameter; everything else is fixed at 0. Mirrors
// mnlfa-crosssectional-final.stan's sparse-DIF mechanism exactly.
//
// If Stage B finds no significant DIF anywhere, ml = mn = 0 and this
// degenerates back to the impact-only model (Stage A) -- Stan handles
// zero-length parameter vectors fine.

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
  matrix[nobs, kdif] xdif;

  matrix[p, kdif] l_pattern;  // 0/1: which item x covariate loading-DIF terms are free
  matrix[p, kdif] n_pattern;  // 0/1: same for intercept-DIF
  int<lower=0> ml;            // sum(l_pattern)
  int<lower=0> mn;            // sum(n_pattern)

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
  vector[kimp] b_mu_int;
  vector[kimp] b_mu_slp;
  real<lower=0> phi_int;
  real<lower=0> phi_slp;

  real z_cor;

  matrix[2, ni] fac_dist;

  matrix[d, ni] fac_eti_raw;
  real<lower=0> eti_sd;

  array[p] ordered[k_max - 1] tau;

  vector[ml] l_dif_free;
  vector[mn] n_dif_free;
}

transformed parameters {
  vector[p] lp;
  lp[1] = 1;
  lp[2:p] = lp_free;

  vector<lower=0>[2] phi_eta;
  phi_eta[1] = phi_int;
  phi_eta[2] = phi_slp;

  real<lower=-1, upper=1> rho = tanh(z_cor);
  matrix[2, 2] L_Omega;
  L_Omega[1, 1] = 1;
  L_Omega[2, 1] = rho;
  L_Omega[1, 2] = 0;
  L_Omega[2, 2] = sqrt(1 - square(rho));

  matrix[2, ni] fac_gr;
  for (k in 1:ni) {
    vector[2] mu_eta;
    mu_eta[1] = ximp[k] * b_mu_int;
    mu_eta[2] = mu_slp + ximp[k] * b_mu_slp;
    fac_gr[, k] = mu_eta + diag_pre_multiply(phi_eta, L_Omega) * fac_dist[, k];
  }

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

  mu_slp   ~ normal(0, sigma_f);
  b_mu_int ~ normal(0, sigma_f);
  b_mu_slp ~ normal(0, sigma_f);
  phi_int  ~ normal(0, sigma_f);
  phi_slp  ~ normal(0, sigma_f);

  z_cor ~ normal(0, sigma_cor);
  to_vector(fac_dist) ~ normal(0, 1);

  eti_sd ~ normal(0, sigma_f);
  to_vector(fac_eti_raw) ~ normal(0, 1);

  l_dif_free ~ normal(0, sigma_di);
  n_dif_free ~ normal(0, sigma_di);

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

    real nu  = np[it] + xdif[j] * n_dif[it]';
    real lam = lp[it] * exp(xdif[j] * l_dif[it]');

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
