// lmnlfa-growth-sigmoid-final.stan
// Stage C of the staged SIGMOIDAL longitudinal MNLFA: growth-curve impact
// free again (b_mu_logk, b_mu_alpha estimated), combined with ONLY the
// item x covariate DIF terms retained from Stage B's screening. Mirrors
// lmnlfa-growth-final.stan's sparse-DIF mechanism exactly -- see
// lmnlfa-growth-sigmoid-impact.stan's header for the identification,
// recentering, and rate-prior rationale (all unchanged here).

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

  matrix[p, kdif] l_pattern;
  matrix[p, kdif] n_pattern;
  int<lower=0> ml;
  int<lower=0> mn;

  real<lower=0> sigma_l;
  real<lower=0> sigma_nu;
  real<lower=0> sigma_cor;
  real<lower=0> sigma_f;
  real<lower=0> sigma_k;
  real<lower=0> sigma_di;
}

parameters {
  vector<lower=0>[p] lp;
  vector[p] np;

  real mu_logk;
  vector[kimp] b_mu_logk;
  real<lower=0> phi_logk;

  real mu_alpha;
  vector[kimp] b_mu_alpha;
  real<lower=0> phi_alpha;

  real z_cor;

  matrix[2, ni] fac_dist;

  matrix[d, ni] fac_eti_raw;
  real<lower=0> eti_sd;

  array[p] ordered[k_max - 1] tau;

  vector[ml] l_dif_free;
  vector[mn] n_dif_free;
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
    mu_eta[1] = mu_logk + ximp[k] * b_mu_logk;
    mu_eta[2] = mu_alpha + ximp[k] * b_mu_alpha;
    fac_gr[, k] = mu_eta + diag_pre_multiply(phi_gr, L_Omega) * fac_dist[, k];
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
  lp ~ normal(1, sigma_l);
  np ~ normal(0, sigma_nu);

  mu_logk    ~ normal(0, sigma_k);
  b_mu_logk  ~ normal(0, sigma_k);
  phi_logk   ~ normal(0, sigma_k);

  mu_alpha   ~ normal(0, sigma_f);
  b_mu_alpha ~ normal(0, sigma_f);
  phi_alpha  ~ normal(0, sigma_f);

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
