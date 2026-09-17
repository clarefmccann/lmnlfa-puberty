// lmnlfa-growth-sigmoid-impact.stan
// Stage A of the staged SIGMOIDAL longitudinal MNLFA: growth-factor impact
// ONLY (race, WHtR on the rate and inflection-age means), no DIF at all --
// mirrors lmnlfa-growth-impact.stan's role exactly, just with a bounded
// logistic growth curve instead of a linear one:
//   eta_j = 1 + 4 / (1 + exp(-k_i * (age_c[j] - alpha_i)))
// Fixed floor = 1, fixed ceiling = 5 (not estimated -- these are the
// requested bounds). k_i (rate) and alpha_i (inflection age) are
// person-specific, correlated random effects.
//
// IDENTIFICATION: unlike every other model in this project, there is NO
// marker-item constraint here (lp is a fully free, positive vector over
// all p items). This is deliberate, not an oversight: eta_j's own scale
// and location are already pinned by the FIXED constants 1 and 4 in the
// formula above (sig_j = inv_logit(...) is bounded [0,1] by construction,
// with no free parameter able to rescale it), which resolves the usual
// IRT scale+location indeterminacy on its own. Fixing lp[1] = 1 on top of
// that would be redundant, not additionally identifying.
//
// RECENTERING: eta_j as computed above lives in ~[1, 5], not the ~0-
// centered scale every other item-parameter prior in this project
// (np ~ normal(0, sigma_nu), tau ~ normal(0, 1.5)) was calibrated for.
// Rather than re-deriving those priors, (eta_j - 3) is what actually
// enters the item response equations -- same effective spread, existing
// priors stay valid. The reported/saved eta_j (see fac_gr, not shown
// directly -- k_i/alpha_i are what's saved; eta itself is a function of
// age, not a stored parameter) is still the true, interpretable
// [1,5]-ish quantity.
//
// RATE PRIOR: k_i = exp(fac_gr[1, i]) needs its own, TIGHTER prior scale
// (sigma_k, separate from the general sigma_f used for phi_int/phi_slp
// elsewhere in this project) -- reusing a wide generic prior here risks
// numerically extreme rates (near-step-function or near-flat sigmoids)
// that are likely to cause real identifiability problems, distinct from
// (and probably worse than) the convergence issues the linear model's
// own tuning already had to work through.
//
// Age-as-DIF, WHtR dual-role (baseline for impact, occasion-specific for
// DIF), and the wobble term are otherwise unchanged from
// lmnlfa-growth-impact.stan -- see that file's header for the full
// rationale.

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
  matrix[ni, kimp] ximp;  // growth-curve impact covariates: race_c1-3, whtr_c

  real<lower=0> sigma_l;
  real<lower=0> sigma_nu;
  real<lower=0> sigma_cor;
  real<lower=0> sigma_f;   // prior scale for alpha (inflection age) terms
  real<lower=0> sigma_k;   // prior scale for log-rate terms (tighter than sigma_f)
}

parameters {
  vector<lower=0>[p] lp;   // ALL item loadings free (no marker item -- see header)
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

  // fac_gr[1, k] = log-rate, fac_gr[2, k] = inflection age (age_c scale)
  matrix[2, ni] fac_gr;
  for (k in 1:ni) {
    vector[2] mu_eta;
    mu_eta[1] = mu_logk + ximp[k] * b_mu_logk;
    mu_eta[2] = mu_alpha + ximp[k] * b_mu_alpha;
    fac_gr[, k] = mu_eta + diag_pre_multiply(phi_gr, L_Omega) * fac_dist[, k];
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
    // true [1,5]-ish latent puberty level, plus occasion-specific wobble
    real eta_true = 1 + 4 * sig + fac_eti_raw[ti, pe] * eti_sd;
    real eta_j = eta_true - 3; // recentered for the measurement model

    if (is_binary[it] == 1) {
      y[j] ~ bernoulli_logit(np[it] + lp[it] * eta_j);
    } else {
      y[j] ~ ordered_logistic(np[it] + lp[it] * eta_j, head(tau[it], k_item[it] - 1));
    }
  }
}
