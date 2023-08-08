data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // population-level effects
  real Intercept;  // temporary intercept for centered predictors
  real<lower=0> sigma;  // dispersion parameter
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // standardized group-level effects
  
  // Manually added
   vector[N_1] r_1_1;  // actual group-level effects
}
//transformed parameters {
//  vector[N_1] r_1_1;  // actual group-level effects
//  real lprior = 0;  // prior contributions to the log posterior
//  r_1_1 = (sd_1[1] * (z_1[1]));
//  lprior += student_t_lpdf(Intercept | 3, 2, 2.5);
//  lprior += student_t_lpdf(sigma | 3, 0, 2.5)
//    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
//  lprior += student_t_lpdf(sd_1 | 3, 0, 2.5)
//    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
//}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    mu += Intercept;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n];
      
    }
    target += normal_id_glm_lpdf(Y | Xc, mu, b, sigma);
  }
  // priors including constants
  // target += lprior;
  // target += std_normal_lpdf(z_1[1]);

  target += student_t_lpdf(Intercept | 3, 2, 2.5);
  target += student_t_lpdf(sigma | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += student_t_lpdf(sd_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
    
  // Manually added
  target += normal_lpdf(r_1_1 | 0, sd_1[1]);

  // manual soft sum to zero constraint
  // target += normal_lpdf(sum(r_1_1) | 0, 0.0001 *N_1);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
}