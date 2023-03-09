data {
  int<lower=1> K;  // number of clusters, known at the time of post-processing
  int<lower=1> p;  // number of food items
  int<lower=1> d;  // number of consumption levels
  int<lower=1> n;  // number of subjects
  int<lower=1> q;  // number of covariates in probit regression
  
  array[n, p] int X;                // categorical food data
  array[n] int<lower=0, upper=1> y; // binary outcome data
  array[n, q] real V;               // covariate matrix excluding class assignment
  vector<lower=0>[n] weights;       // individual-level survey weights
  
  vector[K] alpha;         // hyperparameter for pi prior
  array[K] vector[d] eta;  // hyperparameter for theta prior
  vector[q] mu0;           // hyperparameter for mean of xi prior
  cov_matrix[q] Sig0;      // hyperparameter for covariance of xi prior
}
parameters {
  simplex[K] pi;                 // cluster probabilities
  array[p, K] simplex[d] theta;  // cluster-specific item consumption probabilities
  array[K] vector[q] xi;         // regression coefficients
}
transformed parameters {
  array[n] vector[K] log_cond_c; // log p(c_i=k| -)
  for (i in 1:n) {
    for (k in 1:K) {
      log_cond_c[i, k] = log(pi[k]) + bernoulli_lpmf(y[i] | Phi(to_row_vector(V[i, ]) * xi[k]));
      for (j in 1:p) {
        log_cond_c[i, k] = log_cond_c[i, k] + log(theta[j, k, X[i, j]]);
      }
    }
  }
}
model {
  pi ~ dirichlet(alpha);  // prior for pi
  for (j in 1:p) {        // prior for theta
    for (k in 1:K) {
      theta[j, k] ~ dirichlet(eta[k]);
    }
  }
  for (k in 1:K) {
    xi[k] ~ multi_normal(mu0, Sig0);
  }

  for (i in 1:n) {  // increment log-probability
    target += weights[i] * log_sum_exp(log_cond_c[i]);
  }
}

