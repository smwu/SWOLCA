data {
  int K;  // number of clusters, known at the time of post-processing
  int p;  // number of food items
  int d;  // number of consumption levels
  int n;  // number of subjects
  int q;  // number of covariates in probit regression
  array[n, p] int X;  // categorical food data
  array[n] int y;  // binary outcome data
  //real V[n, q]; // covariate design matrix for probit regression. Categorical variables listed before continuous
  real V_k[n, q, K]; // covariate matrix where all units are assigned to k
}
parameters {
  simplex[K] pi;  // cluster probabilities
  real<lower=0, upper=1> theta[p, K, d];  // cluster-specific item consumption probabilities
  real xi[q];  // regression coefficients
  
}
model {
  real log_lik[K];
  for (i in 1:n) {
    for (k in 1:K) {
      real log_mult_k;  // multinomial log-pmf for cluster k
      for (j in 1:p) {
        vector[d] theta_k = to_vector(theta[j,k,]);  // theta[j,k,] vector of consumption probabilities
        log_mult_k += categorical_lpmf(X[i,j] | theta_k);
      }
      log_lik[k] = log(pi[k]) + log_mult_k + 
        bernoulli_lpmf(y[i] | Phi(to_row_vector(V_k[i, ,k]) * to_vector(xi)));
    }
    target += log_sum_exp(log_lik);
  }
}
generated quantities {
  real log_lik[K];
}
