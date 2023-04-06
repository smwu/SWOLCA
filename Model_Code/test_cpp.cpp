#include <RcppArmadillo.h>
#include <RcppTN.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppTN)]]

// Single draw from Categorical distribution
// [[Rcpp::export]]
int rcat_cpp(rowvec probs) {
  int num_categs = probs.size();
  IntegerVector draw(num_categs);
  rmultinom(1, probs.begin(), num_categs, draw.begin());
  int categ = which_max(draw);
  return categ; 
}

// log-sum-exp trick from https://github.com/helske/seqHMM/blob/main/src/logSumExp.cpp
// [[Rcpp::export]]
double logSumExp_cpp(const arma::rowvec& x) {
  int maxi = x.index_max();
  double maxv = x(maxi);
  if (!(maxv > -arma::datum::inf)) {
    return -arma::datum::inf;
  }
  double cumsum = 0.0;
  for (unsigned int i = 0; i < x.n_elem; i++) {
    if ((i != maxi) && (x(i) > -arma::datum::inf)) {
      cumsum += exp(x(i) - maxv);
    }
  }
  return maxv + log1p(cumsum);
}

// Update c
// [[Rcpp::export]]
vec update_c_OLCA(vec& c_all, const int& n, const int& K, const int& p, 
             const cube& theta, const mat& x_mat, const vec& pi) {
  mat log_cond_c(n, K);        // Individual log-likelihood for each class
  mat pred_class_probs(n, K);  // Posterior class membership probabilities
  
  // Calculate posterior class membership, p(c_i=k|-), for each class k and
  // update class assignments
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < K; k++) {
      // Calculate theta component of individual log-likelihood for class k
      double log_theta_comp_k = 0.0;
      for (int j = 0; j < p; j++) {
        // Subtract 1 from exposure value due to 0-based indexing
        log_theta_comp_k += log(theta(j, k, x_mat(i, j) - 1));
      }
      // Individual log-likelihood for class k
      log_cond_c(i, k) = log(pi(k)) + log_theta_comp_k;
    }
    // Calculate p(c_i=k|-) = p(x,y,c_i=k) / p(x,y)
    pred_class_probs.row(i) = exp(log_cond_c.row(i) - logSumExp_cpp(log_cond_c.row(i)));
    // Update class assignment using the posterior probabilities
    // Be careful of 0-based indexing
    c_all(i) = rcat_cpp(pred_class_probs.row(i)) + 1;
  }
  return c_all;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
library(R.matlab)
library(stringr)
library(fastDummies)
library(LaplacesDemon)
library(gtools)
set.seed(1)
wd = "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/"
data_dir = "Data/"
scen_samp = 101
iter_pop = 1
samp_n = 1
data_path = paste0(wd, data_dir, "simdata_scen", scen_samp, "_iter", iter_pop,
                    "_samp", samp_n, ".mat")   # Input dataset
data_vars = readMat(data_path)$sim.data
names(data_vars) = str_replace_all(dimnames(data_vars)[[1]], "[.]", "_")

K = 3
n = 10
p = 5
d = 4
S = 2
q = 2
x_mat = data_vars$X_data[1:n, 1:p]
y_all = data_vars$Y_data[1:n]
z_all = rnorm(n)
z_all = ifelse(y_all == 1, abs(z_all), -abs(z_all))
s_all = data_vars$true_Si[1:n]
V = as.matrix(dummy_cols(data.frame(x = factor(s_all, levels = 1:S)),
                          remove_selected_columns = TRUE))
kappa = sum(data_vars$sample_wt[1:n]) / n
w_all = c(data_vars$sample_wt[1:n] / kappa)
c_all = data_vars$true_Ci[1:n]

alpha = rep(1, 3) / 3
eta = rep(1, d)
mu0 = Sig0 = vector("list", K)
for (k in 1:K) {
  mu0[[k]] = rnorm(n = q)
  Sig0[[k]] = diag(rinvgamma(n = q, shape = 3.5, scale = 6.25), nrow = q, ncol = q)
}

pi = rdirichlet(1, alpha)
theta = data_vars$true_global_thetas[1:p, , ]
xi = matrix(data_vars$true_xi, nrow = K, ncol = q, byrow = FALSE)
loglik = numeric(n)
*/

/*** R
# print(pi)
# print(theta)
# print(c_all)
# print(z_all)
# print(x_mat)
# print(V)
# print(y_all)
# print(xi)
update_c_OLCA(c_all, n, K, p, theta, x_mat, pi)
c_all
*/
