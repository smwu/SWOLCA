#include <RcppArmadillo.h>
#include <RcppTN.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppTN)]]

// [[Rcpp::export]]
arma::field<arma::cube> create_4d_array(int p, int K, int d, int M, cube values) {
  arma::field<arma::cube> array_4d(M);
  array_4d.fill(arma::cube(p,K,d)); // fill all with cubes
  array_4d[0] = values;
  double b = array_4d[0](0,0,0);
  Rcout << "b: " << b;
  return array_4d;
}

// [[Rcpp::export]]
arma::field<arma::cube> posterior_field(int n, int L, Rcpp::NumericVector N, int TIME) {
  arma::field<arma::cube> A(TIME);
  A.fill(arma::cube(n, L, max(N), arma::fill::zeros));
  double b = A[1](1,1,1);
  Rcout << "b: " << b;
  return A;
}

// Draw from multivariate Normal distribution
// [[Rcpp::export]]
mat mvrnorm_cpp2(const int& n, const vec& mu, const mat& sigma) {
  int ncols = sigma.n_cols;
  mat z = randn(n, ncols) * chol(sigma);
  return mu.t() + z;
}

// Draw from multivariate Normal distribution
// [[Rcpp::export]]
mat mvrnorm_cpp3(const int& n, const rowvec& mu, const mat& sigma) {
  Environment pkg = Environment::namespace_env("LaplacesDemon");
  Function f("rmvn");
  int ncols = sigma.n_cols;
  // NumericMatrix temp = f(_["n"] = 1, _["mu"] = mu, _["Sigma"] = sigma);
  mat temp = as<arma::mat>(f(1, _["mu"] = mu, _["Sigma"] = sigma));
  return temp;
}

// Draw from truncated Normal distribution
// [[Rcpp::export]]
double rtruncnorm_cpp(const int& n, const double& a, const double& b,
                      const double& mean, const double& sd) {
  Environment pkg = Environment::namespace_env("truncnorm");
  Function f("rtruncnorm");
  double rtn = as<double>(f(_["n"] = n, _["a"] = a, _["b"] = b, _["mean"] = mean, 
                            _["sd"] = sd));
  return rtn;
}

// Draw from truncated Normal distribution
// [[Rcpp::export]]
double rtn1_cpp(const double& lin_pred, const double& sd, const double& a,
                const double& b) {
  return RcppTN::rtn1(lin_pred, sd, a, b);
}

/*** R
# # posterior_field(3, 4, c(1,2,3), 10)
# values <- array(1:6, dim=c(2,3,1))
# arr <- create_4d_array(2,3,1,5, values)
# dim(arr)
# arr
# arr[1,1]
# arr[1,1][[1]][1,,]
# mu0 <- Sig0 <- vector("list", 3)
# for (k in 1:3) {
#   # MVN(0,1) hyperprior for prior mean of xi
#   mu0[[k]] <- rnorm(n = 2)
#   # InvGamma(3.5, 6.25) hyperprior for prior variance of xi. Assume uncorrelated 
#   # components and mean variance 2.5 for a weakly informative prior on xi
#   Sig0[[k]] <- diag(rinvgamma(n = 2, shape = 3.5, scale = 6.25), nrow = 2, ncol = 2)
# }
# library(LaplacesDemon)
# set.seed(1)
# mvrnorm_cpp3(1, mu0[[1]], Sig0[[1]])
# set.seed(1)
# rmvn(1, mu0[[1]], Sig0[[1]])
set.seed(1)
rtruncnorm_cpp(1, 0, Inf, 0.6, 1)
rtruncnorm_cpp(1, -Inf, 0, 0.6, 1)
set.seed(1)
rtn1_cpp(0.6, 1, 0, Inf)
rtn1_cpp(0.6, 1, -Inf, 0)
library(truncnorm)
set.seed(1)
rtruncnorm(1, 0, Inf, 0.6, 1)
rtruncnorm(1, -Inf, 0, 0.6, 1)

# library(microbenchmark)
# microbenchmark(
#   rtruncnorm_cpp(1, 0, Inf, 0.6, 1),
#   rtn1_cpp(0.6, 1, 0, Inf)
#   # mvrnorm_cpp3(1, mu0[[1]], Sig0[[1]]),
#   # mvrnorm_cpp2(1, mu0[[1]], Sig0[[1]]),
#   # rmvn(1, mu0[[1]], Sig0[[1]])
# )
*/