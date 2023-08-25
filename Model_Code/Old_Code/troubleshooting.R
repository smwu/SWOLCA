### Test R vs. Rcpp implementations

library(plyr)
library(dplyr)
library(LaplacesDemon)
library(truncnorm)
library(stringr)
library(R.matlab)
library(fastDummies)
library(matrixStats)
library(Matrix)
library(gtools)
library(e1071)
library(rstan)
library(survey)

set.seed(1)
scen_samp <- 111211
iter_pop <- 1
samp_n <- 3
# Define directories
wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/"
data_dir <- "Data/"
res_dir <- "Results/"
model_dir <- "Model_Code/"
model <- "wsOFMM"

# Define paths
# data_path <- paste0(wd, data_dir, "simdata_scen", scen_samp, "_iter", iter_pop,
#                     "_samp", samp_n, ".mat")   # Input dataset
data_path <- paste0(wd, data_dir, "simdata_scen", scen_samp, "_iter", iter_pop,
                    "_samp", samp_n, ".RData")  # Input dataset
res_path <- paste0(wd, res_dir, model, "_results_R20000_scen", scen_samp, 
                   "_samp", samp_n, ".RData")  # Output file
adj_path <- paste0(wd, res_dir, model, "_results_R20000_adjR_scen", scen_samp, 
                   "_samp", samp_n, ".RData")      # Adjusted output file
stan_path <- paste0(wd, model_dir, "WSOLCA_main.stan")  # Stan file
covs <- "true_Si"

load(data_path)
data_vars <- sim_data

# Obtain dimensions
n <- dim(data_vars$X_data)[1]        # Number of individuals
p <- dim(data_vars$X_data)[2]        # Number of exposure items
d <- max(apply(data_vars$X_data, 2,  # Number of exposure categories 
               function(x) length(unique(x))))  # CHANGE TO ADAPT TO ITEM
# Obtain data
x_mat <- data_vars$X_data            # Categorical exposure matrix, nxp
y_all <- c(data_vars$Y_data)         # Binary outcome vector, nx1

#======= R
s_all <- data_vars$true_Si           # Stratifying variable, nx1
S <- length(unique(s_all))           # Number of strata
s_mat <- dummy_cols(data.frame(s = factor(s_all)),  # Stratifying variable as dummy columns
                    remove_selected_columns = TRUE)
q1 <- S                               # Number of regression covariates excluding class assignment
V1 <- as.matrix(s_mat)                # Regression design matrix without class assignment, nxq

#======= Rcpp
s_all <- data_vars[[covs]]    # Stratifying variable, nx1
# Stratifying variable as dummy columns
s_mat <- dummy_cols(data.frame(s = factor(s_all)),  
                    remove_selected_columns = TRUE)
V2 <- as.matrix(s_mat)              # Regression design matrix without class assignment, nxq
q2 <- ncol(V2)                       # Number of regression covariates excluding class assignment

identical(q1, q2)
identical(V1, V2)

#======================= PASSED!

V <- V1
q <- q1

# Obtain normalized weights
kappa <- sum(data_vars$sample_wt) / n   # Weights norm. constant. If sum(weights)=N, this is 1/(sampl_frac)
w_all <- c(data_vars$sample_wt / kappa) # Weights normalized to sum to n, nx1

K_fixed <- 3
eta <- rep(1, d)                 # Hyperparameter for prior for theta


#====== R
# 'init_OLCA' initializes priors and variables for the OLCA model given hyperparameters
# Inputs:
#   alpha: Vector of hyperparameters for pi. Kx1
#   eta: Vector of hyperparameters for theta. dx1
#   n: Number of individuals
#   K: Number of classes
#   p: Number of exposure items
#   d: Number of exposure categories
# Outputs: `OLCA_params` list with the following items:
#   pi: Vector parameter pi for class membership probabilities. Kx1
#   theta: Array parameter theta for item category probabilities. pxKxd
#   c_all: Vector of random initial class assignments. nx1
init_OLCA <- function(alpha, eta, n, K, p, d) {
  # Prior for pi
  pi <- c(rdirichlet(n = 1, alpha = alpha))
  
  # Initialize class assignment, c, for individuals
  c_all <- rcat(n = n, p = pi)
  
  # Prior for theta
  theta <- array(0, dim = c(p, K, d))
  for (j in 1:p) {
    for (k in 1:K) {
      theta[j, k, ] <- c(rdirichlet(n = 1, alpha = eta))
    }
  }
  
  OLCA_params <- list(pi = pi, theta = theta, c_all = c_all)
  return(OLCA_params)
}

# `init_probit` initializes priors and variables for the probit regression model
# given hyperparameters 
# Inputs:
#   mu0: List of vectors of mean hyperparameters for xi. List of K qx1 vectors
#   Sig0: List of matrices of variance hyperparameters for xi. List of K qxq matrices
#   K: Number of classes
#   q: Number of regression covariates excluding class assignment
#   n: Number of individuals
#   V: Regression design matrix without class assignment. nxq
#   y_all: Vector of outcomes. nx1
#   c_all: Vector of random initial class assignments. nx1
# Outputs: `probit_params` list with the following items:
#   xi: Matrix parameter xi for probit regression coefficients. Kxq
#   z_all: Vector of latent variables in the probit model. nx1
init_probit <- function(mu0, Sig0, K, q, n, V, y_all, c_all) {
  # Initialize variables
  xi <- matrix(NA, nrow = K, ncol = q)
  z_all <- lin_pred <- numeric(n)
  
  # Prior for xi. Same prior for each class
  for (k in 1:K) {
    xi[k, ] <- rmvn(n = 1, mu = mu0[[k]], Sigma = Sig0[[k]])
  }
  
  # Initialize probit model latent variable, z, following Albert and Chib (1993)
  # Linear predictor using covariate values and class assignment for each individual
  for (i in 1:n) {
    lin_pred[i] <- V[i, ] %*% xi[c_all[i], ]  # V[i]*xi[c_i]
  }
  # For cases, z_i ~ TruncNormal(low=0, high=Inf, mean=V*xi[c_i], var=1)
  z_all[y_all == 1] <- rtruncnorm(n = 1, a = 0, b = Inf, 
                                  mean = lin_pred[y_all == 1], sd = 1)
  # For controls, z_i ~ TruncNormal(low=-Inf, high=0, mean=V*xi[c_i], var=1)
  z_all[y_all == 0] <- rtruncnorm(n = 1, a = -Inf, b = 0, 
                                  mean = lin_pred[y_all == 0], sd = 1)
  # Control extreme values
  z_all[z_all == Inf] <- qnorm(1 - (1e-10))
  z_all[z_all == -Inf] <- qnorm(1e-10)
  
  probit_params <- list(xi = xi, z_all = z_all)
  return(probit_params)
}
set.seed(1)

# Initialize OLCA model using fixed number of classes
alpha <- rep(1, K_fixed) / K_fixed  # Hyperparameter for prior for pi
# Obtain pi, theta, c_all
OLCA_params1 <- init_OLCA(alpha = alpha, eta = eta, n = n, K = K_fixed, p = p, 
                         d = d)

# Initialize probit model using fixed number of classes
# Initialize hyperparameters for xi
mu0 <- Sig0 <- vector("list", K_fixed)
for (k in 1:K_fixed) {
  # MVN(0,1) hyperprior for prior mean of xi
  mu0[[k]] <- rnorm(n = q)
  # InvGamma(3.5, 6.25) hyperprior for prior variance of xi. Assume uncorrelated 
  # components and mean variance 2.5 for a weakly informative prior on xi
  Sig0[[k]] <- diag(rinvgamma(n = q, shape = 3.5, scale = 6.25), nrow = q, ncol = q)
}
# Obtain xi, z_all                        
probit_params1 <- init_probit(mu0 = mu0, Sig0 = Sig0, K = K_fixed, q = q, n = n, 
                             V = V, y_all = y_all, c_all = OLCA_params1$c_all)

#====== Rcpp
# Source R helper functions
source(paste0(wd, model_dir, "helper_functions.R"))
# Source Rcpp functions
Rcpp::sourceCpp(paste0(wd, model_dir, "main_Rcpp_functions.cpp"))

set.seed(1)
# Initialize OLCA model using fixed number of classes
# alpha <- rep(2, K_fixed) # Hyperparameter for prior for pi
  alpha <- rep(1, K_fixed) / K_fixed  # Hyperparameter for prior for pi
# Obtain pi, theta, c_all
OLCA_params2 <- init_OLCA(alpha = alpha, eta = eta, n = n, K = K_fixed, p = p, 
                         d = d)

# Initialize probit model using fixed number of classes
# Initialize hyperparameters for xi
mu0 <- Sig0 <- vector("list", K_fixed)
for (k in 1:K_fixed) {
  # MVN(0,1) hyperprior for prior mean of xi
  mu0[[k]] <- rnorm(n = q)
  # InvGamma(3.5, 6.25) hyperprior for prior variance of xi. Assume uncorrelated 
  # components and mean variance 2.5 for a weakly informative prior on xi
  Sig0[[k]] <- diag(rinvgamma(n = q, shape = 3.5, scale = 6.25), nrow = q, ncol = q)
}
# Obtain xi, z_all                        
probit_params2 <- init_probit(mu0 = mu0, Sig0 = Sig0, K = K_fixed, q = q, n = n, 
                             V = V, y_all = y_all, c_all = OLCA_params2$c_all)

identical(OLCA_params1, OLCA_params2)
identical(probit_params1, probit_params2)

#======================= PASSED!

OLCA_params <- OLCA_params1
probit_params <- probit_params1

n_runs <- 200
burn <- 100
thin <- 5

#===== R
run_MCMC <- function(OLCA_params, probit_params, n_runs, burn, thin, K, p, d, n, 
                     q, w_all, x_mat, y_all, V, alpha, eta, mu0, Sig0) {
  n_storage <- ceiling(n_runs / thin)  # Number of MCMC iterations to store
  
  # Initialize variables
  pi_MCMC <- matrix(NA, nrow = n_storage, ncol = K)
  theta_MCMC <- array(NA, dim = c(n_storage, p, K, d))
  xi_MCMC <- array(NA, dim = c(n_storage, K, q))
  c_all_MCMC <- z_all_MCMC <- matrix(NA, nrow = n_storage, ncol = n)
  loglik_MCMC <- matrix(NA, nrow = n_storage, ncol = n)
  
  # Initialized values
  pi <- OLCA_params$pi
  theta <- OLCA_params$theta
  c_all <- OLCA_params$c_all
  xi <- probit_params$xi
  z_all <- probit_params$z_all
  
  alpha_post <- numeric(K)                            # Posterior parameters for pi
  log_cond_c <- matrix(NA, nrow = n, ncol = K)        # Individual log-likelihood for each class
  pred_class_probs <- matrix(NA, nrow = n, ncol = K)  # Posterior class membership probabilities
  eta_post <- numeric(d)                              # Posterior parameters for theta
  loglik <- numeric(n)                                # Individual log-likelihood
  lin_pred <- numeric(n)                              # Linear predictor, V*xi
  
  # Update parameters and variables
  for (m in 1:n_runs) {
    #============== Update pi ==================================================
    for (k in 1:K) {
      # Add sum of normalized weights for those assigned to class k, equiv. to
      # weighted number of individuals assigned to each class
      alpha_post[k] <- alpha[k] + sum(w_all[c_all == k])
    }
    pi <- c(rdirichlet(n = 1, alpha = alpha_post))
    
    #============== Update c ===================================================
    # Calculate posterior class membership, p(c_i=k|-), for each class k and 
    # update class assignments
    for (i in 1:n) {
      for (k in 1:K) {
        # Calculate theta component of individual log-likelihood for class k
        log_theta_comp_k <- 0
        for (j in 1:p) {
          log_theta_comp_k <- log_theta_comp_k + log(theta[j, k, x_mat[i, j]])
        }
        # Individual log-likelihood for class k
        log_cond_c[i, k] <- log(pi[k]) + log_theta_comp_k +
          log(dnorm(z_all[i], mean = V[i, ] %*% xi[k, ])) +  # probit component
          log(y_all[i]*(z_all[i] > 0) + (1 - y_all[i])*(z_all[i] <= 0))
      }
      # Calculate p(c_i=k|-) = p(x,y,c_i=k) / p(x,y)
      pred_class_probs[i, ] <- exp(log_cond_c[i, ] - logSumExp(log_cond_c[i, ]))
      # Update class assignment using the posterior probabilities
      c_all[i] <- rcat(n = 1, p = pred_class_probs[i, ])
    }
    
    #============== Update theta ===============================================
    for (j in 1:p) {
      for (k in 1:K) {
        for (r in 1:d) {
          # Add sum of normalized weights for those assigned to class k with x_ij = r
          eta_post[r] <- eta[r] + sum(w_all[(c_all == k) & (x_mat[,j] == r)])
        }
        theta[j, k, ] <- c(rdirichlet(n = 1, alpha = eta_post))
      }
    }
    
    #============== Update xi ==================================================
    W_tilde <- sparseMatrix(i = 1:n, j = 1:n, x = w_all)  # Sparse diagonal normalized weight matrix
    V_sparse <- as(V, "sparseMatrix")                     # Sparse design matrix without class assignments
    for (k in 1:K) {
      # Sparse diagonal matrix subsetting to obs in class k
      C_k <- sparseMatrix(i = 1:n, j = 1:n, x = (c_all == k))
      
      # Draw xi from conditional posterior distribution
      Sig_post <- solve(Sig0[[k]]) + as.matrix(t(V) %*% C_k %*% W_tilde %*% V)
      mu_right <- (solve(Sig0[[k]]) %*% mu0[[k]]) + 
        (as.matrix(t(V) %*% C_k %*% W_tilde %*% z_all))
      mu_post <- c(solve(Sig_post) %*% mu_right)
      xi[k, ] <- rmvn(n = 1, mu = mu_post, Sigma = solve(Sig_post))
    }
    
    #============== Update z ===================================================
    # Linear predictor using covariate values and class assignment for each individual
    for (i in 1:n) {
      lin_pred[i] <- V[i, ] %*% xi[c_all[i], ]  # V[i]*xi[c_i]
    }
    # Probit model latent variable z, following Albert and Chib (1993)
    # For cases, z_i ~ TruncNormal(low=0, high=Inf, mean=V*xi[c_i], var=1)
    z_all[y_all == 1] <- rtruncnorm(n = 1, a = 0, b = Inf, 
                                    mean = lin_pred[y_all == 1], sd = 1)
    # For controls, z_i ~ TruncNormal(low=-Inf, high=0, mean=V*xi[c_i], var=1)
    z_all[y_all == 0] <- rtruncnorm(n = 1, a = -Inf, b = 0, 
                                    mean = lin_pred[y_all == 0], sd = 1)
    # Control extreme values
    z_all[z_all == Inf] <- qnorm(1 - (1e-10))
    z_all[z_all == -Inf] <- qnorm(1e-10)
    
    
    #============== Update individual log-likelihood  ==========================
    for (i in 1:n) {
      c_i <- c_all[i]
      # Calculate theta component of individual log-likelihood
      log_theta_comp <- 0
      for (j in 1:p) {
        log_theta_comp <- log_theta_comp + log(theta[j, c_i, x_mat[i, j]])
      }
      loglik[i] <- log(pi[c_i]) + log_theta_comp +
        log(dnorm(z_all[i], mean = V[i, ] %*% xi[c_i, ])) + 
        log(y_all[i]*(z_all[i] > 0) + (1 - y_all[i])*(z_all[i] <= 0))
    }
    
    #============== Store posterior values based on thinning  ==================
    if (m %% thin == 0) {
      m_thin <- m / thin
      pi_MCMC[m_thin, ] <- pi
      theta_MCMC[m_thin, , , ] <- theta
      xi_MCMC[m_thin, , ] <- xi
      c_all_MCMC[m_thin, ] <- c_all
      z_all_MCMC[m_thin, ] <- z_all
      loglik_MCMC[m_thin] <- sum(loglik)
    }
    
    #============== Relabel classes every 10 iterations to encourage mixing ====
    if (m %% 10 == 0) {
      new_order <- permute(1:K)      # New permuted labels
      new_c_all <- numeric(K)        # Class assignments with new labels
      for (k in 1:K) {
        new_c_all[c_all == k] <- new_order[k]
      }
      c_all <- new_c_all             # Relabel class assignments
      pi <- pi[new_order]            # Relabel class probabilities
      theta <- theta[, new_order, ]  # Relabel item category probabilities
      xi <- xi[new_order, ]          # Relabel probit coefficients
      print(paste0("Iteration ", m, " completed!"))
    }
  }
  
  # Discard burn-in
  warmup <- ceiling(burn / thin)
  pi_MCMC <- pi_MCMC[-(1:warmup), ]
  theta_MCMC <- theta_MCMC[-(1:warmup), , , ]
  xi_MCMC <- xi_MCMC[-(1:warmup), , ]
  c_all_MCMC <- c_all_MCMC[-(1:warmup), ]
  z_all_MCMC <- z_all_MCMC[-(1:warmup), ]
  loglik_MCMC <- loglik_MCMC[-(1:warmup), ]
  
  MCMC_out <- list(pi_MCMC = pi_MCMC, theta_MCMC = theta_MCMC, xi_MCMC = xi_MCMC,
                   c_all_MCMC = c_all_MCMC, z_all_MCMC = z_all_MCMC, 
                   loglik_MCMC = loglik_MCMC)
  return(MCMC_out)
}

set.seed(1)
MCMC_out1 <- run_MCMC(OLCA_params = OLCA_params, probit_params = probit_params, 
                     n_runs = n_runs, burn = burn, thin = thin, K = K_fixed, 
                     p = p, d = d, n = n, q = q, w_all = w_all, x_mat = x_mat, 
                     y_all = y_all, V = V, alpha = alpha, eta = eta, 
                     Sig0 = Sig0, mu0 = mu0)

#===== Rcpp
# Source R helper functions
source(paste0(wd, model_dir, "helper_functions.R"))
# Source Rcpp functions
Rcpp::sourceCpp(paste0(wd, model_dir, "main_Rcpp_functions.cpp"))

set.seed(1)
MCMC_out2 <- run_MCMC_Rcpp(OLCA_params = OLCA_params, probit_params = probit_params, 
                          n_runs = n_runs, burn = burn, thin = thin, K = K_fixed, 
                          p = p, d = d, n = n, q = q, w_all = w_all, x_mat = x_mat, 
                          y_all = y_all, V = V, alpha = alpha, eta = eta, 
                          Sig0 = Sig0, mu0 = mu0)

identical(MCMC_out1, MCMC_out2)

#==================== FAILED

K <- K_fixed

#==================== Unpacking run_MCMC
n_storage <- ceiling(n_runs / thin)  # Number of MCMC iterations to store

# Initialize variables
pi_MCMC <- matrix(NA, nrow = n_storage, ncol = K)
theta_MCMC <- array(NA, dim = c(n_storage, p, K, d))
xi_MCMC <- array(NA, dim = c(n_storage, K, q))
c_all_MCMC <- z_all_MCMC <- matrix(NA, nrow = n_storage, ncol = n)
loglik_MCMC <- numeric(n_storage)
loglik <- numeric(n)  # Individual log-likelihood

# Initialized values
pi <- OLCA_params$pi
theta <- OLCA_params$theta
c_all <- as.double(OLCA_params$c_all)  # allows updating by reference in rcpparmadillo
xi <- probit_params$xi
z_all <- probit_params$z_all

#============== Update pi ==================================================
#====== R
set.seed(1)
alpha_post <- numeric(K)                            # Posterior parameters for pi
log_cond_c <- matrix(NA, nrow = n, ncol = K)        # Individual log-likelihood for each class
pred_class_probs <- matrix(NA, nrow = n, ncol = K)  # Posterior class membership probabilities
eta_post <- numeric(d)                              # Posterior parameters for theta
loglik <- numeric(n)                                # Individual log-likelihood
lin_pred <- numeric(n)                              # Linear predictor, V*xi
for (k in 1:K) {
  # Add sum of normalized weights for those assigned to class k, equiv. to
  # weighted number of individuals assigned to each class
  alpha_post[k] <- alpha[k] + sum(w_all[c_all == k])
}
pi1 <- c(rdirichlet(n = 1, alpha = alpha_post))
pi1

#====== Rcpp
set.seed(1)
update_pi(pi = pi, w_all = w_all, c_all = c_all, K = K, alpha = alpha)
pi

identical(pi1, pi)
#=============== PASSED!

pi <- pi1

#============== Update c ===================================================
#====== R
set.seed(1)
c_all1 <- c_all
# Calculate posterior class membership, p(c_i=k|-), for each class k and 
# update class assignments
for (i in 1:n) {
  for (k in 1:K) {
    # Calculate theta component of individual log-likelihood for class k
    log_theta_comp_k <- 0
    for (j in 1:p) {
      log_theta_comp_k <- log_theta_comp_k + log(theta[j, k, x_mat[i, j]])
    }
    # Individual log-likelihood for class k
    log_cond_c[i, k] <- log(pi[k]) + log_theta_comp_k +
      log(dnorm(z_all[i], mean = V[i, ] %*% xi[k, ])) +  # probit component
      log(y_all[i]*(z_all[i] > 0) + (1 - y_all[i])*(z_all[i] <= 0))
  }
  # Calculate p(c_i=k|-) = p(x,y,c_i=k) / p(x,y)
  pred_class_probs[i, ] <- exp(log_cond_c[i, ] - logSumExp(log_cond_c[i, ]))
  # Update class assignment using the posterior probabilities
  c_all1[i] <- rcat(n = 1, p = pred_class_probs[i, ])
}

#====== Rcpp
set.seed(1)
update_c(c_all = c_all, n = n, K = K, p = p, theta = theta, 
         x_mat = x_mat, pi = pi, z_all = z_all, V = V, xi = xi, 
         y_all = y_all)

identical(c_all1, c_all)
#=============== PASSED!
c_all <- c_all1

#============== Update theta ===============================================
#===== R
set.seed(1)
theta1 <- theta

for (j in 1:p) {
  for (k in 1:K) {
    for (r in 1:d) {
      # Add sum of normalized weights for those assigned to class k with x_ij = r
      eta_post[r] <- eta[r] + sum(w_all[(c_all == k) & (x_mat[,j] == r)])
    }
    theta1[j, k, ] <- c(rdirichlet(n = 1, alpha = eta_post))
  }
}

#====== Rcpp
set.seed(1)
update_theta(theta = theta, p = p, K = K, d = d, eta = eta, 
             w_all = w_all, c_all = c_all, x_mat = x_mat)

identical(theta1, theta)
#=============== PASSED!

theta <- theta1

#============== Update xi ==================================================
#====== R
set.seed(1)
xi1 <- xi
W_tilde <- sparseMatrix(i = 1:n, j = 1:n, x = w_all)  # Sparse diagonal normalized weight matrix
V_sparse <- as(V, "sparseMatrix")                     # Sparse design matrix without class assignments
for (k in 1:K) {
  # Sparse diagonal matrix subsetting to obs in class k
  C_k <- sparseMatrix(i = 1:n, j = 1:n, x = (c_all == k))
  
  # Draw xi from conditional posterior distribution
  Sig_post <- solve(Sig0[[k]]) + as.matrix(t(V) %*% C_k %*% W_tilde %*% V)
  mu_right <- (solve(Sig0[[k]]) %*% mu0[[k]]) + 
    (as.matrix(t(V) %*% C_k %*% W_tilde %*% z_all))
  mu_post <- c(solve(Sig_post) %*% mu_right)
  xi1[k, ] <- rmvn(n = 1, mu = mu_post, Sigma = solve(Sig_post))
}

#====== Rcpp
set.seed(1)
update_xi(xi = xi, n = n, K = K, w_all = w_all, c_all = c_all, 
          z_all = z_all, V = V, y_all = y_all, mu0 = mu0, Sig0 = Sig0)

all.equal(xi1, xi)
#================ FAILED
#================ Passed after troubleshooting
# Changed to call rmvn function from Rcpp
Rcpp::sourceCpp(paste0(wd, model_dir, "main_Rcpp_functions.cpp"))
# Now all equal up to machine precision

xi <- xi1

#============== Update z ===================================================
#===== R
set.seed(3)
z_all1 <- z_all
# Linear predictor using covariate values and class assignment for each individual
for (i in 1:n) {
  lin_pred[i] <- V[i, ] %*% xi[c_all[i], ]  # V[i]*xi[c_i]
}
# Probit model latent variable z, following Albert and Chib (1993)
# For cases, z_i ~ TruncNormal(low=0, high=Inf, mean=V*xi[c_i], var=1)
z_all1[y_all == 1] <- rtruncnorm(n = 1, a = 0, b = Inf, 
                                mean = lin_pred[y_all == 1], sd = 1)
# For controls, z_i ~ TruncNormal(low=-Inf, high=0, mean=V*xi[c_i], var=1)
z_all1[y_all == 0] <- rtruncnorm(n = 1, a = -Inf, b = 0, 
                                mean = lin_pred[y_all == 0], sd = 1)
# Control extreme values
z_all1[z_all1 == Inf] <- qnorm(1 - (1e-10))
z_all1[z_all1 == -Inf] <- qnorm(1e-10)

#===== Rcpp
set.seed(3)
update_z(z_all = z_all, n = n, V = V, xi = xi, c_all = c_all, 
         y_all = y_all)

identical(z_all1, z_all)
#================ FAILED
#================ Can't get identical but results somewhat similar
# Changed to call rtruncnorm function from Rcpp
Rcpp::sourceCpp(paste0(wd, model_dir, "main_Rcpp_functions.cpp"))
# Now all equal up to machine precision

z_all <- z_all1

#============== Update individual log-likelihood  ==========================
#====== R 
set.seed(1)
loglik1 <- loglik

for (i in 1:n) {
  c_i <- c_all[i]
  # Calculate theta component of individual log-likelihood
  log_theta_comp <- 0
  for (j in 1:p) {
    log_theta_comp <- log_theta_comp + log(theta[j, c_i, x_mat[i, j]])
  }
  loglik1[i] <- log(pi[c_i]) + log_theta_comp +
    log(dnorm(z_all[i], mean = V[i, ] %*% xi[c_i, ])) + 
    log(y_all[i]*(z_all[i] > 0) + (1 - y_all[i])*(z_all[i] <= 0))
}

#====== Rcpp
set.seed(1)
update_loglik(loglik = loglik, n = n, p = p, c_all = c_all, 
              theta = theta, x_mat = x_mat, pi = pi, 
              z_all = z_all, V = V, xi = xi, y_all = y_all)

identical(loglik1, loglik)

#================ PASSED!

loglik <- loglik1

MCMC_out <- MCMC_out1

#====== R
post_process <- function(MCMC_out, p, d, q) {
  # Get median number of classes with >= 5% of individuals, over all iterations
  M <- dim(MCMC_out$pi_MCMC)[1]  # Number of stored MCMC iterations
  K_med <- round(median(rowSums(MCMC_out$pi_MCMC >= 0.05)))
  
  # Cluster individuals into reduced number of classes using agglomerative clustering
  # Calculate pairwise distance matrix using Hamming distance: proportion of 
  # iterations where two individuals have differing class assignments
  distMat <- hamming.distance(t(MCMC_out$c_all_MCMC))
  dendrogram <- hclust(as.dist(distMat), method = "complete") # Hierarchical clustering dendrogram
  red_c_all <- cutree(dendrogram, k = K_med)                  # Group individuals into K_med classes
  # For each iteration, relabel new classes using the most common old class assignment
  relabel_red_classes <- matrix(NA, nrow = M, ncol = K_med)   # Apply new classes to each iteration
  for (k in 1:K_med) {
    relabel_red_classes[, k] <- apply(as.matrix(MCMC_out$c_all_MCMC[, red_c_all == k]), 
                                      1, get_mode)
  }
  
  # Reduce and reorder parameter estimates using new classes
  pi <- matrix(NA, nrow = M, ncol = K_med)
  theta <- array(NA, dim = c(M, p, K_med, d))
  xi <- array(NA, dim = c(M, K_med, q))
  for (m in 1:M) {
    iter_order <- relabel_red_classes[m, ]
    pi_temp <- MCMC_out$pi_MCMC[m, iter_order]
    pi[m, ] <- pi_temp / sum(pi_temp)
    theta[m, , , ] <- MCMC_out$theta_MCMC[m, , iter_order, ]
    xi[m, , ] <- MCMC_out$xi_MCMC[m, iter_order, ]
  }
  
  post_MCMC_out <- list(K_med = K_med, pi = pi, theta = theta, xi = xi)
  
  return(post_MCMC_out)  
  # plot(dendrogram, labels = FALSE)
  # rect.hclust(dendrogram, k = K_med)
}

set.seed(1)
post_MCMC_out1 <- post_process(MCMC_out = MCMC_out, p = p, d = d, q = q)

#====== Rcpp
source(paste0(wd, model_dir, "helper_functions.R"))
set.seed(1)
post_MCMC_out2 <- post_process(MCMC_out = MCMC_out, p = p, d = d, q = q)

identical(post_MCMC_out1, post_MCMC_out2)

#================ PASSED!

post_MCMC_out <- post_MCMC_out1

#===== R
analyze_results <- function(MCMC_out, post_MCMC_out, n, p, V, y_all, x_mat) {
  
  #============== Identify unique classes using modal exposure categories ======
  # Posterior median estimate for theta across iterations
  theta_med_temp <- apply(post_MCMC_out$theta, c(2, 3, 4), median)
  # Posterior modal exposure categories for each exposure item and reduced class
  theta_modes <- apply(theta_med_temp, c(1, 2), which.max)
  # Identify unique classes
  unique_classes <- which(!duplicated(theta_modes, MARGIN = 2))
  # Number of unique classes
  K_red <- length(unique_classes) 
  
  #============== Use new classes to adjust and re-normalize posterior samples =
  # Combine duplicated classes and re-normalize pi to sum to 1
  M <- dim(post_MCMC_out$pi)[1]                # Number of iterations
  pi_red <- post_MCMC_out$pi[, unique_classes] # Initialize pi for unique classes
  if (K_red < dim(post_MCMC_out$pi)[2]) {  # Check if there are duplicated classes
    for (k in 1:K_red) {
      # Find duplicated modal theta patterns
      dupes_k <- apply(theta_modes, 2, function(x)  
        identical(x,theta_modes[, unique_classes[k]]))
      # Combine class proportions for all duplicated patterns together
      pi_red[, k] <- apply(as.matrix(post_MCMC_out$pi[, dupes_k]), 1, sum)  
    }
  }
  # Re-normalize to ensure pi sums to 1 for each iteration
  pi_red = pi_red / rowSums(pi_red)  
  
  # Get posterior parameter samples for unique classes for theta and xi
  theta_red <- post_MCMC_out$theta[, , unique_classes, ]
  theta_red <- aaply(theta_red, c(1, 2, 3), function(x) x / sum(x)) # Re-normalize
  xi_red <- post_MCMC_out$xi[, unique_classes, ]
  
  #============== Posterior median estimates ===================================
  pi_med <- apply(pi_red, 2, median, na.rm = TRUE)
  pi_med <- pi_med / sum(pi_med)  # Re-normalize
  theta_med <- apply(theta_red, c(2, 3, 4), median, na.rm = TRUE)
  theta_med <- aaply(theta_med, c(1, 2), function(x) x / sum(x))  # Re-normalize
  xi_med <- apply(xi_red, c(2, 3), median, na.rm = TRUE)
  
  #============== Update c using unique classes and posterior estimates ========
  z_all <- MCMC_out$z_all_MCMC[M, ]
  c_all <- MCMC_out$c_all_MCMC[M, ]                      # Final class assignments
  pred_class_probs <- matrix(NA, nrow = n, ncol = K_red) # Posterior class membership probabilities
  log_cond_c <- matrix(NA, nrow = n, ncol = K_red)       # Individual log-likelihood for each class
  Phi_med_all_c <- pnorm(V %*% t(xi_med))  # Outcome probabilities for all classes
  Phi_med <- numeric(n)                    # Initialize individual outcome probabilities
  # Calculate posterior class membership, p(c_i=k|-), for each class k
  for (i in 1:n) {
    for (k in 1:K_red) {
      # Calculate theta component of individual log-likelihood assuming class k
      log_theta_comp_k <- 0
      for (j in 1:p) {
        log_theta_comp_k <- log_theta_comp_k + log(theta_med[j, k, x_mat[i, j]])
      }
      # Individual log-likelihood for class k
      log_cond_c[i, k] <- log(pi_med[k]) + log_theta_comp_k +
        log(dnorm(z_all[i], mean = V[i, ] %*% xi_med[k, ])) +  # probit component
        log(y_all[i]*(z_all[i] > 0) + (1 - y_all[i])*(z_all[i] <= 0))
    }
    # Calculate p(c_i=k|-) = p(x,y,c_i=k) / p(x,y)
    pred_class_probs[i, ] <- exp(log_cond_c[i, ] - logSumExp(log_cond_c[i, ]))
    # Update class assignment using the posterior probabilities
    c_all[i] <- rcat(n = 1, p = pred_class_probs[i, ])
    # Calculate outcome probabilities P(Y=1|-) using updated class assignment
    Phi_med[i] <- Phi_med_all_c[i, c_all[i]] 
  }
  
  #============== Update individual log-likelihood  ============================
  loglik_med <- numeric(n)  # Individual log-likelihood
  for (i in 1:n) {
    c_i <- c_all[i]
    # Calculate theta component of individual log-likelihood
    log_theta_comp <- 0
    for (j in 1:p) {
      log_theta_comp <- log_theta_comp + log(theta_med[j, c_i, x_mat[i, j]])
    }
    # Calculate individual log-likelihood using median estimates
    loglik_med[i] <- log(pi_med[c_i]) + log_theta_comp +
      log(dnorm(z_all[i], mean = V[i, ] %*% xi_med[c_i, ])) + 
      log(y_all[i]*(z_all[i] > 0) + (1 - y_all[i])*(z_all[i] <= 0))
  }
  
  analysis <- list(K_red = K_red, pi_red = pi_red, theta_red = theta_red, 
                   xi_red = xi_red, pi_med = pi_med, theta_med = theta_med, 
                   xi_med = xi_med, Phi_med = Phi_med, c_all = c_all, 
                   pred_class_probs = pred_class_probs, loglik_med = loglik_med)
  return(analysis)
}

set.seed(1)
analysis1 <- analyze_results(MCMC_out = MCMC_out, post_MCMC_out = post_MCMC_out, 
                            n = n, p = p, V = V, y_all = y_all, x_mat = x_mat)

#===== Rcpp
source(paste0(wd, model_dir, "helper_functions.R"))

set.seed(1)
analysis2 <- analyze_results(MCMC_out = MCMC_out, post_MCMC_out = post_MCMC_out, 
                            n = n, p = p, V = V, y_all = y_all, x_mat = x_mat)

identical(analysis1, analysis2)

#================ PASSED!

analysis <- analysis1

mod_stan <- stan_model(stan_path)

#===== R
var_adjust <- function(mod_stan, analysis, K, p, d, n, q, x_mat, y_all, V, w_all, 
                       s_all, clus_id_all) {
  #=============== Run Stan model ==============================================
  # Define data for Stan model
  alpha <- rep(1, K) / K            # Hyperparameter for prior for pi
  eta <- rep(1, d)                  # Hyperparameter for prior for theta
  mu0 <- Sig0 <- vector("list", K)  # Hyperparameters for xi
  for (k in 1:K) {
    mu0[[k]] <- rnorm(n = q)
    Sig0[[k]] <- diag(rinvgamma(n = q, shape = 3.5, scale = 6.25), nrow = q, ncol = q)
  }
  data_stan <- list(K = K, p = p, d = d, n = n, q = q, X = x_mat, y = y_all, 
                    V = V, weights = w_all, alpha = alpha, eta = eta, mu0 = mu0, 
                    Sig0 = Sig0)
  
  # Stan parameters of interest
  par_stan <- c('pi', 'theta', 'xi')  # subset of parameters interested in
  
  # Run Stan model
  # Stan will pass warnings from calling 0 chains, but will still create an 
  # out_stan object for the 'grad_log_prob()' method
  out_stan <- sampling(object = mod_stan, data = data_stan, pars = par_stan,
                       chains = 0, iter = 0, refresh = 0)
  
  #=============== Convert to unconstrained parameters =========================
  # Convert params from constrained space to unconstrained space
  unc_par_hat <- unconstrain_pars(out_stan, list("pi" = analysis$pi_med,
                                                 "theta" = analysis$theta_med,
                                                 "xi" = analysis$xi_med))
  # Get posterior MCMC samples in unconstrained space for all parameters
  M <- dim(analysis$pi_red)[1]
  unc_par_samps <- lapply(1:M, unconstrain, stan_model = out_stan, K = K, 
                          pi = analysis$pi_red, theta = analysis$theta_red, 
                          xi = analysis$xi_red)
  unc_par_samps <- matrix(unlist(unc_par_samps), byrow = TRUE, nrow = M)
  
  #=============== Post-processing adjustment in unconstrained space ===========
  # Estimate Hessian
  H_hat <- -1*optimHess(unc_par_hat, gr = function(x){grad_log_prob(out_stan, x)})
  
  # Create survey design
  svy_data <- data.frame(s = s_all, 
                         x = x_mat,
                         y = y_all, 
                         wts = w_all,
                         clus = clus_id_all)
  svydes <- svydesign(ids = ~clus, strata = ~s, weights = ~wts, data = svy_data)
  # create svrepdesign
  svyrep <- as.svrepdesign(design = svydes, type = "mrbbootstrap", 
                           replicates = 100)
  
  rep_temp <- withReplicates(design = svyrep, theta = grad_par, 
                             stan_mod = mod_stan, stan_data = data_stan, 
                             par_stan = par_stan, u_pars = unc_par_hat)
  J_hat <- vcov(rep_temp)
  
  # compute adjustment
  H_inv <- solve(H_hat)
  V1 <- H_inv %*% J_hat %*% H_inv
  
  # Check for issues with negative diagonals
  if (min(diag(V1)) < 0) {
    print("V1 has negative variances")
  }
  if (min(diag(H_inv)) < 0) {
    print("H_inv has negative variances")
  }
  # If matrices are not p.d. due to rounding issues, convert to nearest p.d. 
  # matrix using method proposed in Higham (2002)
  if (min(Re(eigen(V1)$values)) < 0) { 
    V1_pd <- nearPD(V1)
    R1 <- chol(V1_pd$mat)
    print(paste0("V1: absolute eigenvalue difference to nearest p.d. matrix: ", 
                 sum(abs(eigen(V1)$values - eigen(V1_pd$mat)$values))))
  } else {
    R1 <- chol(V1)
  }
  if (min(Re(eigen(H_inv)$values)) < 0) {
    H_inv_pd <- nearPD(H_inv)
    R2_inv <- chol(H_inv_pd$mat)
    print(paste0("H_inv: absolute eigenvalue difference to nearest p.d. matrix: ", 
                 sum(abs(eigen(H_inv)$values - eigen(H_inv_pd$mat)$values))))
  } else {
    R2_inv <- chol(H_inv)
  }
  # Obtain the variance adjustment matrix
  R2 <- solve(R2_inv)
  R2R1 <- R2 %*% R1
  
  # Apply variance adjustment to parameters
  par_adj <- apply(unc_par_samps, 1, DEadj, par_hat = unc_par_hat, R2R1 = R2R1, 
                   simplify = FALSE)
  par_adj <- matrix(unlist(par_adj), byrow=TRUE, nrow=M)
  
  #=============== Convert adjusted to constrained space =======================
  # Constrained adjusted parameters for all MCMC samples
  pi_red_adj <- matrix(NA, nrow=M, ncol=K)
  theta_red_adj <- array(NA, dim=c(M, p, K, d))
  xi_red_adj <- array(NA, dim=c(M, K, q))
  for (i in 1:M) {
    constr_pars <- constrain_pars(out_stan, par_adj[i,])
    pi_red_adj[i, ] <- constr_pars$pi
    theta_red_adj[i,,,] <- constr_pars$theta
    xi_red_adj[i,,] <- constr_pars$xi
  }
  
  #=============== Output adjusted parameters ==================================
  pi_med_adj <- apply(pi_red_adj, 2, median)
  theta_med_adj <- apply(theta_red_adj, c(2,3,4), median)
  xi_med_adj <- apply(xi_red_adj, c(2,3), median)
  
  # Update Phi_med using adjusted xi estimate
  c_all <- analysis$c_all
  Phi_med_all_c <- pnorm(V %*% t(xi_med_adj))  # Outcome probabilities for all classes
  Phi_med_adj <- numeric(n)                    # Initialize individual outcome probabilities
  # Calculate posterior class membership, p(c_i=k|-), for each class k
  for (i in 1:n) {
    # Calculate outcome probabilities P(Y=1|-) using updated class assignment
    Phi_med_adj[i] <- Phi_med_all_c[i, c_all[i]] 
  }
  
  analysis_adj <- list(pi_red_adj = pi_red_adj, theta_red_adj = theta_red_adj, 
                       xi_red_adj = xi_red_adj, pi_med_adj = pi_med_adj, 
                       theta_med_adj = theta_med_adj, xi_med_adj = xi_med_adj, 
                       Phi_med_adj = Phi_med_adj, c_all = c_all,
                       pred_class_probs = analysis$pred_class_probs,
                       log_lik_med = analysis$loglik_med)
  return(analysis_adj)
}

set.seed(1)
analysis_adj1 <- var_adjust(mod_stan = mod_stan, analysis = analysis, 
                           K = analysis$K_red, p = p, d = d, n = n, q = q, 
                           x_mat = x_mat, y_all = y_all, V = V, w_all = w_all, 
                           s_all = s_all, clus_id_all = clus_id_all)

#===== R
source(paste0(wd, model_dir, "helper_functions.R"))

set.seed(1)
clus_id_all <- data_vars$cluster_id  # Cluster indicators, nx1
analysis_adj2 <- var_adjust(mod_stan = mod_stan, analysis = analysis, 
                           K = analysis$K_red, p = p, d = d, n = n, q = q, 
                           x_mat = x_mat, y_all = y_all, V = V, w_all = w_all, 
                           s_all = s_all, clus_id_all = clus_id_all)

identical(analysis_adj1, analysis_adj2)

#================ FAILED but very similar
#================ After adding in cluster id to R code, PASSED!


#================ Try whole function
#===== R
set.seed(3)
results_adj1 <- WSOLCA_main(data_path = data_path, res_path = res_path,
                           adj_path = adj_path, stan_path = stan_path, 
                           save_res = FALSE, n_runs = 2000, burn = 1000, 
                           thin = 5)

#===== Rcpp
# Source R helper functions
source(paste0(wd, model_dir, "helper_functions.R"))
# Source Rcpp functions
Rcpp::sourceCpp(paste0(wd, model_dir, "main_Rcpp_functions.cpp"))
set.seed(3)
results_adj2 <- WSOLCA_main_Rcpp(data_path = data_path, adapt_path = adapt_path,
                                 adj_path = adj_path, stan_path = stan_path, 
                                 save_res = FALSE, n_runs = 2000, 
                                 burn = 1000, thin = thin, covs = "true_Si")

