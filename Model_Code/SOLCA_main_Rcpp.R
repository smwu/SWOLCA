#===================================================
## Supervised Overfitted Latent Class Model
## Programmer: SM Wu   
## Data: Simulations       
#===================================================

# Read in command line arguments
args <- commandArgs(trailingOnly = TRUE)
scen_samp <- args[[1]]  # Simulation scenario
iter_pop <- args[[2]]   # Population iteration
samp_n <- args[[3]]     # Sample number

# Load libraries
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
library(Rcpp)
library(RcppArmadillo)
library(RcppTN)

#========================= HELPER FUNCTIONs ====================================

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

# `run_MCMC_Rcpp` runs the Gibbs sampler MCMC algorithm to obtain posterior samples
# Inputs:
#   OLCA_params: output from 'init_OLCA' containing 'pi', 'c_all', 'theta'
#   probit_params: output from 'init_probit' containing 'xi', 'z_all'
#   n_runs: Number of MCMC iterations
#   burn: Burn-in period
#   thin: Thinning factor
#   K: Number of classes
#   p: Number of exposure items
#   d: Number of exposure categories
#   n: Number of individuals
#   q: Number of regression covariates excluding class assignment
#   w_all: Weights normalized to sum to n. nx1
#   x_mat: Categorical exposure matrix. nxp
#   y_all: Vector of outcomes. nx1
#   V: Regression design matrix without class assignment. nxq
#   alpha: Vector of hyperparameters for pi. Kx1
#   eta: Vector of hyperparameters for theta. dx1
#   mu0: List of vectors of mean hyperparameters for xi. List of K qx1 vectors
#   Sig0: List of matrices of variance hyperparameters for xi. List of K qxq matrices
# Outputs: `MCMC_out` list containing the following items:
#   pi_MCMC: Matrix of posterior samples for pi. (n_iter)xK
#   theta_MCMC: Array of posterior samples for theta. (n_iter)xpxKxd
#   xi_MCMC: Array of posterior samples for xi. (n_iter)xKxq
#   c_all_MCMC: Matrix of posterior samples for c_all. (n_iter)xn
#   z_all_MCMC: Matrix of posterior samples for z_all. (n_iter)xn
#   loglik_MCMC: Vector of posterior samples for log-likelihood. (n_iter)x1
run_MCMC_Rcpp <- function(OLCA_params, probit_params, n_runs, burn, thin, K, p, d, n, 
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
    pi <- update_pi(pi = pi, w_all = w_all, c_all = c_all, K = K, alpha = alpha)
    # print(pi[1:10])
    c_all <- update_c(c_all = c_all, n = n, K = K, p = p, theta = theta, 
                      x_mat = x_mat, pi = pi, z_all = z_all, V = V, xi = xi, 
                      y_all = y_all)
    # print(c_all[1:10])
    theta <-update_theta(theta = theta, p = p, K = K, d = d, eta = eta, 
                         w_all = w_all, c_all = c_all, x_mat = x_mat)
    # print(theta[1:5,1:10,1])
    xi <- update_xi(xi = xi, n = n, K = K, w_all = w_all, c_all = c_all, 
                    z_all = z_all, V = V, y_all = y_all, mu0 = mu0, Sig0 = Sig0)
    # print(xi[1:3,])
    z_all <- update_z(z_all = z_all, n = n, V = V, xi = xi, c_all = c_all, 
                      y_all = y_all)
    # print(z_all[1:10])
    loglik <- update_loglik(loglik = loglik, n = n, p = p, c_all = c_all, 
                            theta = theta, x_mat = x_mat, pi = pi, 
                            z_all = z_all, V = V, xi = xi, y_all = y_all)
    # print(loglik[1:10])
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

# `get_mode` obtains the modal value given an input vector
# Inputs:
#   v: Input vector
# Outputs most common value found in input vector `v`
get_mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# 'post_process' conducts post-processing to cluster individuals into a 
# reduced number of classes and reorder posterior parameter samples according to
# the reduced number of classes
# Inputs:
#   MCMC_out: Output from 'run_MCMC' containing 'pi_MCMC', 'theta_MCMC', 'xi_MCMC', 'c_all_MCMC', 'z_all_MCMC', and 'loglik_MCMC'
#   p: Number of exposure items
#   d: Number of exposure categories
#   q: Number of regression covariates excluding class assignment
# Outputs: `post_MCMC_out` list containing the following items:
#   K_med: Median, across iterations, of number of classes with >= 5% of individuals
#   pi: Matrix of reduced and relabeled posterior samples for pi. (n_iter)x(K_med)
#   theta: Array of reduced and relabeled posterior samples for theta. (n_iter)xpx(K_med)xd
#   xi: Array of reduced and relabeled posterior samples for xi. (n_iter)x(K_med)xq
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

# `analyze_results` obtains posterior parameter samples and estimates
# Inputs:
#   MCMC_out: Output from 'run_MCMC' containing 'pi_MCMC', 'theta_MCMC', 'xi_MCMC', 
# 'c_all_MCMC', 'z_all_MCMC', 'loglik_MCMC'
#   post_MCMC_out: output from 'post_process' containing 'K_med', 'pi', 'theta', 'xi'
#   n: Number of individuals
#   p: Number of exposure items
#   V: Regression design matrix without class assignment. nxq
#   y_all: Vector of outcomes. nx1
#   x_mat: Categorical exposure matrix. nxp
# Outputs: `analysis` list containing the following items:
#   K_red: Number of unique classes
#   pi_red: Matrix of final posterior samples for pi. Mx(K_red)
#   theta_red: Array of final posterior samples for theta. Mxpx(K_red)xd
#   xi_red: Array of final posterior samples for xi. Mx(K_red)xq
#   pi_med: Vector of posterior median estimates for pi. (K_red)x1
#   theta_med: Array of posterior median estimates for theta. px(K_red)xd
#   xi_med: Matrix of posterior median estimates for xi. (K_red)xq
#   Phi_med: Vector of final individual outcome probabilities. nx1
#   c_all: Vector of final individual class assignments. nx1
#   pred_class_probs: Matrix of individual posterior class probabilities. nx(K_red)
#   loglik_med: Vector of final indiviudal log-likehoods. nx1
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


#========================= MAIN FUNCTION =======================================

# 'SOLCA_main_Rcpp' runs the SOLCA model and saves and returns results
# Inputs:
#   data_path: String path for input dataset
#   res_path: String path for output file
#   save_res: Boolean specifying if results should be saved. Default = TRUE
#   n_runs: Number of MCMC iterations
#   burn: Burn-in period
#   thin: Thinning factor
# Outputs: Saves and returns list `res` containing:
#   analysis: List of posterior model results
#   runtime: Total runtime for model
#   data_vars: Input dataset
#   MCMC_out: List of MCMC output
# Also saved 'analysis' MCMC output prior to variance adjustment
SOLCA_main_Rcpp <- function(data_path, res_path, save_res = TRUE, n_runs, 
                            burn, thin) {
  start_time <- Sys.time()
  
  #================= Read in data ==============================================
  print("Read in data")
  # CHANGE TO RDATA
  data_vars <- readMat(data_path)$sim.data
  names(data_vars) <- str_replace_all(dimnames(data_vars)[[1]], "[.]", "_")
  
  # Obtain dimensions
  n <- dim(data_vars$X_data)[1]        # Number of individuals
  p <- dim(data_vars$X_data)[2]        # Number of exposure items
  d <- max(apply(data_vars$X_data, 2,  # Number of exposure categories 
                 function(x) length(unique(x))))  # CHANGE TO ADAPT TO ITEM
  # Obtain data
  x_mat <- data_vars$X_data            # Categorical exposure matrix, nxp
  y_all <- c(data_vars$Y_data)         # Binary outcome vector, nx1
  s_all <- data_vars$true_Si           # Stratifying variable, nx1
  S <- length(unique(s_all))           # Number of strata
  s_mat <- dummy_cols(data.frame(s = factor(s_all)),  # Stratifying variable as dummy columns
                      remove_selected_columns = TRUE)
  q <- S                               # Number of regression covariates excluding class assignment
  V <- as.matrix(s_mat)                # Regression design matrix without class assignment, nxq
  
  # Set normalized weights to 1
  w_all <- rep(1, n)
  
  #================= ADAPTIVE SAMPLER ==========================================
  print("Adaptive sampler")
  #================= Initialize priors and variables for OLCA model ============
  K_max <- 30                      # Upper limit for number of classes
  alpha <- rep(1, K_max) / K_max   # Hyperparameter for prior for pi
  eta <- rep(1, d)                 # Hyperparameter for prior for theta
  # Obtain pi, theta, c_all
  OLCA_params <- init_OLCA(alpha = alpha, eta = eta, n = n, K = K_max, p = p, 
                           d = d)
  
  #================= Initialize priors and variables for probit model ==========
  # Initialize hyperparameters for xi
  mu0 <- Sig0 <- vector("list", K_max)
  for (k in 1:K_max) {
    # MVN(0,1) hyperprior for prior mean of xi
    mu0[[k]] <- rnorm(n = q)
    # InvGamma(3.5, 6.25) hyperprior for prior variance of xi. Assume uncorrelated 
    # components and mean variance 2.5 for a weakly informative prior on xi
    Sig0[[k]] <- diag(rinvgamma(n = q, shape = 3.5, scale = 6.25), nrow = q, ncol = q)
  }
  # Obtain xi, z_all                        
  probit_params <- init_probit(mu0 = mu0, Sig0 = Sig0, K = K_max, q = q, n = n, 
                               V = V, y_all = y_all, c_all = OLCA_params$c_all)
  
  #================= Run adaptive sampler to obtain number of classes ==========
  # Obtain pi_MCMC, theta_MCMC, xi_MCMC, c_all_MCMC, z_all_MCMC, loglik_MCMC
  MCMC_out <- run_MCMC_Rcpp(OLCA_params = OLCA_params, probit_params = probit_params, 
                            n_runs = n_runs, burn = burn, thin = thin, K = K_max, 
                            p = p, d = d, n = n, q = q, w_all = w_all, x_mat = x_mat, 
                            y_all = y_all, V = V, alpha = alpha, eta = eta, 
                            Sig0 = Sig0, mu0 = mu0)
  
  #================= Post-processing for adaptive sampler ======================
  # Get median number of classes with >= 5% of individuals, over all iterations
  M <- dim(MCMC_out$pi_MCMC)[1]  # Number of stored MCMC iterations
  K_med <- round(median(rowSums(MCMC_out$pi_MCMC >= 0.05)))
  # Get number of unique classes for fixed sampler
  K_fixed <- K_med
  print(paste0("K_fixed: ", K_fixed))
  # Reduce memory burden
  rm(OLCA_params, probit_params, MCMC_out)
  
  # # Obtain K_med, pi, theta, xi
  # post_MCMC_out <- post_process(MCMC_out = MCMC_out, p = p, d = d, q = q)
  # 
  # # Identify unique classes using modal exposure categories
  # # Posterior median estimate for theta across iterations
  # theta_med_temp <- apply(post_MCMC_out$theta, c(2, 3, 4), median)
  # # Posterior modal exposure categories for each exposure item and reduced class
  # theta_modes <- apply(theta_med_temp, c(1, 2), which.max)
  # # Identify unique classes
  # unique_classes <- which(!duplicated(theta_modes, MARGIN = 2))
  # 
  # # Get number of unique classes for fixed sampler
  # K_fixed <- length(unique_classes) 
  # 
  # # Reduce memory burden
  # rm(OLCA_params, probit_params, MCMC_out, post_MCMC_out)
  
  #================= FIXED SAMPLER =============================================
  print("Fixed sampler")
  #================= Run fixed sampler to obtain posteriors ====================
  # Initialize OLCA model using fixed number of classes
  alpha <- rep(1, K_fixed) / K_fixed  # Hyperparameter for prior for pi
  # Obtain pi, theta, c_all
  OLCA_params <- init_OLCA(alpha = alpha, eta = eta, n = n, K = K_fixed, p = p, 
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
  probit_params <- init_probit(mu0 = mu0, Sig0 = Sig0, K = K_fixed, q = q, n = n, 
                               V = V, y_all = y_all, c_all = OLCA_params$c_all)
  
  # Run MCMC algorithm using fixed number of classes
  # Obtain pi_MCMC, theta_MCMC, xi_MCMC, c_all_MCMC, z_all_MCMC, loglik_MCMC
  MCMC_out <- run_MCMC_Rcpp(OLCA_params = OLCA_params, probit_params = probit_params, 
                            n_runs = n_runs, burn = burn, thin = thin, K = K_fixed, 
                            p = p, d = d, n = n, q = q, w_all = w_all, x_mat = x_mat, 
                            y_all = y_all, V = V, alpha = alpha, eta = eta, 
                            Sig0 = Sig0, mu0 = mu0)
  
  # Post-processing to recalibrate labels and remove extraneous empty classes
  # Obtain K_med, pi, theta, xi, loglik_MCMC
  post_MCMC_out <- post_process(MCMC_out = MCMC_out, p = p, d = d, q = q)
  
  # Obtain posterior estimates, reduce number of classes, analyze results
  # Obtain K_red, pi_red, theta_red, xi_red, pi_med, theta_med, xi_med, Phi_med, 
  # c_all, pred_class_probs, loglik_med
  analysis <- analyze_results(MCMC_out = MCMC_out, post_MCMC_out = post_MCMC_out, 
                              n = n, p = p, V = V, y_all = y_all, x_mat = x_mat)
  
  runtime <- Sys.time() - start_time
  
  #================= Save and return output ====================================
  res <- list(analysis = analysis, runtime = runtime, 
              data_vars = data_vars, MCMC_out = MCMC_out)
  if (save_res) {
    save(res, file = res_path)
  }
  return(res)
}

#===================== Running main SOLCA function ============================

# Define directories
wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/"
# wd <- "~/Documents/Harvard/Research/Briana/supRPC/wsOFMM/"
data_dir <- "Data/"
res_dir <- "Results/"
model_dir <- "Model_Code/"
model <- "sOFMM"

# # Testing code
# scen_samp <- 201
# iter_pop <- 1
# samp_n <- 1

# Define paths
# REMOVE ITER_POP
data_path <- paste0(wd, data_dir, "simdata_scen", scen_samp, "_iter", iter_pop,
                    "_samp", samp_n, ".mat")   # Input dataset
res_path <- paste0(wd, res_dir, model, "_results_scen", scen_samp, 
                   "_samp", samp_n, ".RData")  # Output file

# Check if results already exist
already_done <- file.exists(res_path)
if (already_done) {
  print(paste0('Scenario ', scen_samp, ' iter ', iter_pop, ' samp ', samp_n,
               ' already exists.'))
} else {
  Rcpp::sourceCpp(paste0(wd, model_dir, "WSOLCA_main_Rcpp.cpp"))
  set.seed(samp_n)
  print(paste0("Running WSOLCA_main for scenario ", scen_samp, ' iter ', 
               iter_pop,' samp ', samp_n))
  results <- SOLCA_main_Rcpp(data_path = data_path, res_path = res_path,
                                  save_res = TRUE, n_runs = 25000, burn = 15000, 
                                  thin = 5)
  print(paste0("Runtime: ", results$runtime))
}

