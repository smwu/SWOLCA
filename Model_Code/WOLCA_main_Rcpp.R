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

# `run_MCMC_Rcpp_WOLCA` runs the Gibbs sampler MCMC algorithm to obtain posterior 
# samples for the two-step unsupervised WOLCA model
# Inputs:
#   OLCA_params: output from 'init_OLCA' containing 'pi', 'c_all', 'theta'
#   n_runs: Number of MCMC iterations
#   burn: Burn-in period
#   thin: Thinning factor
#   K: Number of classes
#   p: Number of exposure items
#   d: Number of exposure categories
#   n: Number of individuals
#   w_all: Weights normalized to sum to n. nx1
#   x_mat: Categorical exposure matrix. nxp
#   alpha: Vector of hyperparameters for pi. Kx1
#   eta: Vector of hyperparameters for theta. dx1
# Outputs: `MCMC_out` list containing the following items:
#   pi_MCMC: Matrix of posterior samples for pi. (n_iter)xK
#   theta_MCMC: Array of posterior samples for theta. (n_iter)xpxKxd
#   c_all_MCMC: Matrix of posterior samples for c_all. (n_iter)xn
run_MCMC_Rcpp_WOLCA <- function(OLCA_params, n_runs, burn, thin, K, p, d, n, 
                               w_all, x_mat, alpha, eta) {
  n_storage <- ceiling(n_runs / thin)  # Number of MCMC iterations to store
  
  # Initialize variables
  pi_MCMC <- matrix(NA, nrow = n_storage, ncol = K)
  theta_MCMC <- array(NA, dim = c(n_storage, p, K, d))
  c_all_MCMC <- z_all_MCMC <- matrix(NA, nrow = n_storage, ncol = n)

  # Initialized values
  pi <- OLCA_params$pi
  theta <- OLCA_params$theta
  c_all <- OLCA_params$c_all

  alpha_post <- numeric(K)                            # Posterior parameters for pi
  eta_post <- numeric(d)                              # Posterior parameters for theta

  # Update parameters and variables
  for (m in 1:n_runs) {
    pi <- update_pi(pi = pi, w_all = w_all, c_all = c_all, K = K, alpha = alpha)
    # print(pi[1:10])
    c_all <- update_c_WOLCA(c_all = c_all, n = n, K = K, p = p, theta = theta, 
                           x_mat = x_mat, pi = pi)
    # print(c_all[1:10])
    theta <-update_theta(theta = theta, p = p, K = K, d = d, eta = eta, 
                         w_all = w_all, c_all = c_all, x_mat = x_mat)
    # print(theta[1:5,1:10,1])

    #============== Store posterior values based on thinning  ==================
    if (m %% thin == 0) {
      m_thin <- m / thin
      pi_MCMC[m_thin, ] <- pi
      theta_MCMC[m_thin, , , ] <- theta
      c_all_MCMC[m_thin, ] <- c_all
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
      print(paste0("Iteration ", m, " completed!"))
    }
  }
  
  # Discard burn-in
  warmup <- ceiling(burn / thin)
  pi_MCMC <- pi_MCMC[-(1:warmup), ]
  theta_MCMC <- theta_MCMC[-(1:warmup), , , ]
  c_all_MCMC <- c_all_MCMC[-(1:warmup), ]

  MCMC_out <- list(pi_MCMC = pi_MCMC, theta_MCMC = theta_MCMC, 
                   c_all_MCMC = c_all_MCMC)
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

# 'post_process_WOLCA' conducts post-processing to cluster individuals into a 
# reduced number of classes and reorder posterior parameter samples according to
# the reduced number of classes
# Inputs:
#   MCMC_out: Output from 'run_MCMC' containing 'pi_MCMC', 'theta_MCMC', 'xi_MCMC', 'c_all_MCMC', 'z_all_MCMC', and 'loglik_MCMC'
#   p: Number of exposure items
#   d: Number of exposure categories
# Outputs: `post_MCMC_out` list containing the following items:
#   K_med: Median, across iterations, of number of classes with >= 5% of individuals
#   pi: Matrix of reduced and relabeled posterior samples for pi. (n_iter)x(K_med)
#   theta: Array of reduced and relabeled posterior samples for theta. (n_iter)xpx(K_med)xd
post_process_WOLCA <- function(MCMC_out, p, d) {
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
  for (m in 1:M) {
    iter_order <- relabel_red_classes[m, ]
    pi_temp <- MCMC_out$pi_MCMC[m, iter_order]
    pi[m, ] <- pi_temp / sum(pi_temp)
    theta[m, , , ] <- MCMC_out$theta_MCMC[m, , iter_order, ]
  }
  
  post_MCMC_out <- list(K_med = K_med, pi = pi, theta = theta)
  
  return(post_MCMC_out)  
  # plot(dendrogram, labels = FALSE)
  # rect.hclust(dendrogram, k = K_med)
}

# `analyze_results_WOLCA` obtains posterior parameter samples and estimates
# Inputs:
#   MCMC_out: Output from 'run_MCMC_Rcpp_WOLCA' containing 'pi_MCMC', 
# 'theta_MCMC', 'c_all_MCMC'
#   post_MCMC_out: output from 'post_process_WOLCA' containing 'K_med', 'pi', 'theta'
#   n: Number of individuals
#   p: Number of exposure items
#   x_mat: Categorical exposure matrix. nxp
# Outputs: `analysis` list containing the following items:
#   K_red: Number of unique classes
#   pi_red: Matrix of final posterior samples for pi. Mx(K_red)
#   theta_red: Array of final posterior samples for theta. Mxpx(K_red)xd
#   pi_med: Vector of posterior median estimates for pi. (K_red)x1
#   theta_med: Array of posterior median estimates for theta. px(K_red)xd
#   c_all: Vector of final individual class assignments. nx1
#   pred_class_probs: Matrix of individual posterior class probabilities. nx(K_red)
analyze_results_WOLCA <- function(MCMC_out, post_MCMC_out, n, p, x_mat) {
  
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

  #============== Posterior median estimates ===================================
  pi_med <- apply(pi_red, 2, median, na.rm = TRUE)
  pi_med <- pi_med / sum(pi_med)  # Re-normalize
  theta_med <- apply(theta_red, c(2, 3, 4), median, na.rm = TRUE)
  theta_med <- aaply(theta_med, c(1, 2), function(x) x / sum(x))  # Re-normalize

  #============== Update c using unique classes and posterior estimates ========
  c_all <- MCMC_out$c_all_MCMC[M, ]                      # Final class assignments
  pred_class_probs <- matrix(NA, nrow = n, ncol = K_red) # Posterior class membership probabilities
  log_cond_c <- matrix(NA, nrow = n, ncol = K_red)       # Individual log-likelihood for each class
  # Calculate posterior class membership, p(c_i=k|-), for each class k
  for (i in 1:n) {
    for (k in 1:K_red) {
      # Calculate theta component of individual log-likelihood assuming class k
      log_theta_comp_k <- 0
      for (j in 1:p) {
        log_theta_comp_k <- log_theta_comp_k + log(theta_med[j, k, x_mat[i, j]])
      }
      # Individual log-likelihood for class k
      log_cond_c[i, k] <- log(pi_med[k]) + log_theta_comp_k
    }
    # Calculate p(c_i=k|-) = p(x,y,c_i=k) / p(x,y)
    pred_class_probs[i, ] <- exp(log_cond_c[i, ] - logSumExp(log_cond_c[i, ]))
    # Update class assignment using the posterior probabilities
    c_all[i] <- rcat(n = 1, p = pred_class_probs[i, ])
  }
  
  analysis <- list(K_red = K_red, pi_red = pi_red, theta_red = theta_red, 
                   pi_med = pi_med, theta_med = theta_med, c_all = c_all, 
                   pred_class_probs = pred_class_probs)
  
  return(analysis)
}


#========================= MAIN FUNCTION =======================================

# 'WOLCA_main_Rcpp' runs the WOLCA model and saves and returns results
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
WOLCA_main_Rcpp <- function(data_path, res_path, save_res = TRUE, n_runs, 
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
  
  # Obtain normalized weights
  kappa <- sum(data_vars$sample_wt) / n   # Weights norm. constant. If sum(weights)=N, this is 1/(sampl_frac)
  w_all <- c(data_vars$sample_wt / kappa) # Weights normalized to sum to n, nx1
  
  #================= ADAPTIVE SAMPLER ==========================================
  print("Adaptive sampler")
  #================= Initialize priors and variables for OLCA model ============
  K_max <- 30                      # Upper limit for number of classes
  alpha <- rep(1, K_max) / K_max   # Hyperparameter for prior for pi
  eta <- rep(1, d)                 # Hyperparameter for prior for theta
  # Obtain pi, theta, c_all
  OLCA_params <- init_OLCA(alpha = alpha, eta = eta, n = n, K = K_max, p = p, 
                           d = d)
  
  #================= Run adaptive sampler to obtain number of classes ==========
  # Obtain pi_MCMC, theta_MCMC, c_all_MCMC
  MCMC_out <- run_MCMC_Rcpp_WOLCA(OLCA_params = OLCA_params,  
                            n_runs = n_runs, burn = burn, thin = thin, K = K_max, 
                            p = p, d = d, n = n, w_all = w_all, x_mat = x_mat, 
                            alpha = alpha, eta = eta)
  
  #================= Post-processing for adaptive sampler ======================
  # Get median number of classes with >= 5% of individuals, over all iterations
  M <- dim(MCMC_out$pi_MCMC)[1]  # Number of stored MCMC iterations
  K_med <- round(median(rowSums(MCMC_out$pi_MCMC >= 0.05)))
  # Get number of unique classes for fixed sampler
  K_fixed <- K_med
  print(paste0("K_fixed: ", K_fixed))
  # Reduce memory burden
  rm(OLCA_params, MCMC_out)
  
  # # Obtain K_med, pi, theta, xi
  # post_MCMC_out <- post_process_WOLCA(MCMC_out = MCMC_out, p = p, d = d)
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
  
  # Run MCMC algorithm using fixed number of classes
  # Obtain pi_MCMC, theta_MCMC, c_all_MCMC
  MCMC_out <- run_MCMC_Rcpp_WOLCA(OLCA_params = OLCA_params,  
                                 n_runs = n_runs, burn = burn, thin = thin, K = K_fixed, 
                                 p = p, d = d, n = n, w_all = w_all, x_mat = x_mat, 
                                 alpha = alpha, eta = eta)
  
  # Post-processing to recalibrate labels and remove extraneous empty classes
  # Obtain K_med, pi, theta, xi, loglik_MCMC
  post_MCMC_out <- post_process_WOLCA(MCMC_out = MCMC_out, p = p, d = d)
  
  # Obtain posterior estimates, reduce number of classes, analyze results
  # Obtain K_red, pi_red, theta_red, xi_red, pi_med, theta_med, xi_med, Phi_med, 
  # c_all, pred_class_probs, loglik_med
  analysis <- analyze_results_WOLCA(MCMC_out = MCMC_out, post_MCMC_out = post_MCMC_out, 
                                  n = n, p = p, x_mat = x_mat)
  
  #================= Fit probit model ==========================================
  V_ref <- data.frame(s = factor(s_all), c = factor(analysis$c_all), y = y_all)
  fit <- glm(y ~ s * c, data = V_ref, family = binomial(link = "probit"))
  coefs <- fit$coefficients
  ci <- confint(fit)
  
  # Convert format to match WSOLCA and SOLCA
  xi_med <- xi_med_lb <- xi_med_ub <- matrix(NA, nrow = analysis$K_red, ncol = S)
  for (k in 1:analysis$K_red) {
    for (s in 1:S) {
      xi_med[k, s] <- coefs[1] + (k != 1) * coefs[S + (k-1)] + (s != 1) * coefs[s] + 
        (k != 1) * (s != 1) * coefs[S + (analysis$K_red-1) + (k-1)]
      xi_med_lb[k, s] <- ci[1, 1] + (k != 1) * ci[S + (k-1), 1] + 
        (s != 1) * ci[s, 1] + 
        (k != 1) * (s != 1) * ci[S + (analysis$K_red-1) + (k-1), 1]
      xi_med_ub[k, s] <- ci[1, 2] + (k != 1) * ci[S + (k-1), 2] + 
        (s != 1) * ci[s, 2] + 
        (k != 1) * (s != 1) * ci[S + (analysis$K_red-1) + (k-1), 2]
    }
  }
  
  analysis$xi_med <- xi_med
  analysis$xi_med_lb <- xi_med_lb
  analysis$xi_med_ub <- xi_med_ub
  
  runtime <- Sys.time() - start_time
  
  #================= Save and return output ====================================
  res <- list(analysis = analysis, runtime = runtime, 
              data_vars = data_vars, MCMC_out = MCMC_out)
  if (save_res) {
    save(res, file = res_path)
  }
  return(res)
}

#===================== Running main WOLCA function ============================

# Define directories
wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/"
# wd <- "~/Documents/Harvard/Research/Briana/supRPC/wsOFMM/"
data_dir <- "Data/"
res_dir <- "Results/"
model_dir <- "Model_Code/"
model <- "wOFMM"

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
  print(paste0("Running WOLCA_main for scenario ", scen_samp, ' iter ', 
               iter_pop,' samp ', samp_n))
  results <- WOLCA_main_Rcpp(data_path = data_path, res_path = res_path,
                             save_res = TRUE, n_runs = 25000, burn = 15000, 
                             thin = 5)
  print(paste0("Runtime: ", results$runtime))
}

