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

#========================= MAIN FUNCTION =======================================

# 'WOLCA_main_Rcpp' runs the WOLCA model and saves and returns results
# Inputs:
#   data_path: String path for input dataset
#   adapt_path: String path for adaptive sampler file
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
WOLCA_main_Rcpp <- function(data_path, adapt_path, res_path, save_res = TRUE, 
                            n_runs, burn, thin, covs = NULL) {
  start_time <- Sys.time()
  
  #================= Read in data ==============================================
  print("Read in data")
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
  if (!is.null(covs)) {  # Other covariates are included in the probit model  
    s_all <- data_vars[[covs]]    # Stratifying variable, nx1
    # Stratifying variable as dummy columns
    s_mat <- dummy_cols(data.frame(s = factor(s_all)),  
                        remove_selected_columns = TRUE)
    V <- as.matrix(s_mat)              # Regression design matrix without class assignment, nxq
    q <- ncol(V)                       # Number of regression covariates excluding class assignment
  } else {  # Only latent class is included in the probit model
    s_all <- NULL  # No stratifying variable
    V <- matrix(1, nrow = n) 
    q <- 1
  }
  
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
                                  n_runs = round(n_runs/2), burn = round(burn/2), 
                                  thin = thin, K = K_max, p = p, d = d, n = n, 
                                  w_all = w_all, x_mat = x_mat, alpha = alpha, 
                                  eta = eta)
  
  #================= Post-processing for adaptive sampler ======================
  # Get median number of classes with >= 5% of individuals, over all iterations
  M <- dim(MCMC_out$pi_MCMC)[1]  # Number of stored MCMC iterations
  K_med <- round(median(rowSums(MCMC_out$pi_MCMC >= 0.05)))
  # Get number of unique classes for fixed sampler
  K_fixed <- K_med
  print(paste0("K_fixed: ", K_fixed))
  # Save adaptive output
  adapt_MCMC <- list(MCMC_out = MCMC_out, K_fixed = K_fixed)
  if (save_res) {
    save(adapt_MCMC, file = adapt_path)
  }
  # Reduce memory burden
  rm(OLCA_params, MCMC_out)
  
  
  #================= FIXED SAMPLER =============================================
  print("Fixed sampler")
  #================= Run fixed sampler to obtain posteriors ====================
  # Initialize OLCA model using fixed number of classes
  alpha <- rep(2, K_fixed) # Hyperparameter for prior for pi
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
  if (!is.null(s_all)) {  # Include stratifying variable 
    # V_ref <- data.frame(s = factor(s_all), c = factor(analysis$c_all), y = y_all)
    # fit <- glm(y ~ s * c, data = V_ref, family = binomial(link = "probit"))
    svy_data <- data.frame(x = x_mat,
                           y = y_all, 
                           s = factor(s_all),
                           c = factor(analysis$c_all),
                           wts = w_all)
    svydes <- svydesign(ids = ~1, strata = ~s, weights = ~wts, data = svy_data)
    fit <- svyglm(y ~ s * c, design = svydes, family = quasibinomial(link = "probit"))
  } else {  # No stratifying variable 
    # V_ref <- data.frame(c = factor(analysis$c_all), y = y_all)
    # fit <- glm(y ~ c, data = V_ref, family = binomial(link = "probit"))
    svy_data <- data.frame(x = x_mat,
                           y = y_all, 
                           c = factor(analysis$c_all),
                           wts = w_all)
    svydes <- svydesign(ids = ~1, weights = ~wts, data = svy_data)
    fit <- svyglm(y ~ c, design = svydes, family = quasibinomial(link = "probit"))
  }
  coefs <- fit$coefficients
  ci <- confint(fit)
  
  # Convert format to match WSOLCA and SOLCA
  xi_med <- xi_med_lb <- xi_med_ub <- matrix(NA, nrow = analysis$K_red, ncol = q)
  if (!is.null(s_all)) {  # Include stratifying variable 
    for (k in 1:analysis$K_red) {
      for (s in 1:q) {
        xi_med[k, s] <- coefs[1] + (k != 1) * coefs[q + (k-1)] + (s != 1) * coefs[s] + 
          (k != 1) * (s != 1) * coefs[q + (analysis$K_red-1) + (k-1)]
        xi_med_lb[k, s] <- ci[1, 1] + (k != 1) * ci[q + (k-1), 1] + 
          (s != 1) * ci[s, 1] + 
          (k != 1) * (s != 1) * ci[q + (analysis$K_red-1) + (k-1), 1]
        xi_med_ub[k, s] <- ci[1, 2] + (k != 1) * ci[q + (k-1), 2] + 
          (s != 1) * ci[s, 2] + 
          (k != 1) * (s != 1) * ci[q + (analysis$K_red-1) + (k-1), 2]
      }
    }
  } else {  # Stratifying variable not included in probit model
    for (k in 1:analysis$K_red) {
      xi_med[k, 1] <- coefs[1] + (k != 1) * coefs[q + (k-1)] 
      xi_med_lb[k, 1] <- ci[1, 1] + (k != 1) * ci[q + (k-1), 1] 
      xi_med_ub[k, 1] <- ci[1, 2] + (k != 1) * ci[q + (k-1), 2]
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



#===================== RUN MAIN WOLCA FUNCTION =================================

# Define directories
wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/"
# wd <- "~/Documents/Harvard/Research/Briana/supRPC/wsOFMM/"
data_dir <- "Data/"
res_dir <- "Results/"
model_dir <- "Model_Code/"
model <- "wOFMM"

# # Testing code
# scen_samp <- 111111
# iter_pop <- 1
# samp_n <- 1
# n_runs <- 100
# burn <- 50
# thin <- 5
# save_res <- FALSE
# covs <- "true_Si"

# Define paths
data_path <- paste0(wd, data_dir, "simdata_scen", scen_samp, "_iter", iter_pop,
                    "_samp", samp_n, ".RData")   # Input dataset
adapt_path <- paste0(wd, res_dir, model, "_adapt_scen", scen_samp, 
                     "_samp", samp_n, ".RData")  # Output file
res_path <- paste0(wd, res_dir, model, "_results_wt_scen", scen_samp, 
                   "_samp", samp_n, ".RData")  # Output file

# Check if results already exist
already_done <- file.exists(res_path)
if (already_done) {
  print(paste0('Scenario ', scen_samp, ' iter ', iter_pop, ' samp ', samp_n,
               ' already exists.'))
} else {
  n_runs <- 20000
  burn <- 10000
  thin <- 5
  save_res <- TRUE
  covs <- "true_Si"

  # Source R helper functions
  source(paste0(wd, model_dir, "helper_functions.R"))
  # Source Rcpp functions
  Rcpp::sourceCpp(paste0(wd, model_dir, "main_Rcpp_functions.cpp"))
  # Set seed
  set.seed(samp_n)
  # Run model
  print(paste0("Running WOLCA_main for scenario ", scen_samp, ' iter ', 
               iter_pop,' samp ', samp_n))
  results <- WOLCA_main_Rcpp(data_path = data_path, adapt_path = adapt_path,
                             res_path = res_path,
                             save_res = save_res, n_runs = n_runs, 
                             burn = burn, thin = thin, covs = covs)
  print(paste0("Runtime: ", results$runtime))
}

