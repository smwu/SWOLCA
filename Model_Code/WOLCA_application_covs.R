#===================================================
## Weighted Overfitted Latent Class Model
## Programmer: SM Wu   
## Data: NHANES Application with Covariates
## Date Updated: 2023/07/15
#===================================================

# Load libraries
library(plyr)
library(dplyr)
library(tidyr)
library(forcats)
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
library(ggpubr)

#========================= MAIN FUNCTION =======================================

# 'WOLCA_app_covs_Rcpp' runs the WOLCA model with covariates and saves and 
# returns results
# Inputs:
#   data_vars: Output list from "process_data" function 
#   adapt_path: String path for adaptive sampler file
#   res_path: String path for output file
#   save_res: Boolean specifying if results should be saved. Default = TRUE
#   n_runs: Number of MCMC iterations
#   burn: Burn-in period
#   thin: Thinning factor
#   K_known: Number of latent classes. If NULL (default), adaptive sampler runs
# Outputs: Saves and returns list `res` containing:
#   analysis: List of posterior model results
#   runtime: Total runtime for model
#   data_vars: Input dataset
#   MCMC_out: List of full MCMC output
#   post_MCMC_out: List of MCMC output after relabeling
#   K_MCMC: Adaptive sampler MCMC output for K
# Also saves list `adapt_MCMC` containing:
#   MCMC_out: List of full MCMC output
#   K_fixed: Number of classes to use for fixed sampler; output from adaptive sampler
#   K_MCMC: Adaptive sampler MCMC output for K
WOLCA_app_covs_Rcpp <- function(data_vars, adapt_path, res_path, save_res = TRUE, 
                            n_runs, burn, thin, K_known = NULL) {
  start_time <- Sys.time()
  
  #================= Read in data ==============================================
  print("Read in data")
  x_mat <- data_vars$x_mat
  y_all <- data_vars$y_all
  s_all <- data_vars$s_all
  clus_id_all <- data_vars$clus_id_all
  sample_wt <- data_vars$sample_wt
  V <- data_vars$V
  
  # Obtain dimensions
  n <- dim(x_mat)[1]        # Number of individuals
  p <- dim(x_mat)[2]        # Number of exposure items
  d <- max(apply(x_mat, 2,  # Number of exposure categories 
                 function(x) length(unique(x))))  # CHANGE TO ADAPT TO ITEM
  q <- ncol(V)                         # Number of regression covariates excluding class assignment
  
  # Obtain normalized weights
  kappa <- sum(sample_wt) / n   # Weights norm. constant. If sum(weights)=N, this is 1/(sampl_frac)
  w_all <- c(sample_wt / kappa) # Weights normalized to sum to n, nx1
  
  #================= ADAPTIVE SAMPLER ==========================================
  if (!is.null(K_known)) {
    K_fixed <- K_known
  } else {
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
                                    n_runs = n_runs, burn = burn, thin = thin, 
                                    K = K_max, p = p, d = d, n = n, w_all = w_all, 
                                    x_mat = x_mat, alpha = alpha, eta = eta)
    
    #================= Post-processing for adaptive sampler ======================
    # Get median number of classes with >= 5% of individuals, over all iterations
    M <- dim(MCMC_out$pi_MCMC)[1]  # Number of stored MCMC iterations
    K_MCMC <- rowSums(MCMC_out$pi_MCMC >= 0.05)
    K_med <- round(median(K_MCMC))
    # Get number of unique classes for fixed sampler
    K_fixed <- K_med
    print(paste0("K_fixed: ", K_fixed))
    # Save adaptive output
    adapt_MCMC <- list(MCMC_out = MCMC_out, K_fixed = K_fixed, K_MCMC = K_MCMC)
    if (save_res) {
      save(adapt_MCMC, file = adapt_path)
    }
    # Reduce memory burden
    rm(OLCA_params, MCMC_out)
  }  
  
  #================= FIXED SAMPLER =============================================
  print("Fixed sampler")
  #================= Run fixed sampler to obtain posteriors ====================
  # Initialize OLCA model using fixed number of classes
  alpha <- rep(1, K_fixed) / K_fixed  # Hyperparameter for prior for pi
  eta <- rep(1, d)                 # Hyperparameter for prior for theta
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
  svy_data <- data.frame(s = factor(s_all),
                         clus = clus_id_all,
                         x = x_mat,
                         y = y_all, 
                         wts = w_all,
                         c = factor(analysis$c_all),
                         age_cat = data_vars$V_data$age_cat,
                         racethnic = data_vars$V_data$racethnic,
                         smoker = data_vars$V_data$smoker,
                         physactive = data_vars$V_data$physactive)
  svydes <- svydesign(ids = ~clus, strata = ~s, weights = ~wts, data = svy_data)
  fit <- svyglm(y ~ c + age_cat + racethnic + smoker + physactive + 
                  c:age_cat + c:racethnic + c:smoker + c:physactive, 
                design = svydes, family = quasibinomial(link = "probit")) 
  coefs <- fit$coefficients
  ci <- confint(fit)
  
  # If zero/negative residual df, manually calculate the Wald confidence interval 
  # using a t-distribution with degrees of freedom from the survey design. 
  # Best if no cluster-level covariates in the regression model
  if (all(is.na(ci))) {
    ci <- manual_CI(model_object = fit, svy_df = degf(svydes), ci = 0.95)[, -1]
  }

  # Convert format to match WSOLCA and SOLCA
  xi_med <- xi_med_lb <- xi_med_ub <- matrix(NA, nrow = analysis$K_red, ncol = q)
  xi_med[1, ] <- coefs[c(1, (analysis$K_red + 1:(q-1)))]
  xi_med_lb[1, ] <- ci[c(1, (analysis$K_red + 1:(q-1))), 1]
  xi_med_ub[1, ] <- ci[c(1, (analysis$K_red + 1:(q-1))), 2]
  for (k in 2:analysis$K_red) {
    xi_med[k, ] <- coefs[c(k, (k + (q-1)) + (analysis$K_red-1) * (1:(q-1)))] + 
      xi_med[1, ]
    xi_med_lb[k, ] <- ci[c(k, (k + (q-1)) + (analysis$K_red-1) * (1:(q-1))), 1] + 
      xi_med_lb[1, ]
    xi_med_ub[k, ] <- ci[c(k, (k + (q-1)) + (analysis$K_red-1) * (1:(q-1))), 2] + 
      xi_med_ub[1, ]
  }
  
  analysis$xi_med <- xi_med
  analysis$xi_med_lb <- xi_med_lb
  analysis$xi_med_ub <- xi_med_ub
  analysis$fit <- fit
  
  runtime <- Sys.time() - start_time
  
  #================= Save and return output ====================================
  res <- list(analysis = analysis, runtime = runtime, 
              data_vars = data_vars, V = V, MCMC_out = MCMC_out, 
              post_MCMC_out = post_MCMC_out)
  if (is.null(K_known)) {
    res$K_MCMC <- K_MCMC
  }
  if (save_res) {
    save(res, file = res_path)
  }
  return(res)
}



#===================== RUN MAIN WOLCA FUNCTION =================================

# Define directories
wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/SWOLCA/"
# wd <- "~/Documents/Harvard/Research/Briana/supRPC/SWOLCA/"
data_dir <- "Data/"
res_dir <- "Results/July6/"
model_dir <- "Model_Code/"
model <- "wOFMM"

# Define paths
data_path <- paste0(wd, data_dir, "nhanes1518_adult_low_f_12jul2023.csv")   # Input dataset
adapt_path <- paste0(wd, res_dir, model, "_adapt_nhanesNOEDUC_NOLEG", ".RData")  # Adaptive output file
res_path <- paste0(wd, res_dir, model, "_results_adj_nhanesNOEDUC_NOLEG", ".RData")  # Output file

# Check if results already exist
already_done <- file.exists(res_path)
if (already_done) {
  print(paste0('NHANES results already exist.'))
} else {
  # Source application preparation helper functions
  source(paste0(wd, model_dir, "app_helper_functions.R"))
  # Read and process data
  data_vars <- process_data(data_path = data_path,
                            covs = c("age_cat", "racethnic", "smoker", "physactive"),
                            formula = "~ age_cat + racethnic + smoker + physactive")
  
  # Source R helper functions
  source(paste0(wd, model_dir, "helper_functions.R"))
  # Source Rcpp functions
  Rcpp::sourceCpp(paste0(wd, model_dir, "main_Rcpp_functions.cpp"))

  # Set seed and run model
  set.seed(20230225)
  print(paste0("Running WOLCA_application..."))
  results <- WOLCA_app_covs_Rcpp(data_vars = data_vars, 
                                 adapt_path = adapt_path,
                                 res_path = res_path, save_res = TRUE, 
                                 n_runs = 20000, burn = 10000, thin = 5, 
                                 K_known = NULL)
  print(paste0("Runtime: ", results$runtime))
}



#===================== PLOT OUTPUT =============================================

load(res_path)
age_categs <- c("[20,40)", "[40,60)", ">=60")
educ_categs <- c("Some College", "HS/GED", "<HS")
racethnic_categs <- c("NH White", "NH Black", "NH Asian", "Hispanic/Latino", "Other/Mixed")
smoker_categs <- c("Non-Smoker", "Smoker")
physactive_categs <- c("Inactive", "Active")

# Reorder classes
new_order <- c(3, 2, 5, 4, 1)
# new_order <- c(5, 2, 4, 1, 3)
res <- reorder_classes(res = res, model = "wOFMM", new_order = new_order)

plot_theta_modes(res, model = "wOFMM")
plot_theta_probs(res, model = "wOFMM")
plot_Phi_line(res, model = "wOFMM", age_categs = age_categs, 
              racethnic_categs = racethnic_categs,
              educ_categs = educ_categs, smoker_categs = smoker_categs, 
              physactive_categs = physactive_categs, ymax = 0.9)
plot_pi_boxplots(res, model = "wOFMM")

# Output reference cell coefficients table for xi 
convert_to_ref_wolca(xi_med = res$analysis$xi_med, 
                     xi_med_lb = res$analysis$xi_med_lb,
                     xi_med_ub = res$analysis$xi_med_ub,
                     age_categs = age_categs, racethnic_categs = racethnic_categs, 
                     smoker_categs = smoker_categs, 
                     physactive_categs = physactive_categs, format = "latex")


