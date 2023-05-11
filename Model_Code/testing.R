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

#================= Effect modifier scenario 11211 ==============================
# Define directories
wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/"
# wd <- "~/Documents/Harvard/Research/Briana/supRPC/wsOFMM/"
data_dir <- "Data/"
res_dir <- "Results/"
model_dir <- "Model_Code/"

# Source R helper functions
source(paste0(wd, model_dir, "helper_functions.R"))
# Source Rcpp functions
Rcpp::sourceCpp(paste0(wd, model_dir, "main_Rcpp_functions.cpp"))

# Testing code
scen_samp <- 11211
iter_pop <- 1
samp_n <- 1

#================= Run WSOLCA without S ========================================
model <- "wsOFMM"

# Define paths
data_path <- paste0(wd, data_dir, "simdata_scen", scen_samp, "_iter", iter_pop,
                    "_samp", samp_n, ".RData")   # Input dataset
res_path <- paste0(wd, res_dir, model, "_results_scen", scen_samp, 
                   "_samp", samp_n, ".RData")  # Output file
adj_path <- paste0(wd, res_dir, model, "_results_adjRcpp_scen", scen_samp, 
                   "_samp", samp_n, ".RData")      # Adjusted output file
stan_path <- paste0(wd, model_dir, "WSOLCA_main.stan")  # Stan file

# Set seed
set.seed(samp_n)
# Run model
print(paste0("Running WSOLCA_main for scenario ", scen_samp, ' iter ', 
             iter_pop,' samp ', samp_n))
results_adj <- WSOLCA_main_Rcpp(data_path = data_path, res_path = res_path,
                                adj_path = adj_path, stan_path = stan_path, 
                                save_res = TRUE, n_runs = 20000, burn = 10000, 
                                thin = 5, covs = NULL)
print(paste0("Runtime: ", results_adj$runtime))

#================= Run SOLCA without S =========================================
model <- "sOFMM"
# Define paths
data_path <- paste0(wd, data_dir, "simdata_scen", scen_samp, "_iter", iter_pop,
                    "_samp", samp_n, ".RData")   # Input dataset
res_path <- paste0(wd, res_dir, model, "_results_scen", scen_samp, 
                   "_samp", samp_n, ".RData")  # Output file
# Set seed
set.seed(samp_n)
# Run model
print(paste0("Running SOLCA_main for scenario ", scen_samp, ' iter ', 
             iter_pop,' samp ', samp_n))
results <- SOLCA_main_Rcpp(data_path = data_path, res_path = res_path,
                           save_res = TRUE, n_runs = 20000, burn = 10000, 
                           thin = 5, covs = NULL)
print(paste0("Runtime: ", results$runtime))

#================ Compare with true xi =========================================
load(data_path)
sim_data$true_Phi_mat
# Marginal P(Y=1|C)
prob_s <- sim_data$N_s / sim_data$N
prob_c <- sim_data$true_pi
prob_y_cond_c <- numeric(sim_data$true_K)
for (k in 1:sim_data$true_K) {
  prob_y_cond_c[k] <- sum(sim_data$true_Phi_mat[k, ] * 
                            sim_data$true_pi_s[, k] * prob_s) / prob_c[k]
}
prob_y_cond_c

load(res_path)
res_MCMC$post_MCMC_out$K_med
median(rowSums(res_MCMC$MCMC_out$pi_MCMC >= 0.05))
res_MCMC$MCMC_out$pi_MCMC[1,]


Rcpp::sourceCpp("Model_Code/test_cpp.cpp")
for (i in 1:10) {
  set.seed(1)
  alpha <- rep(1, 4) / 4
  print(rdirichlet_cpp(alpha))
}

