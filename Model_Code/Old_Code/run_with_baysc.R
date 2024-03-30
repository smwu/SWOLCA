# Read in command line arguments
args <- commandArgs(trailingOnly = TRUE)
scen_samp <- args[[1]]  # Simulation scenario
iter_pop <- args[[2]]   # Population iteration
samp_n <- args[[3]]     # Sample number
covs <- args[[4]]       # Probit covariates
if (covs == 1) {
  covs <- NULL
} else if (covs == 2) {
  covs <- "true_Si"
} else if (covs == 3) {
  covs <- "additional"
}

# Load libraries
library(plyr)
library(dplyr)
library(LaplacesDemon)
library(truncnorm)
library(fastDummies)
library(matrixStats)
library(gtools)
library(e1071)
library(rstan)
library(survey)
library(Rcpp)
library(RcppArmadillo)
library(RcppTN)
library(stats)
library(Matrix)

#===================== RUN MAIN WSOLCA FUNCTION ================================

# Define directories
wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/SWOLCA/"
# wd <- "~/Documents/Github/SWOLCA/"
# data_dir <- "Data/"
# res_dir <- "Results/"
data_dir <- "Data/July6/"
res_dir <- "Results/July6/"
model_dir <- "Model_Code/"
model <- "wsOFMM"

# Testing code
scen_samp <- 111211
iter_pop <- 1
samp_n <- 100

n_runs <- 500
burn <- 250
thin <- 5
covs <- "true_Si"
save_res <- FALSE

# Define paths
data_path <- paste0(wd, data_dir, "simdata_scen", scen_samp, "_iter", iter_pop,
                    "_samp", samp_n, ".RData")   # Input dataset
adapt_path <- paste0(wd, res_dir, model, "_adapt_comb_scen", scen_samp,
                     "_samp", samp_n, ".RData")  # Output file
adj_path <- paste0(wd, res_dir, model, "_results_comb_adjRcpp_scen", scen_samp,
                   "_samp", samp_n, ".RData")      # Adjusted output file
stan_path <- paste0(wd, model_dir, "WSOLCA_main.stan")  # Stan file

# Check if results already exist
already_done <- file.exists(adj_path)
if (already_done) {
  print(paste0('Scenario ', scen_samp, ' iter ', iter_pop, ' samp ', samp_n,
               ' already exists.'))
} else {
  n_runs <- 20000
  burn <- 10000
  thin <- 5
  save_res <- TRUE
  # covs <- "true_Si"

  # Source R helper functions
  source(paste0(wd, model_dir, "helper_functions.R"))
  # Source Rcpp functions
  Rcpp::sourceCpp(paste0(wd, model_dir, "main_Rcpp_functions.cpp"))
  # Set seed
  set.seed(samp_n)

  adapt_seed <- samp_n
  fixed_seed <- 20230629

  # Run model
  print(paste0("Running WSOLCA_main for scenario ", scen_samp, ' iter ',
               iter_pop,' samp ', samp_n))
  results_adj <- WSOLCA_main_Rcpp(data_path = data_path, adapt_path = adapt_path,
                                  adj_path = adj_path, stan_path = stan_path,
                                  save_res = save_res, n_runs = n_runs,
                                  burn = burn, thin = thin, covs = covs)
  print(paste0("Runtime: ", results_adj$res$runtime))
}



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
covs <- args[[4]]       # Probit covariates
if (covs == 1) {
  covs <- NULL
} else if (covs == 2) {
  covs <- "true_Si"
} else if (covs == 3) {
  covs <- "additional"
}

# Load libraries
library(R.matlab)
library(stringr)
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

#===================== RUN MAIN SOLCA FUNCTION =================================

# Define directories
wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/SWOLCA/"
# wd <- "~/Documents/Harvard/Research/Briana/supRPC/SWOLCA/"
data_dir <- "Data/July6/"
res_dir <- "Results/July6/"
model_dir <- "Model_Code/"
model <- "sOFMM"

# # Testing code
# scen_samp <- 111111
# iter_pop <- 1
# samp_n <- 1
# n_runs <- 100
# burn <- 50
# thin <- 5
# covs <- "true_Si"
# save_res <- FALSE

# Define paths
data_path <- paste0(wd, data_dir, "simdata_scen", scen_samp, "_iter", iter_pop,
                    "_samp", samp_n, ".RData")   # Input dataset
# data_path <- paste0(wd, data_dir, "simdata_scen", scen_samp, "_iter", iter_pop,
#                     "_samp", samp_n, ".mat")   # Input dataset
adapt_path <- paste0(wd, res_dir, model, "_adapt_scen", scen_samp,
                     "_samp", samp_n, ".RData")  # Output file
res_path <- paste0(wd, res_dir, model, "_results_scen", scen_samp,
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
  # covs <- "true_Si"
  
  # Source R helper functions
  source(paste0(wd, model_dir, "helper_functions.R"))
  # Source Rcpp functions
  Rcpp::sourceCpp(paste0(wd, model_dir, "main_Rcpp_functions.cpp"))
  # Set seed
  set.seed(samp_n)
  # Run model
  print(paste0("Running SOLCA_main for scenario ", scen_samp, ' iter ',
               iter_pop,' samp ', samp_n))
  results <- SOLCA_main_Rcpp(data_path = data_path, adapt_path = adapt_path,
                             res_path = res_path,
                             save_res = save_res, n_runs = n_runs,
                             burn = burn, thin = thin, covs = covs)
  print(paste0("Runtime: ", results$res$runtime))
}


# Obtain data
x_mat <- data_vars$X_data            # Categorical exposure matrix, nxp
y_all <- c(data_vars$Y_data)         # Binary outcome vector, nx1
if (is.null(covs)) {
  # Probit model only includes latent class C
  # No stratifying variable
  s_all <- NULL
  V <- matrix(1, nrow = n)
  q <- 1
} else if (covs == "true_Si") {
  # Probit model includes C and S: C + S + C:S
  # Stratifying variable, nx1
  s_all <- data_vars[[covs]]
  # Regression design matrix without class assignment, nxq
  V_data <- data.frame(s = as.factor(s_all))
  V <- model.matrix(~ s, V_data)
  # Number of regression covariates excluding class assignment
  q <- ncol(V)
} else if (covs == "additional") {
  # Probit model includes C, S, A (binary), and B (continuous)
  # C + S + A + B + C:S + C:A + C:B
  # Stratifying variable, nx1
  s_all <- data_vars[["true_Si"]]
  a_all <- data_vars[["true_Ai"]]
  b_all <- data_vars[["true_Bi"]]
  # Regression design matrix without class assignment, nxq
  V_data <- data.frame(s = as.factor(s_all), a = as.factor(a_all), b = b_all)
  V <- model.matrix(~ s + a + b, V_data)
  # Number of regression covariates excluding class assignment
  q <- ncol(V)
} else {
  stop("Error: covs must be one of 'true_Si', 'additional', or NULL")
}



#===================================================
## Weighted Overfitted Latent Class Model
## Programmer: SM Wu
## Data: Simulations
#===================================================

# Read in command line arguments
args <- commandArgs(trailingOnly = TRUE)
scen_samp <- args[[1]]  # Simulation scenario
iter_pop <- args[[2]]   # Population iteration
samp_n <- args[[3]]     # Sample number
covs <- args[[4]]       # Probit covariates
if (covs == 1) {
  covs <- NULL
} else if (covs == 2) {
  covs <- "true_Si"
} else if (covs == 3) {
  covs <- "additional"
}

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

#===================== RUN MAIN WOLCA FUNCTION =================================

# Define directories
wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/SWOLCA/"
# wd <- "~/Documents/Harvard/Research/Briana/supRPC/SWOLCA/"
data_dir <- "Data/July6/"
res_dir <- "Results/July6/"
model_dir <- "Model_Code/"
model <- "wOFMM"

# # Testing code
# scen_samp <- 111211
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
  # covs <- "true_Si"
  
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

if (is.null(covs)) {
  # Probit model only includes latent class C
  svy_data <- data.frame(x = x_mat,
                         y = y_all,
                         wts = w_all,
                         c = factor(estimates$c_all),
                         clus = clus_id_all)
  svydes <- svydesign(ids = ~clus, weights = ~wts, data = svy_data)
  fit <- svyglm(y ~ c, design = svydes, family = quasibinomial(link = "probit"))
} else if (covs == "true_Si") {
  # Probit model includes C and S: C + S + C:S
  svy_data <- data.frame(s = factor(s_all),
                         x = x_mat,
                         y = y_all,
                         wts = w_all,
                         c = factor(estimates$c_all),
                         clus = clus_id_all)
  svydes <- svydesign(ids = ~clus, strata = ~s, weights = ~wts, data = svy_data)
  fit <- svyglm(y ~ c * s, design = svydes, family = quasibinomial(link = "probit"))
} else if (covs == "additional") {
  # Probit model includes C, S, A (binary), and B (continuous)
  # C + S + A + B + C:S + C:A + C:B
  svy_data <- data.frame(s = factor(s_all),
                         x = x_mat,
                         y = y_all,
                         wts = w_all,
                         c = factor(estimates$c_all),
                         a = factor(a_all),
                         b = b_all,
                         clus = clus_id_all)
  svydes <- svydesign(ids = ~clus, strata = ~s, weights = ~wts, data = svy_data)
  fit <- svyglm(y ~ c + s + a + b + c:s + c:a + c:b,
                design = svydes, family = quasibinomial(link = "probit"))
}
