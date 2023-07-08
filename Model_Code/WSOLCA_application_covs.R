#===================================================
## Weighted Supervised Overfitted Latent Class Model
## Programmer: SM Wu   
## Data: NHANES Application with Covariates   
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

#========================= Helper functions ====================================
# `process_data` reads in data and converts it into a processed format for the 
# 'WSOLCA_app_covs_Rcpp' function
# Input:
#   data_path: String specifying path for input data
# Output: list 'data_vars' containing the following objects:
#   x_mat
#   y_all
#   s_all
#   clus_id_all
#   sample_wt
#   V
process_data <- function(data_path, covs = "all") {
  # Read in data
  data_vars_all <- read.csv(data_path)
  # Filter for low-income women and complete cases: n = 3028
  data_vars_f_low <- data_vars_all %>% 
    filter(sex == "female", poverty == "At or Below") %>%
    select(SEQN, RIDAGEYR, racethnic, EDU, SDMVPSU, SDMVSTRA, obese, smoker,
           bp.flag, sumrisk, dietwt8yr, bb1:bb29) %>% na.omit()
  
  # Obtain exposure and outcome data
  x_mat <- as.matrix(data_vars_f_low %>% select(bb1:bb29))
  y_all <- data_vars_f_low$bp.flag
  
  # Get stratum IDs
  s_all <- data_vars_f_low$SDMVSTRA  # 59 strata
  # Create unique PSU IDs using stratum IDs
  max_psus <- max(data_vars_f_low$SDMVPSU)
  clus_id_all <- (s_all - 1) * max_psus + data_vars_f_low$SDMVPSU

  # Other covariates to be included in the probit model 
  ##### CHANGE THIS!!!!!!
  if (is.null(covs)) {
    V <- matrix(1, nrow = n) 
    q <- 1
  } else {
    V_data <- data_vars_f_low %>% 
      select(RIDAGEYR, racethnic, EDU, smoker) %>%
      mutate(racethnic = factor(racethnic),
             educ = factor(EDU),
             smoker = factor(smoker), .keep = "unused") %>%
      mutate(
        racethnic = case_when(
          racethnic == "NH White" ~ 1,
          racethnic == "NH Black" ~ 2,
          racethnic == "NH Asian" ~ 3,
          racethnic == "Mexican" | racethnic == "Other Hispanic" ~ 4,
          racethnic == "Other/Mixed" ~ 5),
        educ = case_when(
          educ == "At least some college" ~ 1,
          educ == "HS/GED" ~ 2,
          educ == "less than HS" ~ 3),
        age_cat = case_when(  # See breakdown of 20-44, 45-64, >=65. If too skewed, just use <=45
          20 <= RIDAGEYR & RIDAGEYR <= 29 ~ 1,
          30 <= RIDAGEYR & RIDAGEYR <= 39 ~ 2,
          40 <= RIDAGEYR & RIDAGEYR <= 49 ~ 3,
          50 <= RIDAGEYR & RIDAGEYR <= 59 ~ 4,
          60 <= RIDAGEYR & RIDAGEYR <= 69 ~ 5,
          RIDAGEYR >= 70 ~ 6), .keep = "unused")
    
    # Regression design matrix without class assignment, nxq
    # Exclude stratifying variable as well
    V <- model.matrix(~ age_cat + racethnic + educ + smoker, V_data)
    # V <- as.matrix(other_covs %>% 
    #                  dummy_cols(select_columns = c("age_cat", "racethnic", "educ", "smoker"), 
    #                             remove_selected_columns = TRUE))
    # V <- V[, -c(1, 7, 12, 15)]
  }
  
  # Get sampling weights
  sample_wt <- data_vars_f_low$dietwt8yr
  
  # Return processed data
  data_vars <- list(x_mat = x_mat, y_all = y_all, s_all = s_all, 
                         clus_id_all = clus_id_all, sample_wt = sample_wt, 
                         V = V)
  return(data_vars)
}


#========================= MAIN FUNCTION =======================================

# 'WSOLCA_app_covs_Rcpp' runs the WSOLCA model with covariates and saves and 
# returns results
# Inputs:
#   data_vars: Output list from "process_data" function 
#   res_path: String path for output file
#   adj_path: String path for adjusted output file
#   stan_path: String path for Stan file
#   save_res: Boolean specifying if results should be saved. Default = TRUE
#   n_runs: Number of MCMC iterations
#   burn: Burn-in period
#   thin: Thinning factor
# Outputs: Saves and returns list `res` containing:
#   analysis_adj: List of adjusted posterior model results
#   runtime: Total runtime for model
#   data_vars: Input dataset
# Also saved 'analysis' MCMC output prior to variance adjustment
WSOLCA_app_covs_Rcpp <- function(data_vars, adapt_path, adj_path, stan_path, 
                            save_res = TRUE, n_runs, burn, thin) {
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
  rm(OLCA_params, probit_params, MCMC_out)

  
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
  
  #================= VARIANCE ADJUSTMENT =======================================
  print("Variance adjustment")
  # Create Stan model
  mod_stan <- stan_model(stan_path)
  # Apply variance adjustment for correct coverage
  # Obtain pi_red_adj, theta_red_adj, xi_red_adj, pi_med_adj, theta_med_adj, 
  # xi_med_adj, Phi_med_adj, c_all, pred_class_probs, log_lik_med
  analysis_adj <- var_adjust(mod_stan = mod_stan, analysis = analysis, 
                             K = analysis$K_red, p = p, d = d, n = n, q = q, 
                             x_mat = x_mat, y_all = y_all, V = V, w_all = w_all, 
                             s_all = s_all, clus_id_all = clus_id_all)
  
  runtime <- Sys.time() - start_time
  
  #================= Save and return output ====================================
  res <- list(analysis_adj = analysis_adj, runtime = runtime, 
              data_vars = data_vars, MCMC_out = MCMC_out, 
              post_MCMC_out = post_MCMC_out,
              K_MCMC = adapt_MCMC$K_MCMC)
  if (save_res) {
    save(res, file = adj_path)
  }
  return(res)
}



#===================== RUN MAIN WSOLCA FUNCTION =================================

# Define directories
wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/"
# wd <- "~/Documents/Github/wsOFMM/"
data_dir <- "Data/"
res_dir <- "Results/June22/"
model_dir <- "Model_Code/"
model <- "wsOFMM"

# Define paths
data_path <- paste0(wd, data_dir, "nhanesallage_frscores1118.csv")   # Input dataset
adapt_path <- paste0(wd, res_dir, model, "_adapt_nhanes", ".RData")  # Output file
adj_path <- paste0(wd, res_dir, model, "_results_adj_nhanes", ".RData")  # Adjusted output file
stan_path <- paste0(wd, model_dir, "WSOLCA_main.stan")  # Stan file

# Check if results already exist
already_done <- file.exists(adj_path)
if (already_done) {
  print(paste0('NHANES results already exist.'))
} else {
  # Source R helper functions
  source(paste0(wd, model_dir, "helper_functions.R"))
  # Source Rcpp functions
  Rcpp::sourceCpp(paste0(wd, model_dir, "main_Rcpp_functions.cpp"))
  # Read and process data
  data_vars <- process_data(data_path = data_path)
  # Set seed
  set.seed(20230225)
  print(paste0("Running WSOLCA_application..."))
  results_adj <- WSOLCA_app_covs_Rcpp(data_vars = data_vars, res_path = res_path,
                                 adj_path = adj_path, stan_path = stan_path, 
                                 save_res = TRUE, n_runs = 20000, burn = 10000, 
                                 thin = 5)
  print(paste0("Runtime: ", results_adj$runtime))
}


#===================== PLOT OUTPUT =============================================
load(adj_path)
age_categs <- c("[20,30)", "[30,40)", "[40,50)", "[50,60)", "[60,70)", ">=70")
educ_categs <- c("Some College", "HS/GED", "<HS")
racethnic_categs <- c("NH White", "NH Black", "NH Asian", "Hispanic", "Other/Mixed")
smoker_categs <- c("No", "Yes")
K <- length(res$analysis_adj$pi_med_adj)

#### Plot pi boxplots for each latent class
pi_red <- as.data.frame(res$analysis_adj$pi_red_adj)
colnames(pi_red) <- paste0("pi_", 1:(dim(pi_red)[2]))
pi_red_plot <- pi_red %>% pivot_longer(cols = everything(), names_to = "pi_comp", 
                                       values_to = "value")
pi_red_plot %>% ggplot(aes(x = pi_comp, y = value, fill = pi_comp)) + 
  theme_bw() + scale_fill_brewer(palette="Set2") + 
  geom_boxplot() + 
  ggtitle(as.expression(bquote("Parameter estimation for "~pi~" "))) +
  xlab("Latent Class") + ylab(as.expression(bquote(~pi~"Value"))) + ylim(0,0.5)

#### Plot xi boxplots, separately for each covariate, then combined 
xi_dims <- dim(res$analysis_adj$xi_red_adj)
K <- xi_dims[2]
xi_red <- as.data.frame(t(matrix(res$analysis_adj$xi_red_adj, 
                                 nrow = xi_dims[1], 
                                 ncol = K*xi_dims[3])))
xi_red$Class <- as.character(factor(c(rep(1:K, times = xi_dims[3]))))
xi_red$Covariate_Level <- rep(c(age_categs, educ_categs, racethnic_categs, 
                                smoker_categs), each = K)
xi_red$Covariate <- c(rep("Age", K*6), rep("Education", K*3), 
                      rep("Race/Ethnicity", K*5), rep("Smoking", K*2))
xi_red_plot <- xi_red %>% pivot_longer(cols = -c("Class", "Covariate", "Covariate_Level"), 
                                       names_to = "iter", values_to = "value")
# xi_red_plot %>% ggplot(aes(x = Covariate_Level, y = value, group = Class, fill = Class)) + 
#   theme_bw() + 
#   geom_boxplot() + 
#   facet_grid(.~Covariate, labeller = label_both) +
#   ggtitle(as.expression(bquote("Parameter estimation for "~xi~" ")))
# xi_red_plot %>% ggplot(aes(x = Stratum, y = value, group = Stratum, 
#                            fill = Stratum)) + 
#   theme_bw() + 
#   geom_boxplot() +
#   facet_grid(.~Class, labeller = label_both) + 
#   ggtitle(as.expression(bquote("Parameter estimation for "~xi~" ")))
p3_age <- xi_red_plot %>% filter(Covariate == "Age") %>%
  ggplot(aes(x = Covariate_Level, y = value, fill = Class)) +
  theme_bw() + scale_fill_brewer(palette="Set2") + 
  geom_boxplot() + xlab("Age")
p3_race <- xi_red_plot %>% filter(Covariate == "Race/Ethnicity") %>%
  ggplot(aes(x = Covariate_Level, y = value, fill = Class)) +
  theme_bw() + scale_fill_brewer(palette="Set2") + 
  geom_boxplot() + xlab("Race/Ethnicity")
p3_educ <- xi_red_plot %>% filter(Covariate == "Education") %>%
  ggplot(aes(x = Covariate_Level, y = value, fill = Class)) +
  theme_bw() + scale_fill_brewer(palette="Set2") + 
  geom_boxplot() + xlab("Education")
p3_smoke <- xi_red_plot %>% filter(Covariate == "Smoking") %>%
  ggplot(aes(x = Covariate_Level, y = value, fill = Class)) +
  theme_bw() + scale_fill_brewer(palette="Set2") + 
  geom_boxplot() + xlab("Smoking")
ggarrange(p3_age, p3_race, p3_educ, p3_smoke, common.legend = TRUE,
          legend = "right", ncol = 2, nrow = 2)

#### Plot Phi boxplots separatey for each covariate
Phi_red_plot <- xi_red_plot %>%
  mutate(value = pnorm(value))
Phi_red_plot %>% filter(Covariate == "Age") %>%
  ggplot(aes(x = Covariate_Level, y = value, fill = Class)) +
  theme_bw() + scale_fill_brewer(palette="Set2") + 
  geom_boxplot() + ylab("Hypertension Risk") + xlab("Age")
Phi_red_plot %>% filter(Covariate == "Education") %>%
  ggplot(aes(x = Covariate_Level, y = value, fill = Class)) +
  theme_bw() + scale_fill_brewer(palette="Set2") + 
  geom_boxplot() + ylab("Hypertension Risk") + xlab("Education")
Phi_red_plot %>% filter(Covariate == "Race/Ethnicity") %>%
  ggplot(aes(x = Covariate_Level, y = value, fill = Class)) +
  theme_bw() + scale_fill_brewer(palette="Set2") + 
  geom_boxplot() + ylab("Hypertension Risk") + xlab("Race/Ethnicity")
Phi_red_plot %>% filter(Covariate == "Smoking") %>%
  ggplot(aes(x = Covariate_Level, y = value, fill = Class)) +
  theme_bw() + scale_fill_brewer(palette="Set2") + 
  geom_boxplot() + ylab("Hypertension Risk") + xlab("Smoking Status")


#### Plot Phi line plots, separately for each covariate, then combined 
df_Phi <- data.frame(Stratum = factor(res$data_vars$s_all), 
                     Class = factor(res$analysis_adj$c_all), 
                     Phi = res$analysis_adj$Phi_med)
# df_Phi %>% ggplot(aes(x = Class, y = Phi, col = Stratum, group = Stratum)) + 
#   theme_bw() + 
#   geom_point() + 
#   geom_line() + 
#   ylim(0, 0.5) + 
#   ggtitle("Parameter estimation for conditional outcome probabilities") +
#   ylab("P(Y=1|Class, Stratum)")
df_Phi <- data.frame(pnorm(res$analysis_adj$xi_med_adj))
colnames(df_Phi) <- c(age_categs, educ_categs, racethnic_categs, smoker_categs)
# colnames(df_Phi) <- colnames(res$data_vars$V)
df_Phi$Class <- factor(1:4)
df_Phi <- df_Phi %>% gather("Covariate_Level", "Phi", -Class)
df_Phi$Covariate <- c(rep("Age", K*6), rep("Education", K*3), 
                      rep("Race/Ethnicity", K*5), rep("Smoking", K*2))
# df_Phi %>% ggplot(aes(x = Covariate_Level, y = Phi, col = Class, group = Class)) +
#   theme_bw() + scale_color_brewer(palette="Set2") + 
#   geom_point() + geom_line() + 
#   facet_grid(~Covariate) +
#   ggtitle("Parameter estimation for conditional outcome probabilities") +
#   ylab("P(Y=1|-")
p2_age <- df_Phi %>% filter(Covariate == "Age") %>%
  ggplot(aes(x = Covariate_Level, y = Phi, col = Class, group = Class)) +
  theme_bw() + scale_color_brewer(palette="Set2") + 
  geom_point() + geom_line() + ylim(0,1) + 
  ylab("Hypertension Probability") + xlab("Age") +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10))
p2_educ <- df_Phi %>% filter(Covariate == "Education") %>%
  ggplot(aes(x = Covariate_Level, y = Phi, col = Class, group = Class)) +
  theme_bw() + scale_color_brewer(palette="Set2") + 
  geom_point() + geom_line() +  ylim(0,1) + 
  ylab("Hypertension Probability") + xlab("Education") +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10))
p2_race <- df_Phi %>% filter(Covariate == "Race/Ethnicity") %>%
  mutate(Covariate_Level = factor(Covariate_Level, levels = racethnic_categs)) %>%
  ggplot(aes(x = Covariate_Level, y = Phi, col = Class, group = Class)) +
  theme_bw() + scale_color_brewer(palette="Set2") + 
  geom_point() + geom_line() +  ylim(0,1) + 
  ylab("Hypertension Probability") + xlab("Race/Ethnicity") +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10)) +
  ggtitle("Hypertension Risk by Dietary Pattern and Race/Ethnicity")
p2_smoker <- df_Phi %>% filter(Covariate == "Smoking") %>%
  ggplot(aes(x = Covariate_Level, y = Phi, col = Class, group = Class)) +
  theme_bw() + scale_color_brewer(palette="Set2") + 
  geom_point() + geom_line() +  ylim(0,1) + 
  ylab("Hypertension Probability") + xlab("Smoking") +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10))
ggarrange(p_age, p_educ, p_race, p_smoker, ncol = 4,
          common.legend = TRUE, legend = "top", widths = c(1,0.7,1,0.5))
# df_Phi %>% ggplot(aes(x = Stratum, y = Phi, col = Class, group = Class)) + 
#   theme_bw() + 
#   geom_point() + 
#   geom_line() + 
#   ylim(0, 0.5) + 
#   ggtitle("Parameter estimation for conditional outcome probabilities") +
#   ylab("P(Y=1|Class, Stratum")

#### Plot Phi line plot OLD VERSION
df_Phi <- data.frame(Stratum = factor(res$data_vars$s_all), 
                     Class = factor(res$analysis_adj$c_all), 
                     Phi = res$analysis_adj$Phi_med)
df_Phi %>% ggplot(aes(x = Class, y = Phi, col = Stratum, group = Stratum)) + 
  theme_bw() + 
  geom_point() + 
  geom_line() + 
  ylim(0, 0.5) + 
  ggtitle("Parameter estimation for conditional outcome probabilities") +
  ylab("P(Y=1|Class, Stratum)")
df_Phi %>% ggplot(aes(x = Stratum, y = Phi, col = Class, group = Class)) + 
  theme_bw() + 
  geom_point() + 
  geom_line() + 
  ylim(0, 0.5) + 
  ggtitle("Parameter estimation for conditional outcome probabilities") +
  ylab("P(Y=1|Class, Stratum)")

#### Plot theta modes and probabilities
# Plot modal consumption levels
est_item_probs <- res$analysis_adj$theta_med_adj
mode_item_probs <- as.data.frame(apply(est_item_probs, c(1, 2), which.max))
food_items <- c("Citrus, melon, berries", "Other fruit", "Fruit juice", 
                "Dark green vegs", "Tomatoes", "Other red/orange vegs", 
                "Potatoes", "Other starchy vegs", "Other vegs", "Legumes (vegs)", 
                "Whole grains", "Refined grains", "Meat (unspec)", "Cured meats", 
                "Organ meat", "Poultry", "Seafood (high n3)", "Seafood (low n3)",
                "Eggs", "Soybean products", "Nuts and seeds", "Legumes (protein)", 
                "Milk", "Yogurt", "Cheese", "Oils", "Solid fats", "Added sugar", 
                "Alcoholic drinks")
class_names <- paste0("Class ", 1:(dim(est_item_probs)[2]))
rownames(mode_item_probs) <- food_items
colnames(mode_item_probs) <- class_names
mode_item_probs$Item <- rownames(mode_item_probs)
mode_plot <- mode_item_probs %>% gather("Class", "Level", -Item) 
mode_plot %>% ggplot(aes(x=Class, y=factor(Item, levels = rev(food_items)), 
                         fill=Level)) + 
  geom_tile(color="black") + 
  geom_text(aes(label = Level), col="white", cex=3) +
  scale_fill_gradient(trans="reverse") + 
  theme(legend.position="none") +
  ylab("Item") + xlab("Latent Class") +
  ggtitle("Estimated modal item \nconsumption levels") + 
  theme(text = element_text(size = 15))

# Plot probabilities of consumption levels
dimnames(est_item_probs)[[1]] <- food_items
dimnames(est_item_probs)[[2]] <- class_names
lcmodel <- reshape2::melt(est_item_probs, level=2)
colnames(lcmodel) <- c("Item", "Class", "Level", "Probability")
lcmodel %>%
  ggplot(aes(x = Probability, y = factor(Item, levels = rev(food_items)), 
             fill = factor(Level))) + 
  geom_bar(stat = "identity", position = "stack") + 
  facet_grid(. ~ Class) + 
  scale_fill_brewer(type="seq", palette="Greys") + 
  theme_bw() + 
  labs(x="Item consumption probabilities", y = "Food items",
       fill ="Item \nconsumption \nlevels") + 
  theme( axis.text.x=element_text(size=7),
         panel.grid.major.x=element_blank())
# axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1))


#=================== Plot latent classes by demographic ========================
res_demog <- other_covs %>%
  mutate(racethnic = factor(racethnic),
         smoker = factor(smoker),
         educ = factor(educ),
         age_cat = factor(age_cat))
res_demog$Class <- factor(res$analysis_adj$c_all)
res_demog$y <- res$data_vars$y_all
p_educ <- res_demog %>% ggplot(aes(x = educ, fill = Class, group = Class)) +
  theme_bw() + scale_fill_brewer(palette="Set2") + 
  geom_bar(position = "dodge") + 
  scale_x_discrete(labels = educ_categs) + 
  ggtitle("Dietary Patterns by Education") + xlab("Education") + ylab("Count")
p_age <- res_demog %>% ggplot(aes(x = age_cat, fill = Class, group = Class)) +
  theme_bw() + scale_fill_brewer(palette="Set2") + 
  geom_bar(position = "dodge") + 
  scale_x_discrete(labels = age_categs) + 
  ggtitle("Dietary Patterns by Age") + xlab("Age") + ylab("Count")
p_race <- res_demog %>% ggplot(aes(x = racethnic, fill = Class, group = Class)) +
  theme_bw() + scale_fill_brewer(palette="Set2") + 
  geom_bar(position = "dodge") + 
  scale_x_discrete(labels = racethnic_categs) + 
  ggtitle("Dietary Patterns by Race/Ethnicity") + xlab("Race/Ethnicity") + ylab("Count")
p_smoker <- res_demog %>% ggplot(aes(x = smoker, fill = Class, group = Class)) +
  theme_bw() + scale_fill_brewer(palette="Set2") + 
  geom_bar(position = "dodge") + 
  scale_x_discrete(labels = smoker_categs) + 
  ggtitle("Dietary Patterns by Smoking") + xlab("Smoking Status") + ylab("Count")
ggarrange(p_age, p_educ, p_race, p_smoker, ncol = 4, common.legend = TRUE,
          widths = c(1,0.7,1,0.7))
ggarrange(p_race, p2_race, ncol = 2, common.legend = TRUE, legend = "right")


#===================== Create colored dendrogram ===============================
# from https://atrebas.github.io/post/2019-06-08-lightweight-dendrograms/
dendro_data_k <- function(hc, k) {
  hcdata    <-  ggdendro::dendro_data(hc, type = "rectangle")
  seg       <-  hcdata$segments
  labclust  <-  cutree(hc, k)[hc$order]
  segclust  <-  rep(0L, nrow(seg))
  heights   <-  sort(hc$height, decreasing = TRUE)
  height    <-  mean(c(heights[k], heights[k - 1L]), na.rm = TRUE)
  for (i in 1:k) {
    xi      <-  hcdata$labels$x[labclust == i]
    idx1    <-  seg$x    >= min(xi) & seg$x    <= max(xi)
    idx2    <-  seg$xend >= min(xi) & seg$xend <= max(xi)
    idx3    <-  seg$yend < height
    idx     <-  idx1 & idx2 & idx3
    segclust[idx] <- i
  }
  idx                    <-  which(segclust == 0L)
  segclust[idx]          <-  segclust[idx + 1L]
  hcdata$segments$clust  <-  segclust
  hcdata$segments$line   <-  as.integer(segclust < 1L)
  hcdata$labels$clust    <-  labclust
  hcdata
}

set_labels_params <- function(nbLabels, direction = c("tb", "bt", "lr", "rl"),
                              fan = FALSE) {
  if (fan) {
    angle       <-  360 / nbLabels * 1:nbLabels + 90
    idx         <-  angle >= 90 & angle <= 270
    angle[idx]  <-  angle[idx] + 180
    hjust       <-  rep(0, nbLabels)
    hjust[idx]  <-  1
  } else {
    angle       <-  rep(0, nbLabels)
    hjust       <-  0
    if (direction %in% c("tb", "bt")) { angle <- angle + 45 }
    if (direction %in% c("tb", "rl")) { hjust <- 1 }
  }
  list(angle = angle, hjust = hjust, vjust = 0.5)
}

plot_ggdendro <- function(hcdata, direction = c("lr", "rl", "tb", "bt"),
                          fan = FALSE, scale.color = NULL, branch.size = 1,
                          label.size  = 3, nudge.label = 0.01, expand.y = 0.1) {
  direction <- match.arg(direction) # if fan = FALSE
  ybreaks   <- pretty(segment(hcdata)$y, n = 5)
  ymax      <- max(segment(hcdata)$y)
  ## branches
  p <- ggplot() +
    geom_segment(data         =  segment(hcdata),
                 aes(x        =  x,
                     y        =  y,
                     xend     =  xend,
                     yend     =  yend,
                     linetype =  factor(line),
                     colour   =  factor(clust)),
                 lineend      =  "round",
                 show.legend  =  FALSE,
                 size         =  branch.size)
  ## orientation
  if (fan) {
    p <- p +
      coord_polar(direction = -1) +
      scale_x_continuous(breaks = NULL,
                         limits = c(0, nrow(label(hcdata)))) +
      scale_y_reverse(breaks = ybreaks)
  } else {
    p <- p + scale_x_continuous(breaks = NULL)
    if (direction %in% c("rl", "lr")) {
      p <- p + coord_flip()
    }
    if (direction %in% c("bt", "lr")) {
      p <- p + scale_y_reverse(breaks = ybreaks)
    } else {
      p <- p + scale_y_continuous(breaks = ybreaks)
      nudge.label <- -(nudge.label)
    }
  }
  # labels
  labelParams <- set_labels_params(nrow(hcdata$labels), direction, fan)
  hcdata$labels$angle <- labelParams$angle
  p <- p +
    geom_text(data        =  label(hcdata),
              aes(x       =  x,
                  y       =  y,
                  label   =  label,
                  colour  =  factor(clust),
                  angle   =  angle),
              vjust       =  labelParams$vjust,
              hjust       =  labelParams$hjust,
              nudge_y     =  ymax * nudge.label,
              size        =  label.size,
              show.legend =  FALSE)
  # colors and limits
  if (!is.null(scale.color)) {
    p <- p + scale_color_manual(values = scale.color)
  }
  ylim <- -round(ymax * expand.y, 1)
  p    <- p + expand_limits(y = ylim)
  p
}

library(ggdendro)
MCMC_out <- res$MCMC_out
M <- dim(MCMC_out$pi_MCMC)[1] 
K_med <- round(median(rowSums(MCMC_out$pi_MCMC >= 0.05)))
distMat <- hamming.distance(t(MCMC_out$c_all_MCMC))
dendrogram <- hclust(as.dist(distMat), method = "complete")
dendro_k <- dendro_data_k(dendrogram, K_med)
plot_ggdendro(dendro_k, direction = "tb") + 
  theme_bw() + scale_color_brewer(palette="Set2") + 
  geom_hline(yintercept=1750, linetype = "dashed") +
  xlab("") + ylab("")



# Testing
# results_adj <- WSOLCA_app_Rcpp(data_path = data_path, res_path = res_path,
#                                 adj_path = adj_path, stan_path = stan_path,
#                                 save_res = FALSE, n_runs = 60, burn = 30,
#                                 thin = 3)
