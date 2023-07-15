#===================================================
## Weighted Supervised Overfitted Latent Class Model
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

# 'WSOLCA_app_covs_Rcpp' runs the WSOLCA model with covariates and saves and 
# returns results
# Inputs:
#   data_vars: Output list from "process_data" function 
#   adapt_path: String path for adaptive sampler output file
#   adj_path: String path for adjusted output file
#   stan_path: String path for Stan file
#   save_res: Boolean specifying if results should be saved. Default = TRUE
#   n_runs: Number of MCMC iterations
#   burn: Burn-in period
#   thin: Thinning factor
#   K_known: Number of latent classes. If NULL (default), adaptive sampler runs
# Outputs: Saves and returns list `res` containing:
#   analysis_adj: List of posterior model results
#   runtime: Total runtime for model
#   data_vars: Input dataset
#   MCMC_out: List of full MCMC output
#   post_MCMC_out: List of MCMC output after relabeling
#   K_MCMC: Adaptive sampler MCMC output for K
# Also saves list `adapt_MCMC` containing:
#   MCMC_out: List of full MCMC output
#   K_fixed: Number of classes to use for fixed sampler; output from adaptive sampler
#   K_MCMC: Adaptive sampler MCMC output for K
WSOLCA_app_covs_Rcpp <- function(data_vars, adapt_path, adj_path, stan_path, 
                            save_res = TRUE, n_runs, burn, thin, K_known = NULL) {
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
              post_MCMC_out = post_MCMC_out)
  if (is.null(K_known)) {
    res$K_MCMC <- K_MCMC
  }
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
res_dir <- "Results/July6/"
model_dir <- "Model_Code/"
model <- "wsOFMM"

# Define paths
data_path <- paste0(wd, data_dir, "nhanes1518_adult_low_f_12jul2023.csv")   # Input dataset
adapt_path <- paste0(wd, res_dir, model, "_adapt_nhanesNOEDUC_NOLEG", ".RData")  # Output file
adj_path <- paste0(wd, res_dir, model, "_results_adj_nhanesNOEDUC_NOLEG", ".RData")  # Adjusted output file
stan_path <- paste0(wd, model_dir, "WSOLCA_main.stan")  # Stan file

# Check if results already exist
already_done <- file.exists(adj_path)
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
  print(paste0("Running WSOLCA_application..."))
  results_adj <- WSOLCA_app_covs_Rcpp(data_vars = data_vars, 
                                      adapt_path = adapt_path,
                                      adj_path = adj_path, 
                                      stan_path = stan_path, save_res = TRUE, 
                                      n_runs = 20000, burn = 10000, thin = 5,
                                      K_known = NULL)
  print(paste0("Runtime: ", results_adj$runtime))
}


#===================== PLOT OUTPUT =============================================

load(adj_path)
age_categs <- c("[20,40)", "[40,60)", ">=60")
educ_categs <- c("Some College", "HS/GED", "<HS")
racethnic_categs <- c("NH White", "NH Black", "NH Asian", "Hispanic/Latino", "Other/Mixed")
# racethnic_categs <- c("NH White", "NH Black", "NH Asian", "Mex-Amer", "Other Hispanic")
smoker_categs <- c("No", "Yes")
physactive_categs <- c("Inactive", "Active")
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

### Plot marginal model xi
#### Plot xi boxplots, separately for each covariate, then combined 
# ~ age_cat + racethnic + educ + smoker + Phys_Active
xi_dims <- dim(res$analysis_adj$xi_red_adj)
K <- xi_dims[2]
xi_red <- as.data.frame(t(matrix(res$analysis_adj$xi_red_adj, 
                                 nrow = xi_dims[1], 
                                 ncol = K*xi_dims[3], byrow = FALSE)))
# Class variable: 1-K, 1-K, 1-K...
xi_red$Class <- as.character(factor(c(rep(1:K, times = xi_dims[3]))))
# Covariate level: RefxK, [40,60)xK, >=60xK... 
xi_red$Covariate_Level <- rep("Ref", K)
xi_red$Covariate <- c(rep("Ref", K))
xi_red_plot <- xi_red %>% 
  pivot_longer(cols = -c("Class", "Covariate", "Covariate_Level"), 
               names_to = "iter", values_to = "value")
xi_red_plot %>% ggplot(aes(x = Covariate_Level, y = value, group = Class, fill = Class)) +
  theme_bw() +
  geom_boxplot() +
  facet_grid(.~Covariate, labeller = label_both) +
  ggtitle(as.expression(bquote("Parameter estimation for "~xi~" ")))


#### Plot xi boxplots, separately for each covariate, then combined 
# ~ age_cat + racethnic + educ + smoker + Phys_Active
xi_dims <- dim(res$analysis_adj$xi_red_adj)
K <- xi_dims[2]
xi_red <- as.data.frame(t(matrix(res$analysis_adj$xi_red_adj, 
                                 nrow = xi_dims[1], 
                                 ncol = K*xi_dims[3], byrow = FALSE)))
# Class variable: 1-K, 1-K, 1-K...
xi_red$Class <- as.character(factor(c(rep(1:K, times = xi_dims[3]))))
# Covariate level: RefxK, [40,60)xK, >=60xK... 
xi_red$Covariate_Level <- rep(c("Ref", age_categs[-1], racethnic_categs[-1], 
                                # educ_categs[-1], 
                                smoker_categs[-1],
                                physactive_categs[-1]), each = K)
xi_red$Covariate <- c(rep("Ref", K),
                      rep("Age", K*(length(age_categs) - 1)), 
                      rep("Race/Ethnicity", K*(length(racethnic_categs) - 1)), 
                      # rep("Education", K*(length(educ_categs) - 1)), 
                      rep("Smoking", K*(length(smoker_categs) - 1)),
                      rep("Physical Activity", K*(length(physactive_categs) - 1)))
  # xi_red$Covariate_Level <- rep(c(age_categs, educ_categs, racethnic_categs, 
  #                                 smoker_categs), each = K)
  # xi_red$Covariate <- c(rep("Age", K*6), rep("Education", K*3), 
  #                       rep("Race/Ethnicity", K*5), rep("Smoking", K*2))
xi_red_plot <- xi_red %>% 
  pivot_longer(cols = -c("Class", "Covariate", "Covariate_Level"), 
               names_to = "iter", values_to = "value")
# xi_red_plot %>% ggplot(aes(x = Stratum, y = value, group = Stratum, 
#                            fill = Stratum)) + 
#   theme_bw() + 
#   geom_boxplot() +
#   facet_grid(.~Class, labeller = label_both) + 
#   ggtitle(as.expression(bquote("Parameter estimation for "~xi~" ")))
p3_ref <- xi_red_plot %>% filter(Covariate == "Ref") %>%
  ggplot(aes(x = Covariate_Level, y = value, fill = Class)) +
  theme_bw() + scale_fill_brewer(palette="Set2") + 
  geom_boxplot() + xlab("Reference")
p3_age <- xi_red_plot %>% filter(Covariate == "Age") %>%
  ggplot(aes(x = Covariate_Level, y = value, fill = Class)) +
  theme_bw() + scale_fill_brewer(palette="Set2") + 
  geom_boxplot() + xlab("Age")
p3_race <- xi_red_plot %>% filter(Covariate == "Race/Ethnicity") %>%
  ggplot(aes(x = Covariate_Level, y = value, fill = Class)) +
  theme_bw() + scale_fill_brewer(palette="Set2") + 
  geom_boxplot() + xlab("Race/Ethnicity")
# p3_educ <- xi_red_plot %>% filter(Covariate == "Education") %>%
#   ggplot(aes(x = Covariate_Level, y = value, fill = Class)) +
#   theme_bw() + scale_fill_brewer(palette="Set2") + 
#   geom_boxplot() + xlab("Education")
p3_smoke <- xi_red_plot %>% filter(Covariate == "Smoking") %>%
  ggplot(aes(x = Covariate_Level, y = value, fill = Class)) +
  theme_bw() + scale_fill_brewer(palette="Set2") + 
  geom_boxplot() + xlab("Smoking")
p3_physactive <- xi_red_plot %>% filter(Covariate == "Physical Activity") %>%
  ggplot(aes(x = Covariate_Level, y = value, fill = Class)) +
  theme_bw() + scale_fill_brewer(palette="Set2") + 
  geom_boxplot() + xlab("Physical Activity")
ggarrange(p3_ref, p3_age, p3_race, p3_smoke, p3_physactive, 
          common.legend = TRUE, legend = "right", ncol = 3, nrow = 2)

#### Plot Phi boxplots separatey for each covariate
Phi_red_plot <- xi_red_plot %>%
  mutate(value = pnorm(value))
Phi_red_plot %>% filter(Covariate == "Age") %>%
  ggplot(aes(x = Covariate_Level, y = value, fill = Class)) +
  theme_bw() + scale_fill_brewer(palette="Set2") + 
  geom_boxplot() + ylab("Hypertension Risk") + xlab("Age")
# Phi_red_plot %>% filter(Covariate == "Education") %>%
#   ggplot(aes(x = Covariate_Level, y = value, fill = Class)) +
#   theme_bw() + scale_fill_brewer(palette="Set2") + 
#   geom_boxplot() + ylab("Hypertension Risk") + xlab("Education")
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
                     Phi = res$analysis_adj$Phi_med_adj)
### Plot marginal Phi
df_Phi %>% ggplot(aes(x = Class, y = Phi, col = Class)) +
  theme_bw() +
  geom_point() +
  ylim(0, 1) +
  ggtitle("Parameter estimation for conditional outcome probabilities") +
  ylab("P(Y=1|Class)")

df_Phi %>% ggplot(aes(x = Class, y = Phi, fill = Class)) +
  theme_bw() +
  geom_violin() + geom_boxplot(width = 0.1) +
  ylim(0, 1) +
  ggtitle("Parameter estimation for conditional outcome probabilities") +
  ylab("Estimated P(Y=1|Class)")

### Plot Phi with additional covariates
df_Phi <- data.frame(pnorm(res$analysis_adj$xi_med_adj))
colnames(df_Phi) <- c("Ref", age_categs[-1], racethnic_categs[-1], 
                      # educ_categs[-1], 
                      smoker_categs[-1],
                      physactive_categs[-1])
  # colnames(df_Phi) <- c(age_categs, educ_categs, racethnic_categs, smoker_categs)
# colnames(df_Phi) <- colnames(res$data_vars$V)
df_Phi$Class <- factor(1:K)
df_Phi <- df_Phi %>% gather("Covariate_Level", "Phi", -Class)
df_Phi$Covariate <- c(rep("Ref", K),
                      rep("Age", K*(length(age_categs) - 1)), 
                      rep("Race/Ethnicity", K*(length(racethnic_categs) - 1)), 
                      # rep("Education", K*(length(educ_categs) - 1)), 
                      rep("Smoking", K*(length(smoker_categs) - 1)),
                      rep("Physical Activity", K*(length(physactive_categs) - 1)))
  # df_Phi$Covariate <- c(rep("Age", K*6), rep("Education", K*3), 
  #                       rep("Race/Ethnicity", K*5), rep("Smoking", K*2))
# df_Phi %>% ggplot(aes(x = Covariate_Level, y = Phi, col = Class, group = Class)) +
#   theme_bw() + scale_color_brewer(palette="Set2") + 
#   geom_point() + geom_line() + 
#   facet_grid(~Covariate) +
#   ggtitle("Parameter estimation for conditional outcome probabilities") +
#   ylab("P(Y=1|-")
p_age <- df_Phi %>% filter(Covariate == "Age") %>%
  ggplot(aes(x = Covariate_Level, y = Phi, col = Class, group = Class)) +
  theme_bw() + scale_color_brewer(palette="Set2") + 
  geom_point() + geom_line() + ylim(0,1) + 
  ylab("Hypertension Probability") + xlab("Age") +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10))
# p_educ <- df_Phi %>% filter(Covariate == "Education") %>%
#   ggplot(aes(x = Covariate_Level, y = Phi, col = Class, group = Class)) +
#   theme_bw() + scale_color_brewer(palette="Set2") + 
#   geom_point() + geom_line() +  ylim(0,1) + 
#   ylab("Hypertension Probability") + xlab("Education") +
#   theme(axis.text=element_text(size=10),
#         axis.title=element_text(size=10))
p_race <- df_Phi %>% filter(Covariate == "Race/Ethnicity") %>%
  mutate(Covariate_Level = factor(Covariate_Level, levels = racethnic_categs)) %>%
  ggplot(aes(x = Covariate_Level, y = Phi, col = Class, group = Class)) +
  theme_bw() + scale_color_brewer(palette="Set2") + 
  geom_point() + geom_line() +  ylim(0,1) + 
  ylab("Hypertension Probability") + xlab("Race/Ethnicity") +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10)) +
  ggtitle("Hypertension Risk by Dietary Pattern and Race/Ethnicity")
p_smoker <- df_Phi %>% filter(Covariate == "Smoking") %>%
  ggplot(aes(x = Covariate_Level, y = Phi, col = Class, group = Class)) +
  theme_bw() + scale_color_brewer(palette="Set2") + 
  geom_point() + geom_line() +  ylim(0,1) + 
  ylab("Hypertension Probability") + xlab("Smoking") +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10))
p_physactive <- df_Phi %>% filter(Covariate == "Physical Activity") %>%
  ggplot(aes(x = Covariate_Level, y = Phi, col = Class, group = Class)) +
  theme_bw() + scale_color_brewer(palette="Set2") + 
  geom_point() + geom_line() +  ylim(0,1) + 
  ylab("Hypertension Probability") + xlab("Physical Activity") +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10))
ggarrange(p_age, p_race, p_smoker, p_physactive, ncol = 3, nrow=2,
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
                "Potatoes", "Other starchy vegs", "Other vegs",  
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
V_data <- data_vars_f_low %>% 
  select(RIDAGEYR, RIDRETH3, DMDEDUC2, smoker, Phys_Active) %>%
  mutate(
    age_cat = factor(case_when(  
      20 <= RIDAGEYR & RIDAGEYR <= 39 ~ 1,
      40 <= RIDAGEYR & RIDAGEYR <= 59 ~ 2,
      RIDAGEYR >= 60 ~ 3)),
    racethnic = factor(case_when(
      RIDRETH3 == 3 ~ 1,  # NH White
      RIDRETH3 == 4 ~ 2,  # NH Black
      RIDRETH3 == 6 ~ 3,  # NH Asian
      RIDRETH3 %in% c(1, 2) ~ 4,  # Mexican-American/Other Hispanic
      RIDRETH3 == 7 ~ 5,  # Other/Mixed
      .default = NA)),
    smoker = factor(smoker),
    Phys_Active = factor(Phys_Active),
    .keep = "unused")
res_demog <- V_data %>%
  mutate(racethnic = factor(racethnic),
         smoker = factor(smoker),
         # educ = factor(educ),
         age_cat = factor(age_cat))
res_demog$Class <- factor(res$analysis_adj$c_all)
res_demog$y <- res$data_vars$y_all
# p_educ <- res_demog %>% ggplot(aes(x = educ, fill = Class, group = Class)) +
#   theme_bw() + scale_fill_brewer(palette="Set2") + 
#   geom_bar(position = "dodge") + 
#   scale_x_discrete(labels = educ_categs) + 
#   ggtitle("Dietary Patterns by Education") + xlab("Education") + ylab("Count")
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
p_physactive <- res_demog %>% ggplot(aes(x = Phys_Active, fill = Class, group = Class)) +
  theme_bw() + scale_fill_brewer(palette="Set2") + 
  geom_bar(position = "dodge") + 
  scale_x_discrete(labels = smoker_categs) + 
  ggtitle("Dietary Patterns by Physical Activity") + xlab("Physical Activity") + 
  ylab("Count")
ggarrange(p_age, p_smoker, p_race, p_physactive, ncol = 2, nrow=2, common.legend = TRUE,
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
