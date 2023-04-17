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

#========================= MAIN FUNCTION =======================================

# 'WSOLCA_app_covs_Rcpp' runs the WSOLCA model with covariates and saves and 
# returns results
# Inputs:
#   data_path: String path for input dataset
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
WSOLCA_app_covs_Rcpp <- function(data_path, res_path, adj_path, stan_path, 
                            save_res = TRUE, n_runs, burn, thin) {
  start_time <- Sys.time()
  
  #================= Read in data ==============================================
  print("Read in data")
  data_vars_all <- read.csv(data_path)
  # Filter for low-income women and complete cases: n = 3028
  data_vars_f_low <- data_vars_all %>% 
    filter(sex == "female", poverty == "At or Below") %>%
    select(SEQN, RIDAGEYR, racethnic, EDU, SDMVPSU, SDMVSTRA, obese, smoker,
           bp.flag, sumrisk, dietwt8yr, bb1:bb29) %>% na.omit()
  
  # Obtain data
  x_mat <- as.matrix(data_vars_f_low %>% select(bb1:bb29))
  y_all <- data_vars_f_low$bp.flag
  s_all <- data_vars_f_low$SDMVPSU
  s_mat <- dummy_cols(data.frame(s = factor(s_all)),  # Stratifying variable as dummy columns
                      remove_selected_columns = TRUE)
  other_covs <- data_vars_f_low %>% 
    select(RIDAGEYR, racethnic, EDU, smoker) %>%
    mutate(age_cent = RIDAGEYR - mean(RIDAGEYR, na.rm = TRUE),
           racethnic = dummy_cols(data.frame(racethnic = factor(racethnic)),
                                  remove_selected_columns = TRUE),
           EDU = dummy_cols(data.frame(EDU = factor(EDU)),
                            remove_selected_columns = TRUE),
           smoker = dummy_cols(data.frame(smoker = factor(smoker)),
                               remove_selected_columns = TRUE))
  # Regression design matrix without class assignment, nxq
  V <- as.matrix(cbind(s_mat, other_covs))                
  sample_wt <- data_vars_f_low$dietwt8yr
  data_vars <- data.frame(x_mat, y_all, s_all, sample_wt)
  
  # Obtain dimensions
  n <- dim(x_mat)[1]        # Number of individuals
  p <- dim(x_mat)[2]        # Number of exposure items
  d <- max(apply(x_mat, 2,  # Number of exposure categories 
                 function(x) length(unique(x))))  # CHANGE TO ADAPT TO ITEM
  S <- length(unique(s_all))           # Number of strata
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
  K_med <- round(median(rowSums(MCMC_out$pi_MCMC >= 0.05)))
  # Get number of unique classes for fixed sampler
  K_fixed <- K_med
  print(paste0("K_fixed: ", K_fixed))
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
  res_MCMC <- list(MCMC_out = MCMC_out, post_MCMC_out = post_MCMC_out,
                   analysis = analysis)
  if (save_res) {
    save(res_MCMC, file = res_path)  # Save MCMC output
  }
  
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
                             s_all = s_all)
  
  runtime <- Sys.time() - start_time
  
  #================= Save and return output ====================================
  res <- list(analysis_adj = analysis_adj, runtime = runtime, 
              data_vars = data_vars, MCMC_out = MCMC_out)
  if (save_res) {
    save(res, file = adj_path)
  }
  return(res)
}



#===================== RUN MAIN WSOLCA FUNCTION =================================

# Define directories
wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/"
# wd <- "~/Documents/Harvard/Research/Briana/supRPC/wsOFMM/"
data_dir <- "Data/"
res_dir <- "Results/"
model_dir <- "Model_Code/"
model <- "wsOFMM"

# Define paths
# REMOVE ITER_POP
data_path <- paste0(wd, data_dir, "nhanesallage_frscores1118.csv")   # Input dataset
res_path <- paste0(wd, res_dir, model, "_results_nhanes", ".RData")  # Output file
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
  # Set seed
  set.seed(20230225)
  print(paste0("Running WSOLCA_application..."))
  results_adj <- WSOLCA_app_Rcpp(data_path = data_path, res_path = res_path,
                                 adj_path = adj_path, stan_path = stan_path, 
                                 save_res = TRUE, n_runs = 25000, burn = 15000, 
                                 thin = 5)
  print(paste0("Runtime: ", results_adj$runtime))
}



#===================== PLOT OUTPUT =============================================

load(adj_path)
#### pi
pi_red <- as.data.frame(res$analysis_adj$pi_red_adj)
colnames(pi_red) <- paste0("pi_", 1:(dim(pi_red)[2]))
pi_red_plot <- pi_red %>% pivot_longer(cols = everything(), names_to = "pi_comp", 
                                       values_to = "value")
pi_red_plot %>% ggplot(aes(x = pi_comp, y = value)) + 
  theme_bw() + 
  geom_boxplot() + 
  ggtitle(as.expression(bquote("Parameter estimation for "~pi~" "))) +
  xlab("Pi Component") + ylab("Value")

#### xi
xi_dims <- dim(res$analysis_adj$xi_red_adj)
xi_red <- as.data.frame(t(matrix(res$analysis_adj$xi_red_adj, 
                                 nrow = xi_dims[1], 
                                 ncol = xi_dims[2]*xi_dims[3])))
xi_red$Class <- factor(c(rep(1:xi_dims[2], times = xi_dims[3])))
xi_red$Stratum <- factor(c(rep(1:xi_dims[3], each = xi_dims[2])))
xi_red_plot <- xi_red %>% pivot_longer(cols = -c("Class", "Stratum"), 
                                       names_to = "iter", values_to = "value")
xi_red_plot %>% ggplot(aes(x = Class, y = value, group = Class, fill = Class)) + 
  theme_bw() + 
  geom_boxplot() + 
  facet_grid(.~Stratum, labeller = label_both) +
  ggtitle(as.expression(bquote("Parameter estimation for "~xi~" ")))
# xi_red_plot %>% ggplot(aes(x = Stratum, y = value, group = Stratum, 
#                            fill = Stratum)) + 
#   theme_bw() + 
#   geom_boxplot() +
#   facet_grid(.~Class, labeller = label_both) + 
#   ggtitle(as.expression(bquote("Parameter estimation for "~xi~" ")))

#### Phi
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

#### theta
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
  geom_text(aes(label = Level), col="white", cex=2) +
  scale_fill_gradient(trans="reverse") + 
  theme(legend.position="none") +
  ylab("Item") +
  ggtitle("Estimated modal item \nconsumption levels")

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


# # Testing
# results_adj <- WSOLCA_app_Rcpp(data_path = data_path, res_path = res_path,
#                                 adj_path = adj_path, stan_path = stan_path, 
#                                 save_res = FALSE, n_runs = 60, burn = 30, 
#                                 thin = 3)
