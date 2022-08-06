#====================================================
# Weighted Supervised Overfitted Finite Mixture Model
# Stephanie Wu
# Revised 06/29/2022
#====================================================

library(R.matlab)
library(extraDistr)
library(MASS)

#=========================
# Example parameter setup
scenario <- 6
sim_n <- 1
samp_n <- 2
in_dir <- paste0(getwd(), "/")
out_dir <- in_dir

k_max <- 50
sp_k <- k_max
alpha <- rep(1/sp_k, k_max)
eta <- rep(1, 4)

mu_0 <- rnorm(n=4, mean=0, sd=1)
Sig_0 <- 1 / rgamma(n=4, shape=5/2, rate=2/5)
Sig_0 <- diag(Sig_0)


#=================================
# Set seed and load simulated data
set.seed(sim_n)
samp_data <- readMat(paste0(in_dir, "simdata_scen", scenario, 
                              "_iter", sim_n, 
                              "_samp", samp_n, ".mat"))

#===================
# Get data variables
data_vars <- wtd_get_data_vars_latent(samp_data)

#===============================================
# Initialize priors and variables for OFMM model
OFMM_params <- wtd_init_OFMM_params_latent(data_vars, k_max, alpha, eta)

#=================================================
# Initialize priors and variables for probit model
probit_params <- init_probit_params_latent(data_vars, k_max, q_dem, 
                                           mu_0, Sig_0, OFMM_params)


# wtd_get_data_vars_latent takes in sample data and outputs relevant variables for
# the wsOFMM model with the latent variable formulation
# Inputs: samp_data structural array with at least the following columns:
#   X_data: food intake data as a matrix
#   Y_data: outcome data as a vector
#   sample_wt: survey weights as a vector
# Outputs: data_vars structural array with the following fields:
#   food: matrix of food intake data
#   n: number of individuals
#   p: number of food items
#   d_max: max number of consumption levels over all items
#   d: vector of max number of levels for each food item; px1
#   y: vector of outcomes; nx1
#   wt: vector of survey weights; nx1
#   wt_kappa: vector of normalized weights; nx1
#   wt_kappa_mat: matrix of normalized weights, replicated across items; nxp
#   lin_idx: vector of linear indices for unique item-response combos; (n*p)x1
wtd_get_data_vars_latent <- function(samp_data) {
  food <- samp_data$X_data
  n <- nrow(food)
  p <- ncol(food)
  d_max <- max(c(food))
  d <- apply(food, 1, max)
  y <- samp_data$Y_data
  wt <- samp_data$sample_wt
  
  kappa <- sum(wt) / n
  wt_kappa <- wt / kappa
  wt_kappa_mat <- matrix(rep(wt_kappa, times=p), nrow=n, byrow=FALSE)
  
  data_vars <- list("food" = food, "n" = n, "p" = p, "d_max" = d_max, "d" = d,
                    "y" = y, "wt" = wt, "kappa" = kappa, "wt_kappa" = wt_kappa, 
                    "wt_kappa_mat" = wt_kappa_mat)
}


wtd_init_OFMM_params_latent <- function(data_vars, k_max, alpha, eta) {
  pi <- rdirichlet(1, alpha)
  
  c_i <- rcat(n, prob=pi)  # Draw from a Categ(pi) random variable n times
  n_ci <- numeric(k_max)
  for (k in 1:k_max) {
    n_ci[k] <- sum(data_vars$wt_kappa[c_i == k])
  }
  
  theta <- array(0, dim=c(data_vars$p, k_max, data_vars$d_max), 
                 dimnames=c("p","k_max", "d_max"))
  for (k in 1:k_max) {
    for (j in 1:data_vars$p) {
      d_j <- data_vars$d[j]
      theta[j, k, 1:d_j] <- rdirichlet(1, eta[1:d_j])
    }
  }
  
}











