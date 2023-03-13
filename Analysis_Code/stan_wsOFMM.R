# Running WSOLCA using Stan
# Author: Stephanie Wu
# Date updated: 2023/03/08

library(survey)
library(tidyverse)
library(R.matlab)
library(plyr)
library(fastDummies)
# # Uncomment these lines if mixture_model.stan has syntax errors
# remove.packages(c("StanHeaders", "rstan"))
# if (file.exists(".RData")) file.remove(".RData")
# install.packages("StanHeaders", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(rstan)
library(abind)
library(label.switching)

set.seed(11152022)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


#setwd("/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/")
setwd("~/Documents/Harvard/Research/Briana/supRPC/wsOFMM")
#setwd("/Users/Stephanie/Documents/GitHub/wsOFMM")


#================= FIXED SAMPLER ===============================================

#================= Create Stan data ============================================
# Read in data
data_dir <- "Data/"
res_dir <- "Results/"
analysis_dir <- "Analysis_Code/"
iter_pop <- 1
scen_samp <- 101
model <- "wsOFMM"
R_seq=1:100
samp_n <- 84

sim_samp_path <- paste0(data_dir, "simdata_scen", scen_samp,"_iter", iter_pop, 
                        "_samp", samp_n, ".mat")
sim_res_path <- paste0(res_dir, model, "_latent_results_scen", scen_samp, 
                       "_iter", iter_pop, "_samp", samp_n, ".mat")
already_done <- file.exists(sim_res_path)

# Load simulated sample data for the iteration
sim_samp <- readMat(sim_samp_path)$sim.data
names(sim_samp) <- str_replace_all(dimnames(sim_samp)[[1]], "[.]", "_")

# Posterior estimates from MCMC sampler
K <- 3
p <- dim(sim_samp$X_data)[2]
d <- 4  # CHANGE TO ADAPT TO ITEM
n <- dim(sim_samp$X_data)[1]
x_mat <- sim_samp$X_data
y_all <- c(sim_samp$Y_data)
s_all <- sim_samp$true_Si
S <- length(unique(s_all))
s_mat <- dummy_cols(data.frame(s = factor(s_all)), 
                    remove_selected_columns = TRUE)
V <- as.data.frame(s_mat)
q <- S  # number of additional covariates other than class assignment
w_all <- c(sim_samp$sample_wt / sum(sim_samp$sample_wt) * n)
c_all <- c(analysis$c_i)

n_chains <- 4
alpha <- rep(1, K)/K
eta <- matrix(1, nrow=K, ncol=d)
mu0 <- rep(0, q)
Sig0 <- diag(rep(1, q), nrow=q, ncol=q)

# Create Stan model
mod_stan <- stan_model(paste0(analysis_dir, "WSOLCA.stan"))

data_stan <- list(K = K, p = p, d = d, n = n, q = q, X = x_mat, y = y_all, 
                  V = V, weights = w_all, alpha = alpha, eta = eta, mu0 = mu0, 
                  Sig0 = Sig0)

#=============== Run Stan model ================================================

# # Stan parameters of interest
# par_stan <- c('pi', 'theta', 'xi')  # subset of parameters interested in

# Run Stan model
# Stan will pass warnings from calling 0 chains, but will still create an 
# out_stan object for the 'grad_log_prob()' method
start_time201 <- Sys.time()
out_stan <- sampling(object = mod_stan, data = data_stan, 
                     chains = 2, iter = 2500, warmup = 1500, thin = 5)
end_time201 <- Sys.time()

print(paste0("Runtime: ", end_time101 - start_time101))

#=============== Check convergence =============================================

print(out_stan, c("pi", "theta", "xi", "xi_prod", "theta_prod"))
summary(out_stan)$summary[, 9:10]
traceplot(out_stan, "pi")
traceplot(out_stan, "xi")
traceplot(out_stan, "theta")
summary(out_stan)$summary[-c(1:123), 9:10]

#=============== Post-hoc relabeling ===========================================
out_stan@model_pars
post_par <- rstan::extract(out_stan, c("pi", "theta", "xi", "pred_class_probs", "lp__"))
post_class_probs <- post_par$pred_class_probs  # M x n x K
post_class_probs[1:5, 15, ]
# post_class_temp <- ((post_class_probs[,,1] > 0.5) * 1) + 1
post_class <- apply(post_class_probs, c(1,2), function(x) which.max(x))  # M x n

M <- dim(post_class)[1]
num_comp_params <- (K + (p*K*d) + (S*K)) / K  # 123
# Initialize mcmc arrays
mcmc <- array(NA, dim = c(M, K, num_comp_params))
# Assign posterior draws to the array
mcmc[, , 1] <- pi
for (r in 1:d) {
  for (j in 1:p) {
    mcmc[, , 1 + (r-1)*p + j] <- theta[, j, , r]
  }
}
for (s in 1:S) {
  mcmc[, , (num_comp_params - S) + s] <- xi[, , s]
}

# Set of selected relabeling algorithms
set <-
  c("PRA",
    "ECR",
    "ECR-ITERATIVE-1",
    "AIC",
    "ECR-ITERATIVE-2",
    "STEPHENS")  # took out data-based

# Find the MAP draw as a pivot
mapindex = which.max(post_par$lp__)

# switch labels
label_switch <- label.switching(method = set, zpivot = post_class[mapindex,],
    z = post_class, K = K, prapivot = mcmc[mapindex, ,], constraint = 1,
    mcmc = mcmc, p = post_class_probs)

label_switch$timings
label_switch$similarity

# Permute posterior based on ECR method
mcmc_permuted <- permute.mcmc(mcmc, label_switch$permutations$ECR)$output # MxKxnum_comp_pars

# Change dimension for each parameter to match order of Stan parameters
change_mcmc_dim <- function(mcmc, num_comp_params, M, p, K, d, S, out_stan) {
  param_names <- out_stan %>% names %>% `[` (1:(num_comp_params*K))
  
  mcmc_full_dim <- array(mcmc[, , 2:(num_comp_params - S)], c(M, K, p, d))
  mcmc_flat <- cbind(matrix(mcmc[, , 1], nrow = M), 
                     matrix(aperm(mcmc_full_dim, c(1,3,2,4)), nrow = M), 
                     matrix(mcmc[, , -(1:(num_comp_params - S))], nrow = M))
  mcmc_flat <- array(mcmc_flat, dim = c(M, 1, ncol(mcmc_flat)), 
                dimnames = list(NULL, NULL, param_names))
  return(mcmc_flat)
}


mcmc_permuted_2 <- change_mcmc_dim(mcmc_permuted, num_comp_params, 
                                 M, p, K, d, S, out_stan)

# Reassess model convergence after relabelling
fit_permuted <- monitor(mcmc_permuted_2, warmup = 0, digits_summary = 3)

par_summary <- as.data.frame(fit_permuted)
par_summary[1:250, c(1, 9, 10)]
par_summary[-(1:250), c(1, 9, 10)]

# Troubleshooting and sanity checks
library(testthat)
# Sanity checks
post_par$pi[1,] == mcmc_flat[1, 1:3]
post_par$theta[1, 5, 1, 1] == mcmc_flat[1, 3 + 5]
post_par$theta[1, 1, 2, 1] == mcmc_flat[1, 3 + 30 + 1]
post_par$theta[1, 1, 1, 2] == mcmc_flat[1, 3 + 30*3 + 1]
post_par$xi[1, 2, 1] == mcmc_flat[1, (3*121) + 2] 
post_par$xi[1, 1, 2] == mcmc_flat[1, (3*121) + 4] 


test_mcmc <- array(NA, dim=c(5, 3, 11))
# for (j in 1:2) {
#   for (k in 1:3) {
#     for (r in 1:4) {
#       test_mcmc[, k, (j-1)*d + r] <- paste0(j,"_",k, "_", r)
#     }
#   }
# }
test_mcmc[,,1] <- rep(paste0("pi_", 1:3), each = 5)
for (r in 1:4) {
  for (k in 1:3) {
    for (j in 1:2) {
      test_mcmc[, k, 1+ (r-1)*2 + j] <- paste0(j,"_",k, "_", r)
    }
  }
}
test_mcmc[, , 10] <- rep(paste0("xi_s1_", 1:3), each = 5)
test_mcmc[, , 11] <- rep(paste0("xi_s2_", 1:3), each = 5)
array(test_mcmc, dim = c(5, 3*11))

mod <- array(test_mcmc[, , -c(1, 10, 11)], c(5, K, 2, 4))
cbind(matrix(test_mcmc[, , 1], nrow = 5), 
      matrix(aperm(mod, c(1,3,2,4)), nrow = 5), 
      matrix(test_mcmc[, , c(10, 11)], nrow = 5))

mod <- array(test_mcmc[, , 2:9], dim=c(5, 3, 2, 4))
mod_flat <- array(aperm(mod, c(1, 3, 2, 4)), dim=c(5, 3*2*4))
test_mcmc_flat <- array(test_mcmc, dim=c(5, 3*11))
test_mcmc_flat[, (K+1):9] <- mod_flat
test_mcmc_flat
# Check matches names output
array(aperm(mod, c(1, 3, 2, 4)), dim=c(5, 3*2*4))[1,]
pars[4:70]
# array(aperm(test_mcmc, c(2, 1, 3)), dim=c(5, 3*8))
# array(test_mcmc, dim=c(5, 3*8))

get_indices <- function(j, r, d) {
  return(1 + (j-1)*d + r)
}
test_that(desc = "Correct indices iter 17", code = {
  expect_that(get_indices(1, 1, 4), equals(2))
  expect_that(get_indices(1, 3, 4), equals(4))  
  expect_that(get_indices(3, 2, 4), equals(11))
  expect_that(get_indices(6, 1, 4), equals(22))
})
get_indices_2 <- function(j, r, p) {
  return(1 + (r-1)*p + j)
}
test_that(desc = "Correct indices iter 17", code = {
  expect_that(get_indices_2(1, 1, 30), equals(2))
  expect_that(get_indices_2(1, 3, 30), equals(62))  
  expect_that(get_indices_2(3, 2, 30), equals(34))
  expect_that(get_indices_2(6, 1, 30), equals(7))
})
abs(rowSums(mcmc[1:5, , 1]) - 1) < 0.00001
abs(apply(mcmc[1:5, , 2:5], c(1,2), sum) - 1) < 0.00001



# Run Stan model
# Stan will pass warnings from calling 0 chains, but will still create an 
# out_stan object for the 'grad_log_prob()' method
out_stan <- sampling(object = mod_stan, data = data_stan, pars = par_stan,
                     chains = 0, iter = 0, refresh = 0)

#=============== Convert to unconstrained parameters =========================

# # get dimension of unconstrained parameter space
# get_num_upars(out_stan)  #278 unconstrained, 369 constrained
# convert params from constrained space to unconstrained space
unc_par_hat <- unconstrain_pars(out_stan, list("pi" = c(analysis$pi_med),
                                               "theta" = analysis$theta_med,
                                               "xi" = xi_mat))
# Unconstrained parameters for all MCMC samples
M <- dim(analysis$theta_red)[1]
unc_par_samps <- lapply(1:M, unconstrain, stan_model = out_stan, K = K, 
                        pi = analysis$pi_red, theta = analysis$theta_red, 
                        xi = analysis$xi_red)
unc_par_samps <- matrix(unlist(unc_par_samps), byrow = TRUE, nrow = M)

#=============== Post-processing adjustment in unconstrained space ===========

# Estimate Hessian
Hhat <- -1*optimHess(unc_par_hat, gr = function(x){grad_log_prob(out_stan, x)})

# Create survey design
svy_data <- data.frame(s = s_all, 
                       x = x_mat,
                       y = y_all, 
                       wts = w_all)
svydes <- svydesign(ids = ~1, strata = ~s, weights = ~wts, data = svy_data)
# create svrepdesign
svyrep <- as.svrepdesign(design = svydes, type = "mrbbootstrap", 
                         replicates = 100)

rep_temp <- withReplicates(design = svyrep, theta = grad_par, 
                           stanmod = mod_stan, standata = data_stan, 
                           par_stan = par_stan, upars = unc_par_hat)
Jhat <- vcov(rep_temp)

# compute adjustment
Hi <- solve(Hhat)
V1 <- Hi %*% Jhat %*% Hi

# Check for issues with negative diagonals
neg_V1_Hi <- numeric(2)
if (min(diag(V1)) < 0) {
  neg_V1_Hi[1] <- 1
  print("V1 has negative variances")
}
if (min(diag(Hi)) < 0) {
  neg_V1_Hi[2] <- 1
  print("Hi has negative variances")
}
# If matrices are not p.d. due to rounding issues, convert to nearest p.d. 
# matrix using method proposed in Higham (2002)
if (min(Re(eigen(V1)$values)) < 0) { 
  V1_pd <- nearPD(V1)
  R1 <- chol(V1_pd$mat)
  print(paste0("V1: absolute eigenvalue difference to nearest p.d. matrix: ", 
               sum(abs(eigen(V1)$values - eigen(nearPD(V1)$mat)$values))))
} else {
  R1 <- chol(V1)
}
if (min(Re(eigen(Hi)$values)) < 0) {
  Hi_pd <- nearPD(Hi)
  R2_inv <- chol(Hi_pd$mat)
  print(paste0("Hi: absolute eigenvalue difference to nearest p.d. matrix: ", 
               sum(abs(eigen(Hi)$values - eigen(nearPD(Hi)$mat)$values))))
} else {
  R2_inv <- chol(Hi)
}
# Obtain the variance adjustment matrix
R2 <- solve(R2_inv)
R2R1 <- R2 %*% R1

# Combine chains together
dim(unc_par_samps) <- c(M, n_chains*length(unc_par_hat))

# Adjust samples
# unc_par_samps: posterior estimates from the MCMC samples in unconstrained space
par_adj <- apply(unc_par_samps, 1, DEadj, par_hat = unc_par_hat, R2R1 = R2R1, 
                 simplify = FALSE)
par_adj <- matrix(unlist(par_adj), byrow=TRUE, nrow=M)

#=============== Convert adjusted to constrained space =======================

# Constrained adjusted parameters for all MCMC samples
pi_red_adj <- matrix(NA, nrow=M, ncol=K)
theta_red_adj <- array(NA, dim=c(M, p, K, d))
xi_red_adj <- array(NA, dim=c(M, K, q))
for (i in 1:M) {
  constr_pars <- constrain_pars(out_stan, par_adj[i,])
  pi_red_adj[i, ] <- constr_pars$pi
  theta_red_adj[i,,,] <- constr_pars$theta
  xi_red_adj[i,,] <- constr_pars$xi
}

#=============== Output adjusted parameters ==================================

pi_med_adj <- apply(pi_red_adj, 2, median)
theta_med_adj <- apply(theta_red_adj, c(2,3,4), median)
# xi_med_adj <- apply(xi_red_adj, c(2,3), median)
xi_med_adj <- apply(xi_red_adj, c(2,3), mean)