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

set.seed(11152022)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#===================== Helper functions ========================================

## 'grad_par()' helper function nested in 'withReplicates()' to obtain gradient.
## Stan will pass warnings from calling 0 chains, but will still create an 
## out_stan object for the 'grad_log_prob()' method
grad_par <- function(pwts, svydata, stanmod, standata, par_stan, upars) {
  standata$weights <- pwts
  out_stan <- sampling(object = stanmod, data = standata, pars = par_stan,
                       chains = 0, iter = 0, refresh = 0)
  gradpar <- grad_log_prob(out_stan, upars)
  return(gradpar)
}

#' Helper function to apply matrix rotation
#' @param par unadjusted parameter estimates
#' @param par_hat unadjusted mean parameter estimates
#' @param R2R1 adjustment matrix
#' @return adjusted parameter estimates
DEadj <- function(par, par_hat, R2R1) {
  par_adj <- (par - par_hat) %*% R2R1 + par_hat
  par_adj <- as.vector(par_adj)
  return(par_adj)
}

# 'unconstrain' converts from constrained space to unconstrained space for one 
# row, given input arrays of MCMC parameter output
# Inputs:
#   i: row index
#   K: number of classes
#   stan_model: stan model
#   pi: MCMC matrix output for pi; M rows
#   theta: MCMC array output for theta; dim 1 length = M
#   xi: MCMC matrix output for xi: M rows
# Output: vector of unconstrained parameters
unconstrain <- function(i, K, stan_model, pi, theta, xi) {
  xi_mat <- matrix(xi[i,], byrow = FALSE, nrow = K)
  upars <- unconstrain_pars(stan_model, list("pi" = pi[i,], 
                                             "theta" = theta[i,,,], 
                                             "xi" = xi_mat))
  return(upars)
}

#================= Main post-processing adjustment function ====================
# 'coverage_adj' applies the post-processing adjustment to one iteration of the
# model parameters and saves the adjusted results in an .RData file
# Inputs:
#   analysis: model output
#   sim_samp: dataset
#   mod_stan: stan model
#   sim_adj_path: path for adjusted file
# Outputs: saves the adjusted results in an .RData file and outputs list with:
#   analysis: model output with adjusted estimates replacing original estimates
#   neg_V1_Hi: vector of length 2 taking on value 1 if any negative variances 
#              V1 and Hi respectively
coverage_adj <- function(analysis, sim_samp, mod_stan, sim_adj_path) {
  
  #================= Create Stan data ==========================================
  
  # Posterior estimates from MCMC sampler
  K <- c(analysis$k_red)
  p <- dim(analysis$theta_med)[1]
  d <- dim(analysis$theta_med)[3]
  M <- dim(analysis$theta_red)[1]
  n <- length(analysis$c_i)
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
  xi_mat <- matrix(analysis$xi_med, byrow = FALSE, nrow = K)  # KxS matrix
  # xi_mat <- matrix(colMeans(analysis$xi_red), byrow = FALSE, nrow = K)  # KxS matrix
  
  n_chains <- 1
  alpha <- rep(1, K)/K
  eta <- matrix(1, nrow=K, ncol=d)
  mu0 <- rep(0, q)
  Sig0 <- diag(rep(1, q), nrow=q, ncol=q)
  
  data_stan <- list(K = K, p = p, d = d, n = n, q = q, X = x_mat, y = y_all, 
                    c = c_all, V = V, weights = w_all, 
                    alpha = alpha, eta = eta, mu0 = mu0, Sig0 = Sig0)
  
  #=============== Run Stan model ==============================================
  
  # Stan parameters of interest
  par_stan <- c('pi', 'theta', 'xi')  # subset of parameters interested in
  
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
  xi_med_adj <- apply(xi_red_adj, c(2,3), median)
  # xi_med_adj <- apply(xi_red_adj, c(2,3), mean)
  
  # Replace estimates with adjusted versions
  analysis$pi_red <- pi_red_adj
  analysis$theta_red <- theta_red_adj
  analysis$xi_red <- xi_red_adj
  analysis$pi_med <- pi_med_adj
  analysis$theta_med <- theta_med_adj
  analysis$xi_med <- xi_med_adj
  
  # Save and return results
  save(analysis, file = sim_adj_path)
  return(list(analysis=analysis, neg_V1_Hi=neg_V1_Hi))
}

# ## Troubleshooting
# quantile(pi_red_adj[,1], c(0.025, 0.975))
# quantile(theta_red_adj[,1,3,4], c(0.025, 0.975))
# c(xi_med_adj)

#=== Function to read data and apply post-processing adjustment for all samples
# 'run_adj_samples' applies the post-processing adjustment to all iterations of
# the simulated data
# Inputs:
#   data_dir: directory where data is stored
#   res_dir: directory where model results are stored
#   analysis_dir: directory where analysis code is stored
#   iter_pop: numeric population iteration, usually 1
#   scen_samp: scenario denoting sampling design
#   model: model type, one of "sOFMM", "wsOFMM", or "wOFMM"
#   R_seq: numeric sequence listing simulation scenarios
#   mod_stan: stan model
run_adj_samples <- function(data_dir, res_dir, analysis_dir, iter_pop, scen_samp,
                            model, R_seq, mod_stan) {
  # Track which iterations have negative matrix diagonals
  # Col 1 is V1, col 2 is Hi
  neg_V1_Hi_var <- matrix(NA, nrow=max(R_seq), ncol=2)
  neg_V1_Hi_var[R_seq, ] <- 0
  
  for (i in 1:length(R_seq)) { 
    samp_n = R_seq[i]
    print(samp_n)
    # File name for adjusted output
    sim_adj_path <- paste0(res_dir, model, "_latent_results_adj_scen", scen_samp,
                           "_iter", iter_pop, "_samp", samp_n, ".RData")
    # sim_adj_path <- paste0(res_dir, "postprocess_mean", model, "_latent_results_adj_scen", scen_samp,
    #                        "_iter", iter_pop, "_samp", samp_n, ".RData")
    
    # Read in Matlab output
    sim_samp_path <- paste0(data_dir, "simdata_scen", scen_samp,"_iter", iter_pop, "_samp", samp_n, ".mat")
    sim_res_path <- paste0(res_dir, model, "_latent_results_scen", scen_samp, "_iter", iter_pop, "_samp", samp_n, ".mat")
    if (!file.exists(sim_samp_path) | !file.exists(sim_res_path)) {
      if (!file.exists(sim_samp_path)) {
        print(paste0("File does not exist: simdata_scen", scen_samp,"_iter", iter_pop, "_samp", samp_n))
      }
      if (!file.exists(sim_res_path)) {
        print(paste0("File does not exist: ", model, "_latent_results_scen", scen_samp, "_iter", iter_pop, "_samp", samp_n, ".mat"))
      }
    } else if (file.exists(sim_adj_path)) { # do nothing if adjusted file exists
      print(paste0("File already exists: ", model, "_latent_results_adj_scen", scen_samp, "_iter", iter_pop, "_samp", samp_n, ".RData"))
    } else {
      # Load simulated sample data for the iteration
      sim_samp <- readMat(sim_samp_path)$sim.data
      names(sim_samp) <- str_replace_all(dimnames(sim_samp)[[1]], "[.]", "_")
      
      # Load model output and extract analysis portion
      output <- readMat(sim_res_path)
      analysis <- output$analysis
      names(analysis) <- str_replace_all(dimnames(analysis)[[1]], "[.]", "_")
      
      adjusted <- coverage_adj(analysis=analysis, sim_samp=sim_samp, 
                               mod_stan=mod_stan, sim_adj_path=sim_adj_path)
      neg_V1_Hi_var[samp_n, ] <- adjusted$neg_V1_Hi
    }
  }
  return(neg_V1_Hi_var)
}


#===== Read in data and apply post-processing adjustment for all samples========

setwd("/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/")
#setwd("~/Documents/Harvard/Research/Briana/supRPC/wsOFMM")
#setwd("/Users/Stephanie/Documents/GitHub/wsOFMM")
data_dir <- "Data/"
res_dir <- "Results/"
analysis_dir <- "Analysis_Code/"
iter_pop <- 1
scen_samp <- 101
model <- "wsOFMM"
R_seq=1:100

#samp_n <- 1
# Create Stan model
mod_stan <- stan_model(paste0(analysis_dir, "mixture_model.stan"))

# SRS: Scenario 101
start_time101 <- Sys.time()
out101 <- run_adj_samples(data_dir=data_dir, res_dir=res_dir, analysis_dir=analysis_dir, 
                          iter_pop=iter_pop, scen_samp=scen_samp, model=model, 
                          R_seq=R_seq, mod_stan = mod_stan)
end_time101 <- Sys.time()
# Stratified: Scenario 201
start_time201 <- Sys.time()
out201 <- run_adj_samples(data_dir=data_dir, res_dir=res_dir, analysis_dir=analysis_dir, 
                          iter_pop=iter_pop, scen_samp=201, model=model, 
                          R_seq=R_seq, mod_stan = mod_stan)
end_time201 <- Sys.time()

print(paste0("Runtime SRS: ", end_time101 - start_time101))
colSums(out101, na.rm = TRUE)
print(out101)
print(paste0("Runtime Strat: ", end_time201 - start_time201))
colSums(out201, na.rm = TRUE)
print(out201)


#=============== MISC EXTRA CODE ===============================================
#temp_par_samps <- rstan::extract(out_stan, pars = par_stan, permuted = FALSE)

#pars <- constrain_pars(out_stan, upars)
#grad_constr <- grad_log_prob(out_stan, upars = upars, adjust_transform = FALSE)
# grad_unconstr <- grad_log_prob(out_stan, upars = upars)
# grad_constr <- constrain_pars(out_stan, grad_unconstr)

# HHat_list <- list()
# for (i in 1500:2000) {
#   temp_upars <- unconstrain_pars(out_stan, list("pi" = c(analysis$pi_red[i,]),
#                                                "theta" = analysis$theta_red[i,,,],
#                                                "xi" = analysis$xi_red[i,]))
#   temp_HHat <- -1*optimHess(temp_upars, gr = function(x){grad_log_prob(out_stan, x)})
#   HHat_list[[i]] <- temp_HHat
# }
# HHat_MCMC <- Reduce("+", HHat_list) / length(HHat_list)
# Hi_MCMC <- solve(HHat_MCMC)


# # aaply: split array by rows (margins=1) and apply function DEadj to get 
# # adjusted version of each piece. Resulting array has same dimensions as input: 
# # (#iter)x(#chains)x(#params)
# # unc_par_samps: posterior estimates from the MCMC samples in unconstrained space
# # unc_par_hat: posterior point estimates in unconstrained space; parameter for DEadj function
# # R2R1: projection matrix to apply adjustment; parameter for DEadj function
# par_adj <- aaply(unc_par_samps, 1, DEadj, par_hat = unc_par_hat, R2R1 = R2R1, 
#                  .drop = FALSE)

# # 'constrain' converts from unconstrained space to constrained space for one 
# # row, given input MCMC matrix of unconstrained adjusted parameters
# # Inputs:
# #   i: row index
# #   stan_model: stan model
# #   upars_mat: MCMC matrix of unconstrained parameters 
# # Output: vector of constrained parameters
# constrain <- function(i, stan_model, upars_mat) {
#   constr_pars <- constrain_pars(stan_model, upars_mat[i,])
#   return(list(pi = constr_pars$pi, theta = constr_pars$theta, xi = constr_pars$xi))
# }
# 
# constr_par_samps <- lapply(1:M, constrain, stan_model = out_stan, upars_mat = par_adj)
