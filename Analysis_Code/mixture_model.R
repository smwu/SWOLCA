##First read in the arguments listed at the command line
args = commandArgs(trailingOnly=TRUE)
##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
  print("No arguments supplied.")
  R_seq <- 1  # sample iteration
}else{
  R_seq <- args[1]
}
print(paste0("Sample: ", R_seq))

library(survey)
library(tidyverse)
library(R.matlab)
library(plyr)
#remove.packages(c("StanHeaders", "rstan"))
#install.packages("StanHeaders", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
#install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(rstan)
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

## 'DEadj()' helper function to apply matrix rotation
DEadj <- function(par, par_hat, R2R1) {
  par_adj <- (par - par_hat) %*% R2R1 + par_hat
  par_adj <- as.vector(par_adj)
  return(par_adj)
}

# 'unconstrain' converts from constrained space to unconstrained space for one 
# row, given input arrays of MCMC parameter output
# Inputs:
#   i: row index
#   stan_model: stan model
#   pi: MCMC matrix output for pi; M rows
#   theta: MCMC array output for theta; dim 1 length = M
#   xi: MCMC matrix output for xi: M rows
# Output: vector of unconstrained parameters
unconstrain <- function(i, stan_model, pi, theta, xi) {
  upars <- unconstrain_pars(stan_model, list("pi" = pi[i,], 
                                             "theta" = theta[i,,,], "xi" = xi[i,]))
  return(upars)
}

#================= Main post-processing adjustment function ====================
coverage_adj <- function(analysis, sim_samp, mod_stan, sim_adj_path) {
  
  #================= Create Stan data ==========================================
  
  # Posterior estimates from MCMC sampler
  K <- c(analysis$k_red)
  p <- dim(analysis$theta_med)[1]
  d <- dim(analysis$theta_med)[3]
  n <- length(analysis$c_i)
  q <- length(analysis$xi_med)
  X <- sim_samp$X_data
  y <- c(sim_samp$Y_data)
  S <- length(unique(sim_samp$true_Si))
  weights <- c(sim_samp$sample_wt / sum(sim_samp$sample_wt) * n)
  n_chains <- 1
  
  alpha <- rep(1, K)/K
  eta <- matrix(1, nrow=K, ncol=d)
  mu0 <- rep(0, q)
  Sig0 <- diag(rep(1, q), nrow=q, ncol=q)
  
  # Change xi to reference cell coding 
  # No need to adjust for label switching
  #xi_ref <- analysis$xi_med[1]
  xi_red_ref <- analysis$xi_red
  for (v in 2:q) {
    if (v < (K+S)) {
      #xi_ref[v] <- analysis$xi_med[v] - analysis$xi_med[1]
      xi_red_ref[, v] <- analysis$xi_red[, v] - analysis$xi_red[, 1]
    } else {
      #xi_ref[v] <- analysis$xi_med[v] - analysis$xi_med[K+S-1] - xi_ref[v - K]
      xi_red_ref[, v] <- analysis$xi_red[, v] - analysis$xi_red[, K+S-1] - xi_red_ref[, v - K]
    }
  }
  xi_ref <- apply(xi_red_ref, 2, median)
  
  # Get posterior MCMC chains as a matrix with dimensions (M)x(# parameters)
  # Converting theta array to vector order: dim j -> dim k -> dim r
  # First 1:p for k=1 and d=1, then 1:p for k=2 and d=1, then 1:p for k=3 and d=1, 
  # then 1:p for k=1 and d=2, then 1:p for k=2 and d=2, then 1:p for k=3 and d=2,... 
  M <- dim(analysis$pi_red)[1]
  theta_vectorized <- analysis$theta_red
  dim(theta_vectorized) <- c(M, p*K*d)
  par_samps <- cbind(analysis$pi_red, theta_vectorized, xi_red_ref)
  colnames(par_samps) <- c(paste0("pi", 1:K), 
                           paste0("theta", 1:p, "_", rep(1:K, each=p), "_", rep(1:d, each=p*K)),
                           paste0("xi", 1:q))
  
  # Get mean posterior parameter estimates
  par_hat_mean <- colMeans(par_samps)
  
  # Get median posterior parameter estimates
  par_hat <- c(analysis$pi_med, c(analysis$theta_med), xi_ref)
  names(par_hat) <- c(paste0("pi", 1:K), 
                      paste0("theta", 1:p, "_", rep(1:K, each=p), "_", rep(1:d, each=p*K)),
                      paste0("xi", 1:q))
  
  # Define probit design matrix
  probit_data <- data.frame(s = factor(sim_samp$true_Si),
                            c = factor(analysis$c_i, levels=1:K))
  V <- model.matrix(~ c * s, data = probit_data)
  # Create array with assigned classes
  V_k <- array(NA, dim=c(K, n, q))
  for (k in 1:K) {
    temp_data <- data.frame(s = factor(sim_samp$true_Si), 
                            c = factor(rep(k, n), levels=1:K))
    V_k[k, ,] <- model.matrix(~ c * s, data = temp_data) 
  }
  
  data_stan <- list(K = K, p = p, d = d, n = n, q = q, X = X, y = y, V_k = V_k,
                    weights = weights, alpha = alpha, eta = eta, mu0 = mu0, Sig0 = Sig0)
  
  #=============== Run Stan model ==============================================
  
  # Stan parameters of interest
  par_stan <- c('pi', 'theta', 'xi')  # subset of parameters interested in
  
  # Run Stan model
  # Stan will pass warnings from calling 0 chains, but will still create an 
  # out_stan object for the 'grad_log_prob()' method
  out_stan <- sampling(object = mod_stan, data = data_stan, pars = par_stan,
                       chains = 0, iter = 0, refresh = 0)
  
  #=============== Convert to unconstrained parameters =========================
  
  # get dimension of unconstrained parameter space
  get_num_upars(out_stan)  #278 unconstrained, 369 constrained
  # convert params from constrained space to unconstrained space
  unc_par_hat <- unconstrain_pars(out_stan, list("pi" = c(analysis$pi_med),
                                                 "theta" = analysis$theta_med,
                                                 "xi" = xi_ref))
  # Unconstrained parameters for all MCMC samples
  unc_par_samps <- lapply(1:M, unconstrain, stan_model = out_stan,
                          pi = analysis$pi_red, theta = analysis$theta_red, xi = xi_red_ref)
  unc_par_samps <- matrix(unlist(unc_par_samps), byrow = TRUE, nrow = M)
  
  #=============== Post-processing adjustment in unconstrained space ===========
  
  # Estimate Hessian
  Hhat <- -1*optimHess(unc_par_hat, gr = function(x){grad_log_prob(out_stan, x)})
  
  # Estimate Jhat = Var(gradient)
  print('gradient evaluation')
  
  # Create survey design
  svy_data <- data.frame(s = sim_samp$true_Si, 
                         x = X,
                         y = y, 
                         wts = weights)
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
  # check issue with p.d. projection
  if (min(diag(V1)) < 0) {
    print("V1 has negative variances")
  }
  if (min(diag(Hi)) < 0) {
    print("Hi has negative variances")
  }
  # if matrices are not p.d. due to rounding issues, convert to nearest p.d. matrix
  # first, using method proposed in Higham (2002)
  if (min(Re(eigen(V1)$values)) < 1e-10) { 
    V1_pd <- nearPD(V1)
    R1 <- chol(V1_pd$mat)
  } else {
    R1 <- chol(V1)
  }
  if (min(Re(eigen(Hi)$values)) < 1e-10) {
    Hi_pd <- nearPD(Hi)
    R2i <- chol(Hi_pd$mat)
  } else {
    R2i <- chol(Hi)
  }
  # R1 <- chol(V1)
  # R2i <- chol(Hi)
  R2 <- solve(R2i)
  R2R1 <- R2 %*% R1
  
  # Combine chains together
  dim(unc_par_samps) <- c(M, n_chains*length(unc_par_hat))
  
  # Adjust samples
  # unc_par_samps: posterior estimates from the MCMC samples in unconstrained space
  par_adj <- apply(unc_par_samps, 1, DEadj, par_hat = unc_par_hat, R2R1 = R2R1, simplify = FALSE)
  par_adj <- matrix(unlist(par_adj), byrow=TRUE, nrow=M)
  
  #=============== Convert adjusted to constrained space =======================
  
  print("preparing output")
  
  # Constrained adjusted parameters for all MCMC samples
  pi_red_adj <- matrix(NA, nrow=M, ncol=K)
  theta_red_adj <- array(NA, dim=c(M, p, K, d))
  xi_red_ref_adj <- matrix(NA, nrow=M, ncol=q)
  for (i in 1:M) {
    constr_pars <- constrain_pars(out_stan, par_adj[i,])
    pi_red_adj[i, ] <- constr_pars$pi
    theta_red_adj[i,,,] <- constr_pars$theta
    xi_red_ref_adj[i, ] <- constr_pars$xi
  }
  
  # Change xi back to factor variable coding
  xi_red_adj <- xi_red_ref_adj
  for (v in 2:q) {
    if (v < (K+S)) {
      xi_red_adj[, v] <- xi_red_ref_adj[, 1] + xi_red_ref_adj[, v] 
    } else {
      xi_red_adj[, v] <- xi_red_adj[, K+S-1] + xi_red_ref_adj[, v - K] + xi_red_ref_adj[, v]
    }
  }
  
  #=============== Output adjusted parameters ==================================
  
  pi_mean_adj <- colMeans(pi_red_adj)
  theta_mean_adj <- apply(theta_red_adj, c(2,3,4), mean)
  xi_mean_adj <- colMeans(xi_red_adj)
  
  # adjusted <- list(pi_red_adj = pi_red_adj, theta_red_adj = theta_red_adj, 
  #                  xi_red_adj = xi_red_adj, 
  #                  pi_mean_adj = pi_mean_adj, theta_mean_adj = theta_mean_adj, 
  #                  xi_mean_adj = xi_mean_adj)
  
  # Replace estimates with adjusted versions
  analysis$pi_red <- pi_red_adj
  analysis$theta_red <- theta_red_adj
  analysis$xi_red <- xi_red_adj
  analysis$pi_med <- pi_mean_adj
  analysis$theta_med <- theta_mean_adj
  analysis$xi_med <- xi_mean_adj
  
  # Save and return results
  save(analysis, file = sim_adj_path)
  return(analysis)
}



#=== Function to read data and apply post-processing adjustment for all samples

run_adj_samples <- function(data_dir, res_dir, analysis_dir, iter_pop, scen_samp,
                            model, R_seq, mod_stan) {
  for (i in 1:length(R_seq)) { 
    samp_n = R_seq[i]
    print(samp_n)
    # File name for adjusted output
    sim_adj_path <- paste0(res_dir, model, "_latent_results_adj_scen", scen_samp, 
                           "_iter", iter_pop, "_samp", samp_n, ".RData")
    
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
    }
  } 
}

#===== Read in data and apply post-processing adjustment for all samples========
setwd("/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/")
#setwd("~/Documents/Harvard/Research/Briana/supRPC/wsOFMM")
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

start_time <- Sys.time()
run_adj_samples(data_dir=data_dir, res_dir=res_dir, analysis_dir=anaysis_dir, 
                iter_pop=iter_pop, scen_samp=scen_samp, model=model, R_seq=R_seq, 
                mod_stan = mod_stan)
end_time <- Sys.time()
print(paste0("Runtime SRS: ", end_time - start_time))

start_time <- Sys.time()
run_adj_samples(data_dir=data_dir, res_dir=res_dir, analysis_dir=anaysis_dir, 
                iter_pop=iter_pop, scen_samp=201, model=model, R_seq=R_seq, 
                mod_stan = mod_stan)
end_time <- Sys.time()
print(paste0("Runtime Strat: ", end_time - start_time))


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
