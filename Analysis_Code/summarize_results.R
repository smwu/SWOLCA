#=============================
## Plotting Simulation Results
## Programmer: SM Wu   
## Last Updated: 2023/04/02
#=============================

setwd("/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/")
library(R.matlab)
library(stringr)
library(abind)
library(gtools)
library(flextable)
library(dplyr)
library(bayesplot)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(gridExtra)
library(knitr)
library(kableExtra)
library(gt)

#============== Helper functions ===============================================
# get_true_params returns the observed simulated population values for the 
# parameters of interest
# Input: sim_pop: Simulated population
# Outputs: A list with the following elements:
#   true_pi: Kx1 vector of true pi values
#   true_theta: pxKxd array of true theta values
#   true_xi: (S+K)x1 vector of true xi values
#   true_Phi_mat: SxK matrix of true outcome proportions P(Y=1|S=s,C=k)
get_true_params <- function(sim_pop) {
  # Get true pi using population data
  true_pi <- tabulate(sim_pop$true_Ci) / length(sim_pop$true_Ci)
  # Get true theta
  theta_dim <- dim(sim_pop$true_global_thetas)
  true_theta <- array(NA, dim=theta_dim)
  for (j in 1:theta_dim[1]) {
    for (k in 1:theta_dim[2]) {
      for (r in 1:theta_dim[3]) {
        true_theta[j,k,r] <- sum((sim_pop$X_data[,j]==r) & 
                                 (sim_pop$true_Ci==k)) / sum(sim_pop$true_Ci==k) 
      }
    }
  }
  # Get true xi
  true_xi <- matrix(sim_pop$true_xi, nrow = theta_dim[2], byrow = FALSE)
  S <- length(unique(sim_pop$true_Si))
  K <- length(unique(sim_pop$true_Ci))
  # Get true Phi
  true_Phi_mat <- matrix(NA, nrow=K, ncol=S)
  for (k in 1:K) {
    for (s in 1:S) {
      true_Phi_mat[k, s] <- sum(sim_pop$Y_data==1 & sim_pop$true_Si==s & 
                                 sim_pop$true_Ci==k) / 
        sum(sim_pop$true_Si==s & sim_pop$true_Ci==k)
    }
  }
  # Get true xi
  true_xi <- qnorm(true_Phi_mat)
  
  return(list(true_pi = true_pi, true_theta = true_theta, true_xi = true_xi, 
              true_Phi_mat = true_Phi_mat))
}


# `get_theta_dist` returns squared Euclidean norm between estimated and true theta
# Inputs:
#   est_theta: Estimated theta 3D array. pxKxd
#   true_theta: True theta 3D array. pxKxd
#   order: Vector of optimal ordering of estimated theta elements. Kx1. Default is 
# NULL, in which case the best ordering is chosen over all possible permutations 
# Outputs: a list with the following elements:
#   theta_dist: Squared Euclidean norm between estimated and true thetas using optimal order
#   order: Optimal ordering to match estimated and true thetas. Kx1
#   est_theta_perm: Re-ordered estimated theta to match order of true theta. pxKxd
get_theta_dist <- function(est_theta, true_theta, order = NULL) {
  if (is.null(order)) {  # If no optimal ordering exists
    # Find all subsets of est_pi with size equal to true_pi
    all_perms <- gtools::permutations(n=dim(est_theta)[2], r=dim(true_theta)[2])
    # Obtain vector of dist (Frobenius norm) between est and true pi, calculated for each permutation
    dist_all_perms <- numeric(nrow(all_perms))
    for (i in 1:nrow(all_perms)) {
      est_theta_perm <- est_theta[,all_perms[i, ],]
      dist_all_perms[i] <- sum((est_theta_perm - true_theta)^2)
    }
    # Lowest dist out of all permutations
    theta_dist <- min(dist_all_perms)
    # Obtain optimal ordering of classes
    order <- all_perms[which.min(dist_all_perms), ]
    est_theta_perm <- est_theta[ , order, ]
  } else {  # Use optimal ordering if it exists
    est_theta_perm <- est_theta[ , order, ]
    theta_dist <- sum((est_theta_perm - true_theta)^2)
  }
  # Return dist, ordering, and reordered estimate
  return(list("theta_dist" = theta_dist, "order" = order, "est_theta_perm" = est_theta_perm))
}

# `get_pi_dist` returns squared Euclidean norm between estimated and true pi
# Inputs:
#   est_pi: Estimated pi vector. Kx1
#   true_pi: True pi vector. Kx1
#   order: Vector of optimal ordering of estimated pi elements. Kx1. Default is 
# NULL, in which case the best ordering is chosen over all possible permutations 
# Outputs: a list with the following elements:
#   pi_dist: Squared Euclidean norm between estimated and true pis using optimal order
#   order: Optimal ordering to match estimated and true pis. Kx1
#   est_pi_perm: Re-ordered estimated pi to match order of true pi. Kx1
get_pi_dist <- function(est_pi, true_pi, order=NULL) {
  if (is.null(order)) {  # If no optimal ordering exists
    # Find all subsets of est_pi with size equal to true_pi
    all_perms <- gtools::permutations(n=length(est_pi), r=length(true_pi))
    # Obtain vector of dist between est and true pi, calculated for each permutation
    dist_all_perms <- numeric(nrow(all_perms))
    for (i in 1:nrow(all_perms)) {
      est_pi_perm <- est_pi[all_perms[i, ]]
      dist_all_perms[i] <- sum((est_pi_perm - true_pi)^2)
    }
    # Lowest dist out of all permutations
    pi_dist <- min(dist_all_perms)
    # Obtain optimal ordering of classes
    order <- all_perms[which.min(dist_all_perms), ]
    est_pi_perm <- est_pi[order]
  } else {  # Use optimal ordering if it exists
    est_pi_perm <- est_pi[order]
    pi_dist <- sum((est_pi_perm - true_pi)^2)
  }
  # Return dist, ordering, and reordered estimate
  return(list("pi_dist" = pi_dist, "order" = order, "est_pi_perm" = est_pi_perm))
}

# `get_xi_dist` returns squared Euclidean norm between estimated and true xi
# Inputs:
#   est_xi: Estimated xi matrix. Kxq
#   true_xi: True xi matrix. Kxq
#   order: Vector of optimal ordering of estimated xi elements, determined in 
# get_theta_dist() or get_pi_dist()
#   Q: Number of covariates in probit model, excluding class assignment
# Outputs: a list with the following elements:
#   xi_dist: Squared Euclidean norm between estimated and true xis using optimal order
#   order: Optimal ordering to match estimated and true xis. Kx1
#   est_xi_perm: Re-ordered estimated xi to match order of true xi. Kxq
get_xi_dist <- function(est_xi, true_xi, order, Q) {
  est_xi_perm <- est_xi[order, ]
  xi_dist <- sum((est_xi_perm - true_xi)^2)
  # Return dist, ordering, and reordered estimate
  return(list("xi_dist" = xi_dist, "order" = order, "est_xi_perm" = est_xi_perm))
}

# `get_Phi_dist` returns squared Euclidean norm between estimated and true Phi,
# Inputs:
#   est_Phi: Estimated Phi vector; nx1
#   true_Phi: True Phi vector; nx1
# Outputs: a list with the following elements:
#   Phi_dist: Squared Euclidean norm between estimated and true Phis, calculated
# as a mean across all individual
get_Phi_dist <- function(est_Phi, true_Phi) {
  Phi_dist <- mean((est_Phi - true_Phi)^2)
  # Return dist
  return(list("Phi_dist" = Phi_dist))
}

# `get_marg_eff` returns the true marginal effect of P(Y=1|C)
# Inputs:
#   true_xi: Matrix of true probit model coefficients including S and C as 
# covariates. Applying 'pnorm' gives conditional effects P(Y=1|S,C). KxS
#   true_pi_s: Matrix of stratum-specific class probabilities P(C=k|S=s). SxK
#   true_pi: Vector of overall class membership probabilities P(C=k). Kx1
#   N_s: Vector of stratum sizes. Sx1
#   N: Population size
# Output: `p_y_cond_k` vector of marginal effect P(Y=1|C=k). Kx1
get_marg_eff <- function(true_xi, true_pi_s, true_pi, N_s, N) {
  S <- ncol(true_xi)
  K <- nrow(true_xi)
  p_s <- N_s / N
  p_y_cond_k <- numeric(K)
  p_y_s_cond_k <- matrix(NA, nrow=K, ncol=S)
  for (k in 1:K) {
    for (s in 1:S) {
      # P(Y=1|S=s, C=k)P(C=k|S=s)P(S=s)/P(C=k)
      p_y_s_cond_k[k, s] <- pnorm(true_xi[k, s]) * true_pi_s[s, k] * p_s[s] / 
        true_pi[k]
    }
    p_y_cond_k[k] <- sum(p_y_s_cond_k[k, ])
  }
  return(p_y_cond_k)
}

#============== Get performance metrics ========================================

# `get_metrics` obtains simulation metrics for posterior parameter estimates 
# across the sample iterations
# Inputs:
#   scen_pop: Three-digit population data scenario
#   scen_samp: Five-digit sample data scenario
#   iter_pop: Population data iteration
#   samp_n_seq: Vector sequence of sample indices
#   model: String indicating type of model, e.g., "wsOFMM", "sOFMM", "wOFMM"
#   coverage: Boolean indicating if coverage intervals should be calculated
#   posthoc: Boolean indicating if a posthoc pi is included in the results
#   plot_pi: Boolean indicating if output for pi plots is needed
#   plot_theta: Boolean indicating if output for theta plots is needed
#   plot_xi: Boolean indicating if output for xi plots is needed
# Outputs: list with the following components:
#   K_bias2: bias2 for number of clusters K
#   pi_bias2: bias2 for cluster probabilities pi
#   theta_bias2: bias2 for item response probabilities theta
#   xi_bias2: bias2 for probit coefficients xi
#   posthoc_pi_bias2: if posthoc = TRUE, bias2 for pi applying weights posthoc
#   pi_all: if plot_pi = TRUE, estimated pi values across samples
#   theta_mode: if plot_theta = TRUE, estimate modal theta values across samples
#   mode_mis: if plot_theta = TRUE, number of modal mismatches compared to true theta
#   xi_all: if plot_xi = TRUE, estimated xi values across samples
#   pi_cover_avg: if coverage = TRUE, coverage of pi estimate
#   theta_cover_avg: if coverage = TRUE, coverage of theta estimate
#   xi_cover_avg: if coverage = TRUE, coverage of xi estimate
#   runtime_avg: Average runtime across iterations
get_metrics <- function(wd, data_dir, res_dir, scen_pop, scen_samp, iter_pop=1, 
                        samp_n_seq, model, plot = TRUE) {
  
  #============== Load data and initialize variables ===========================
  # Load simulated population data
  pop_data_path <- paste0(wd, data_dir, "simdata_scen", scen_pop, "_iter",
                          iter_pop, ".RData")
  load(pop_data_path)
  # pop_data_path <- paste0(wd, data_dir, "simdata_scen", scen_pop, "_iter", 
  #                         iter_pop, ".mat") 
  # sim_pop <- readMat(pop_data_path)$sim.data
  # names(sim_pop) <- str_replace_all(dimnames(sim_pop)[[1]], "[.]", "_")
  
  
  # Obtain true observed population parameters
  true_params <- get_true_params(sim_pop = sim_pop)
  true_K <- as.vector(sim_pop$true_K)
  
  # Initialize variables
  L <- length(samp_n_seq)  # Number of samples
  runtime_all <- numeric(L)
  # Bias squared using posterior median
  K_dist <- pi_dist <- theta_dist <- xi_dist <- Phi_dist <- rep(NA, L) 
  # Posterior variance
  pi_var_all <- theta_var_all <- xi_var_all <- Phi_var_all <- rep(NA, L) 
  # Coverage variables
  pi_cover <- matrix(NA, nrow=L, ncol=length(true_params$true_pi))
  theta_cover <- array(NA, c(L, dim(true_params$true_theta)[c(1,2)]))
  xi_cover <- array(NA, c(L, dim(true_params$true_xi)))
  
  # Initialize plotting structures
  if (plot) {
    pi_all <- matrix(NA, nrow=L, ncol=true_K)
    theta_mode_all <- array(NA, dim=c(L, dim(sim_pop$true_global_patterns)))
    mode_mis_all <- rep(NA, L)
    xi_all <- array(NA, c(L, dim(true_params$true_xi)))
  }
  
  #============== Get performance metrics for each iteration ===================
  for (l in 1:L) { # For each sample iteration
    samp_n = samp_n_seq[l]
    print(samp_n)
    
    # Read in sample data. If file does not exist, move on to next iteration
    sim_data_path <- paste0(wd, data_dir, "simdata_scen", scen_samp, "_iter",
                             iter_pop, "_samp", samp_n, ".RData")
    if (!file.exists(sim_data_path)) {
      print(paste0("File does not exist: simdata_scen", scen_samp,"_iter",
                   iter_pop, "_samp", samp_n))
      next
    }
    load(sim_data_path)
    sim_samp <- sim_data
    
    # sim_data_path <- paste0(wd, data_dir, "simdata_scen", scen_samp, "_iter", 
    #                         iter_pop, "_samp", samp_n, ".mat") 
    # sim_samp <- readMat(sim_data_path)$sim.data
    # names(sim_samp) <- str_replace_all(dimnames(sim_samp)[[1]], "[.]", "_")
    
    
    # Read in results data
    if (model == "wsOFMM") {
      # wsOFMM model includes a variance adjustment
      sim_res_path <- paste0(wd, res_dir, model, "_results_2mod20000_adjRcpp_scen", scen_samp, 
                             "_samp", samp_n, ".RData")
      if (!file.exists(sim_res_path)) {
        print(paste0("File does not exist: ", model, "_results_2mod20000_adjRcpp_scen", 
                     scen_samp, "_samp", samp_n, ".RData"))
        next
      } 
      load(sim_res_path)
      analysis <- res$analysis_adj
      names(analysis) <- str_replace_all(names(analysis), "_adj", "")
      runtime_all[l] <- res$runtime
    } else {  # sOFMM and wOFMM models
      if(model == "sOFMM") {
        sim_res_path <- paste0(wd, res_dir, model, "_results_2mod20000_scen", scen_samp, 
                               "_samp", samp_n, ".RData")
        if (!file.exists(sim_res_path)) {
          print(paste0("File does not exist: ", model, "_results_2mod20000_scen", 
                       scen_samp, "_samp", samp_n, ".mat"))
          next
        } 
      } else {  # wOFMM 
        sim_res_path <- paste0(wd, res_dir, model, "_results_wt_2mod20000_scen", scen_samp, 
                               "_samp", samp_n, ".RData")
        if (!file.exists(sim_res_path)) {
          print(paste0("File does not exist: ", model, "_results_wt_2mod20000_scen", 
                       scen_samp, "_samp", samp_n, ".mat"))
          next
        } 
      }
      load(sim_res_path)
      analysis <- res$analysis
      runtime_all[l] <- res$runtime
    }
      
    S <- length(unique(sim_samp$true_Si))  # Number of strata
    M <- dim(analysis$theta_red)[1]        # Number of MCMC iterations
    p <- dim(analysis$theta_red)[2]        # Number of exposure items
    d <- dim(analysis$theta_red)[4]        # Number of exposure levels
    K <- length(analysis$pi_med)           # Number of classes
    Q <- dim(analysis$xi_med)[2]           # Number of additional covariates
    
    # For the effect modifier scenario, change true_xi to the marginal effect
    if (scen_pop == 1121) {
      true_params$true_xi <- qnorm(as.matrix(get_marg_eff(sim_pop$true_xi, 
                                                    sim_pop$true_pi_s, 
                                                    sim_pop$true_pi, 
                                                    sim_pop$N_s, sim_pop$N), 
                                       ncol = 1))
    }
      
    # If number of classes is incorrect, fill remaining components with 0's
    if (K > true_K) {
      # If there are extra estimated classes, add 0s to true parameters
      extra <- K - true_K
      true_params$true_pi <- c(true_params$true_pi, rep(0, extra))
      true_params$true_xi <- rbind(true_params$true_xi, 
                                   matrix(0, nrow = extra, ncol = Q))
      filler <- array(0, dim=c(dim(analysis$theta_med)[1], extra, 
                               dim(analysis$theta_med)[3]))
      true_params$true_theta <- abind(true_params$true_theta, filler, along = 2)
    } else if (K < true_K) {
      # If there are missing estimated classes, add 0s to estimated parameters
      missing <- true_K - K
      analysis$pi_med <- c(analysis$pi_med, rep(0, missing))
      analysis$xi_med <- rbind(analysis$xi_med, 
                               matrix(0, nrow = missing, ncol = Q))
      filler <- array(0, dim=c(dim(analysis$theta_med)[1], missing, 
                               dim(analysis$theta_med)[3]))
      analysis$theta_med <- abind(analysis$theta_med, filler, along = 2)  
      
      # Add 0's to full MCMC outputs for the missing classes
      analysis$pi_red <- abind(analysis$pi_red, array(0, dim=c(M, missing)), 
                               along=2)
      analysis$theta_red <- abind(analysis$theta_red, 
                                  array(0, dim=c(M, p, missing, d)), along = 3)
      # Special treatment of xi for two-step unsupervised model
      if (model == "wOFMM") {
        analysis$xi_med_lb <- rbind(matrix(analysis$xi_med_lb, nrow = K, ncol = Q), 
                                    matrix(0, nrow = missing, ncol = Q))
        analysis$xi_med_ub <- rbind(matrix(analysis$xi_med_ub, nrow = K, ncol = Q), 
                                    matrix(0, nrow = missing, ncol = Q))
      } else {
        analysis$xi_red <- abind(analysis$xi_red, 
                                 array(0, dim=c(M, missing, Q)), along = 2)
      }
    }
    
    #============== Calculated squared Euclidean norm (bias^2) =================
    ##### Number of classes, K
    K_dist[l] <- (K - true_K)^2
    
    ##### theta: get dist (Eucl norm) and optimal ordering
    theta_perm <- get_theta_dist(est_theta = analysis$theta_med, 
                                 true_theta = true_params$true_theta, 
                                 order=NULL)
    theta_dist[l] <- theta_perm$theta_dist
    order <- theta_perm$order
    # Theta mode consumption levels for each item and class (pxK)
    est_modes <- apply(analysis$theta_med[, order, ], c(1,2), which.max)
    true_modes <- apply(true_params$true_theta[,,], c(1,2), which.max)
    # True modal probabilities for each item and class (pxK)
    true_theta_modal <- apply(true_params$true_theta[,,], c(1,2), max)  
    
    ##### pi 
    pi_perm <- get_pi_dist(est_pi = analysis$pi_med, 
                           true_pi = true_params$true_pi, order = order)
    pi_dist[l] <- pi_perm$pi_dist
    
    ##### xi
    xi_perm <- get_xi_dist(est_xi = analysis$xi_med, 
                           true_xi = true_params$true_xi, order = order, Q = Q)
    xi_dist[l] <- xi_perm$xi_dist
    
    ##### Phi
    Phi_dist[l] <- get_Phi_dist(est_Phi = analysis$Phi_med,
                                     true_Phi = c(sim_samp$true_Phi))$Phi_dist
    
    #============== Calculate coverage and CI widths ===========================
    ##### pi
    # Obtain credible intervals for each of the K true clusters
    pi_CI <- apply(analysis$pi_red[, order], 2, 
                   function(x) quantile(x, c(0.025, 0.975)))
    # Assign 1 if interval covers true value, 0 if not
    pi_cover[l, ] <- ifelse((true_params$true_pi >= pi_CI[1,]) & 
                                   (true_params$true_pi <= pi_CI[2,]), 1, 0)
    # CI width averaged over the components
    pi_var_all[l] <- mean(apply(pi_CI, 2, diff))
    
    ##### theta
    theta_var_temp <- numeric(K)
    for (k in 1:K) {
      # Subset theta for cluster k
      est_theta_k <- analysis$theta_red[,,order[k],]
      # Each row provides the indices for one row of modal probabilities
      modal_idx <- cbind(rep(1:M, each=p), rep(1:p, times=M), 
                         rep(est_modes[, k], times=M))
      # estimated probabilities for the mode for cluster k (Mxp)
      est_theta_k_modal <- matrix(est_theta_k[modal_idx], ncol=p, byrow=TRUE)
      # Obtain credible intervals for each item 
      # Margins of apply are the dimensions that should be preserved
      theta_CI <- apply(est_theta_k_modal, 2, 
                        function(x) quantile(x, c(0.025, 0.975)))
      theta_cover[l, , k] <- ifelse((true_theta_modal[,k] >= theta_CI[1,]) &
                                  (true_theta_modal[,k] <= theta_CI[2,]), 1, 0)
      # CI width measures variation in estimating the modes for each k,
      # averaged over the items
      theta_var_temp[k] <- mean(apply(theta_CI, 2, diff))
    }
    # CI width averaged over the classes
    theta_var_all[l] <- mean(theta_var_temp)
    
    ##### xi
    if (model == "wOFMM") {
      # For two-step model, 'coefCI(fitglm)' is used in the Matlab code to 
      # extract CI from the regression model
      xi_CI <- array(NA, dim = c(2, dim(analysis$xi_med)))
      xi_CI[1, , ] <- analysis$xi_med_lb[order, ]
      xi_CI[2, , ] <- analysis$xi_med_ub[order, ]
    } else {
      # Obtain credible intervals for each component
      # Be careful not to drop size-1 dimension
      xi_CI <- apply(analysis$xi_red[, order, , drop = FALSE], c(2, 3), 
                     function(x) quantile(x, c(0.025, 0.975)))
    }
    # Assign 1 if interval covers true value, 0 if not
    xi_cover[l, , ] <- ifelse((true_params$true_xi >= xi_CI[1,,]) & 
                                   (true_params$true_xi <= xi_CI[2,,]), 1, 0)
    # CI width averaged over the components
    xi_var_all[l] <- mean(apply(xi_CI, c(2,3), diff))
    
    #============== Parameter estimate plot outputs ============================
    ##### theta
    if (plot) {
      if (K != true_K) { # mismatch of number of classes
        # Expand estimated array size if necessary
        if (dim(theta_mode_all)[3] < K) {
          filler <- array(NA, dim=c(L, dim(analysis$theta_med)[1], extra))
          theta_mode_all <- abind(theta_mode_all, filler, along = 3)
        }
      } 
      theta_mode_all[l,,1:true_K] <- est_modes
      # Mode mismatches
      mode_mis_all[l] <- sum(abs(est_modes - sim_pop$true_global_patterns))
    
      ##### pi
      if (K != true_K) { # mismatch of number of classes
        # Expand estimated matrix size if necessary
        if (dim(pi_all)[2] < K) {
          filler <- array(NA, dim=c(L, extra))
          pi_all <- abind(pi_all, filler, along = 2)
        }
      }
      pi_all[l,1:true_K] <- analysis$pi_med[order]
    
      ##### xi
      if (K != true_K) { # mismatch of number of classes
        # Expand estimated array size if necessary
        if (dim(xi_all)[2] < K) {
          filler <- array(NA, dim=c(L, extra, Q))
          xi_all <- abind(xi_all, filler, along = 2)
        }
      }
      xi_all[l, , ] <- analysis$xi_med[order, ] 
    }
  }

  #============== Calculate bias^2 averaged over sample iterations =============
  K_bias2 <- mean(K_dist, na.rm = TRUE)
  pi_bias2 <- mean(pi_dist, na.rm = TRUE)
  theta_bias2 <- mean(theta_dist, na.rm = TRUE)
  xi_bias2 <- mean(xi_dist, na.rm = TRUE)
  Phi_bias2 <- mean(Phi_dist, na.rm = TRUE)
  
  # Calculated CI width, averaged across iterations
  pi_var <- mean(pi_var_all, na.rm = TRUE)
  theta_var <- mean(theta_var_all, na.rm = TRUE)
  xi_var <- mean(xi_var_all, na.rm = TRUE)

  # Calculate class-specific coverage, averaged across iterations
  # Coverage for pi
  pi_cover_avg <- colMeans(pi_cover, na.rm = TRUE)
  # Coverage for theta: average over food items
  theta_cover_avg <- colMeans(colMeans(theta_cover, na.rm = TRUE), na.rm = TRUE)
  # Coverage for xi: average over additional Q covariates
  xi_cover_avg <- rowMeans(colMeans(xi_cover, na.rm = TRUE), na.rm = TRUE)
  
  runtime_avg <- mean(runtime_all, na.rm = TRUE)
  
  #============== Return results ===============================================
  ret_list <- list(K_bias2 = K_bias2, pi_bias2 = pi_bias2, pi_var = pi_var, 
                   theta_bias2 = theta_bias2, theta_var = theta_var, 
                   xi_bias2 = xi_bias2, xi_var = xi_var, 
                   Phi_bias2 = Phi_bias2, pi_cover_avg = pi_cover_avg,
                   theta_cover_avg = theta_cover_avg, xi_cover_avg = xi_cover_avg,
                   runtime_avg = runtime_avg, K_dist = K_dist, pi_dist = pi_dist,
                   theta_dist = theta_dist, xi_dist = xi_dist)
  
  if (plot) {
    ret_list[["pi_all"]] <- pi_all
    theta_mode <- apply(theta_mode_all, c(2,3), function(x) mean(x, na.rm = TRUE))
    mode_mis <- mean(mode_mis_all, na.rm = TRUE)
    ret_list[["theta_mode"]] <- theta_mode
    ret_list[["mode_mis"]] <- mode_mis
    ret_list[["xi_all"]] <- xi_all
    ret_list[["mode_mis_all"]] <- mode_mis_all
  }
  
  return(ret_list)
}

#================ SUMMARIES AND PLOTS ==========================================
wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/"
data_dir <- "Data/"
res_dir <- "Results/"
scen_pop <- 1211
scen_samp <- 121111
iter_pop <- 1
samp_n_seq <- 1:100
L <- length(samp_n_seq)
model <- "wsOFMM"
plot <- TRUE


#================ Baseline scenario ============================================
scenario <- "Baseline"
# Get all metrics and outputs
metrics_SRS_ws <- get_metrics(wd=wd, data_dir=data_dir, res_dir=res_dir, 
                              scen_pop=1, scen_samp=101, iter_pop=1, 
                              samp_n_seq=samp_n_seq, model="wsOFMM")
metrics_SRS_s <- get_metrics(wd=wd, data_dir=data_dir, res_dir=res_dir, 
                             scen_pop=1, scen_samp=101, iter_pop=1, 
                             samp_n_seq=samp_n_seq,  model="sOFMM")
metrics_SRS_unsup <- get_metrics(wd=wd, data_dir=data_dir, res_dir=res_dir, 
                                 scen_pop=1, scen_samp=101, iter_pop=1, 
                                 samp_n_seq=samp_n_seq, model="wOFMM")
metrics_Strat_ws <- get_metrics(wd=wd, data_dir=data_dir, res_dir=res_dir, 
                                scen_pop=scen_pop, scen_samp=scen_samp, iter_pop=1, 
                                samp_n_seq=samp_n_seq, model="wsOFMM")
metrics_Strat_s <- get_metrics(wd=wd, data_dir=data_dir, res_dir=res_dir, 
                               scen_pop=scen_pop, scen_samp=scen_samp, iter_pop=1, 
                               samp_n_seq=samp_n_seq, model="sOFMM")
metrics_Strat_unsup <- get_metrics(wd=wd, data_dir=data_dir, res_dir=res_dir, 
                                   scen_pop=scen_pop, scen_samp=scen_samp, iter_pop=1, 
                                   samp_n_seq=samp_n_seq, model="wOFMM")
# Load simulated population data
load(paste0(wd, data_dir, "simdata_scen", scen_pop,"_iter", iter_pop, ".RData"))

# Obtain true observed population parameters
true_params <- get_true_params(sim_pop = sim_pop)    

#================ Supervised scenario ==========================================
scenario <- "Supervised"
# Get all metrics and outputs
metrics_SRS_ws <- get_metrics(wd=wd, data_dir=data_dir, res_dir=res_dir, 
                              scen_pop=2, scen_samp=102, iter_pop=1, 
                              samp_n_seq=samp_n_seq, model="wsOFMM")
metrics_SRS_s <- get_metrics(wd=wd, data_dir=data_dir, res_dir=res_dir, 
                             scen_pop=2, scen_samp=102, iter_pop=1, 
                             samp_n_seq=samp_n_seq,  model="sOFMM")
metrics_SRS_unsup <- get_metrics(wd=wd, data_dir=data_dir, res_dir=res_dir, 
                                 scen_pop=2, scen_samp=102, iter_pop=1, 
                                 samp_n_seq=samp_n_seq, model="wOFMM")
metrics_Strat_ws <- get_metrics(wd=wd, data_dir=data_dir, res_dir=res_dir, 
                                scen_pop=2, scen_samp=202, iter_pop=1, 
                                samp_n_seq=samp_n_seq, model="wsOFMM")
metrics_Strat_s <- get_metrics(wd=wd, data_dir=data_dir, res_dir=res_dir, 
                               scen_pop=2, scen_samp=202, iter_pop=1, 
                               samp_n_seq=samp_n_seq, model="sOFMM")
metrics_Strat_unsup <- get_metrics(wd=wd, data_dir=data_dir, res_dir=res_dir, 
                                   scen_pop=2, scen_samp=202, iter_pop=1, 
                                   samp_n_seq=samp_n_seq, model="wOFMM")
scen_pop <- 2
# Load simulated population data
sim_pop <- readMat(paste0(data_dir, "simdata_scen", scen_pop,"_iter", 
                          iter_pop, ".mat"))$sim.data
names(sim_pop) <- str_replace_all(dimnames(sim_pop)[[1]], "[.]", "_")
# Obtain true observed population parameters
true_params <- get_true_params(sim_pop = sim_pop)   

#================ TABLE METRICS SUMMARY ========================================
### Create table of metrics with bias and variance
metrics_summ <- as.data.frame(matrix(NA, nrow=3, ncol=12))
colnames(metrics_summ) <- c("Sampling Scheme", "Model", 
                            "K Bias^2", "$\\pi$ Bias^2", "$\\pi$ CI width", 
                            "$\\theta$ Bias^2", "$\\theta$ CI width", 
                            "$\\xi$ Bias^2", "$\\xi$ CI width", 
                            "$\\pi$ Coverage","$\\theta$ Coverage", "$\\xi$ Coverage")
metrics_summ[, 1] <- rep("Stratified", 3)
metrics_summ[, 2] <- rep(c("Unwtd(sOFMM)", "Wtd(wsOFMM)", "Unsup(wOFMM)"), 1)  ## latent versions
output_inds <- 1:7
metrics_summ[1, -c(1,2)] <- c(metrics_Strat_s[output_inds], 
                              mean(metrics_Strat_s$pi_cover_avg), 
                              mean(metrics_Strat_s$theta_cover_avg),
                              mean(metrics_Strat_s$xi_cover_avg))
metrics_summ[2, -c(1,2)] <- c(metrics_Strat_ws[output_inds], 
                              mean(metrics_Strat_ws$pi_cover_avg), 
                              mean(metrics_Strat_ws$theta_cover_avg),
                              mean(metrics_Strat_ws$xi_cover_avg))
metrics_summ[3, -c(1,2)] <- c(metrics_Strat_unsup[output_inds], 
                              mean(metrics_Strat_unsup$pi_cover_avg), 
                              mean(metrics_Strat_unsup$theta_cover_avg),
                              mean(metrics_Strat_unsup$xi_cover_avg))
metrics_summ %>% 
  gt(caption = paste0(scen_samp, " scenario metrics of posterior parameter estimates, averaged over 100 samples")) %>%
  cols_label("$\\pi$ Bias^2" = "π Bias^2", "$\\theta$ Bias^2" = "θ Bias^2", 
             "$\\xi$ Bias^2" = "ξ Bias^2",
             "$\\pi$ CI width" = "π CI width", "$\\theta$ CI width" = "θ CI width", 
             "$\\xi$ CI width" = "ξ CI width", "$\\pi$ Coverage" = "π Coverage",
             "$\\theta$ Coverage" = "θ Coverage", "$\\xi$ Coverage" = "ξ Coverage") %>%
  fmt_number(
    columns = 3:12,
    decimals = 4)

### Create table of metrics with bias and variance (OLD VERSION)
metrics_summ <- as.data.frame(matrix(NA, nrow=6, ncol=12))
colnames(metrics_summ) <- c("Sampling Scheme", "Model", 
                        "K Bias^2", "$\\pi$ Bias^2", "$\\pi$ CI width", 
                        "$\\theta$ Bias^2", "$\\theta$ CI width", 
                        "$\\xi$ Bias^2", "$\\xi$ CI width", 
                        "$\\pi$ Coverage","$\\theta$ Coverage", "$\\xi$ Coverage")
metrics_summ[, 1] <- c(rep("SRS", 3), rep("Stratified", 3))
metrics_summ[, 2] <- rep(c("Unwtd(sOFMM)", "Wtd(wsOFMM)", "Unsup(wOFMM)"), 2)  ## latent versions
output_inds <- 1:7
metrics_summ[1, -c(1,2)] <- c(metrics_SRS_s[output_inds])
metrics_summ[2, -c(1,2)] <- c(metrics_SRS_ws[output_inds])
metrics_summ[3, -c(1,2)] <- c(metrics_SRS_unsup[output_inds])
metrics_summ[4, -c(1,2)] <- c(metrics_Strat_s[output_inds], 
                              mean(metrics_Strat_s$pi_cover_avg), 
                              mean(metrics_Strat_s$theta_cover_avg),
                              mean(metrics_Strat_s$xi_cover_avg))
metrics_summ[5, -c(1,2)] <- c(metrics_Strat_ws[output_inds], 
                              mean(metrics_Strat_ws$pi_cover_avg), 
                              mean(metrics_Strat_ws$theta_cover_avg),
                              mean(metrics_Strat_ws$xi_cover_avg))
metrics_summ[6, -c(1,2)] <- c(metrics_Strat_unsup[output_inds], 
                              mean(metrics_Strat_unsup$pi_cover_avg), 
                              mean(metrics_Strat_unsup$theta_cover_avg),
                              mean(metrics_Strat_unsup$xi_cover_avg))
metrics_summ %>% 
  gt(caption = paste0(scen_samp, " scenario metrics of posterior parameter estimates, averaged over 100 samples")) %>%
  cols_label("$\\pi$ Bias^2" = "π Bias^2", "$\\theta$ Bias^2" = "θ Bias^2", 
             "$\\xi$ Bias^2" = "ξ Bias^2",
             "$\\pi$ CI width" = "π CI width", "$\\theta$ CI width" = "θ CI width", 
             "$\\xi$ CI width" = "ξ CI width", "$\\pi$ Coverage" = "π Coverage",
             "$\\theta$ Coverage" = "θ Coverage", "$\\xi$ Coverage" = "ξ Coverage") %>%
  fmt_number(
    columns = 3:12,
    decimals = 4)
  # fmt_scientific(
  #   columns = 3:9,
  #   decimals = 3
  # ) 

#================ COVERAGE PLOTS ===============================================
plot_pi_cov <- data.frame(
  coverage = c(metrics_SRS_s$pi_cover_avg, metrics_SRS_ws$pi_cover_avg, 
               metrics_SRS_unsup$pi_cover_avg, metrics_Strat_s$pi_cover_avg, 
               metrics_Strat_ws$pi_cover_avg, metrics_Strat_unsup$pi_cover_avg),
  pi = rep(c("pi (k=1)", "pi (k=2)", "pi (k=3)"), times=6),
  model = rep(rep(c("SOLCA (Unweighted)", "WSOLCA (Weighted)", 
                    "WOLCA (Unsupervised)"), each=3), times=2),
  sampling = rep(c("SRS", "Stratified"), each=9)
)
pi_plot <- plot_pi_cov %>%  ggplot(aes(x=pi, y=coverage, color=model)) +
  theme_bw() +
  labs(x = expression(pi), y = "Coverage", color = "Model") +
  geom_point(size=2.5, position = position_dodge(0.3)) + 
  geom_hline(aes(yintercept=0.95)) +
  scale_y_continuous(limits=c(0,1), breaks=sort(c(seq(0, 1, length.out=5), 0.95))) + 
  facet_grid(~sampling) +
  ggtitle(paste0(scenario, " scenario coverage for π over 100 samples"))

plot_theta_cov <- data.frame(
  coverage = c(metrics_SRS_s$theta_cover_avg, metrics_SRS_ws$theta_cover_avg, 
               metrics_SRS_unsup$theta_cover_avg, metrics_Strat_s$theta_cover_avg, 
               metrics_Strat_ws$theta_cover_avg, metrics_Strat_unsup$theta_cover_avg),
  theta = rep(c("theta_1", "theta_2", "theta_3"), times=6),
  model = rep(rep(c("SOLCA (Unweighted)", "WSOLCA (Weighted)", 
                    "WOLCA (Unsupervised)"), each=3), times=2),
  sampling = rep(c("SRS", "Stratified"), each=9)
)
theta_plot <- plot_theta_cov %>%  ggplot(aes(x=theta, y=coverage, color=model)) +
  theme_bw() +
  labs(x = expression(theta), y = "Coverage", color = "Model") +
  geom_point(size=2.5, position = position_dodge(0.3)) + 
  geom_hline(aes(yintercept=0.95)) +
  scale_y_continuous(limits=c(0,1), breaks=sort(c(seq(0, 1, length.out=5), 0.95))) + 
  facet_grid(~sampling) +
  ggtitle(paste0(scenario, " scenario coverage for θ over 100 samples"))

plot_xi_cov <- data.frame(
  coverage = c(metrics_SRS_s$xi_cover_avg, metrics_SRS_ws$xi_cover_avg, 
               metrics_SRS_unsup$xi_cover_avg, metrics_Strat_s$xi_cover_avg, 
               metrics_Strat_ws$xi_cover_avg, metrics_Strat_unsup$xi_cover_avg),
  xi = rep(c("xi_1", "xi_2", "xi_3", "xi_4", "xi_5", "xi_6"), times=6),
  model = rep(rep(c("SOLCA (Unweighted)", "WSOLCA (Weighted)", 
                    "WOLCA (Unsupervised)"), each=6), times=2),
  sampling = rep(c("SRS", "Stratified"), each=18)
)
xi_plot <- plot_xi_cov %>%  ggplot(aes(x=xi, y=coverage, color=model)) +
  theme_bw() +
  theme(legend.position="top") +
  labs(x = expression(xi), y = "Coverage", color = "Model") +
  geom_point(size=2.5, position = position_dodge(0.3)) + 
  geom_hline(aes(yintercept=0.95)) +
  scale_y_continuous(limits=c(0,1), breaks=sort(c(seq(0, 1, length.out=5), 0.95))) + 
  facet_grid(~sampling) +
  ggtitle(paste0(scenario, " scenario coverage for ξ over 100 samples"))

ggarrange(pi_plot, theta_plot, nrow = 1, common.legend = TRUE, legend = "top")
xi_plot

#================ PARAMETER PLOTS ACROSS ITERATIONS ============================

#### Plot pi as grouped boxplot over iterations

sim_pop <- readMat(paste0(data_dir, "simdata_scen", scen_pop,"_iter", 
                          iter_pop, ".mat"))$sim.data
names(sim_pop) <- str_replace_all(dimnames(sim_pop)[[1]], "[.]", "_")

# Obtain true observed population parameters
true_params <- get_true_params(sim_pop)
L <- length(samp_n_seq)
pi_plot_data <- as.data.frame(rbind(metrics_Strat_s$pi_all,  
                                    metrics_Strat_ws$pi_all, 
                                    metrics_Strat_unsup$pi_all))
colnames(pi_plot_data) <- paste0("pi_", 1:ncol(pi_plot_data))
pi_plot_data$Model <- c(rep("sOFMM", times=L), rep("wsOFMM", times=L), 
                        rep("wOFMM", times=L))
pi_plot_data <- pi_plot_data %>% gather("pi_component", "value", -Model)
ggplot(pi_plot_data, aes(x=pi_component, y=value, fill=Model)) +
  theme_bw() +
  geom_boxplot() +
  geom_segment(mapping=aes(x=0.5, xend=1.5, y=true_params$true_pi[1], 
                           yend=true_params$true_pi[1]),color="forestgreen") +
  geom_segment(mapping=aes(x=1.5, xend=2.5, y=true_params$true_pi[2], 
                           yend=true_params$true_pi[2]),color="forestgreen") +
  geom_segment(mapping=aes(x=2.5, xend=3.5, y=true_params$true_pi[3], 
                           yend=true_params$true_pi[3]),color="forestgreen") +
  ggtitle(paste0(scenario, " scenario parameter estimation for π under \nstratified sampling with unequal probabilities, \ndisplayed across 100 samples"))


### Plot xi over 100 iterations

xi_plot_data <- as.data.frame(rbind(t(apply(metrics_Strat_s$xi_all, 1, c)),  
                                    t(apply(metrics_Strat_ws$xi_all, 1, c)), 
                                    t(apply(metrics_Strat_unsup$xi_all, 1, c))))
colnames(xi_plot_data) <- paste0("xi_", 1:ncol(xi_plot_data))
xi_plot_data$Model <- c(rep("sOFMM", times=L), rep("wsOFMM", times=L), 
                        rep("wOFMM", times=L))
xi_plot_data <- xi_plot_data %>% gather("xi_component", "value", -Model)
ggplot(xi_plot_data, aes(x=xi_component, y=value, fill=Model)) +
  theme_bw() +
  geom_boxplot() +
  geom_segment(mapping=aes(x=0.5, xend=1.5, y=true_params$true_xi[1,1], 
                           yend=true_params$true_xi[1,1]),color="forestgreen") +
  geom_segment(mapping=aes(x=1.5, xend=2.5, y=true_params$true_xi[2,1], 
                           yend=true_params$true_xi[2,1]),color="forestgreen") +
  geom_segment(mapping=aes(x=2.5, xend=3.5, y=true_params$true_xi[3,1], 
                           yend=true_params$true_xi[3,1]),color="forestgreen") +
  geom_segment(mapping=aes(x=3.5, xend=4.5, y=true_params$true_xi[1,2], 
                           yend=true_params$true_xi[1,2]),color="forestgreen") +
  geom_segment(mapping=aes(x=4.5, xend=5.5, y=true_params$true_xi[2,2], 
                           yend=true_params$true_xi[2,2]),color="forestgreen") +
  geom_segment(mapping=aes(x=5.5, xend=6.5, y=true_params$true_xi[3,2], 
                           yend=true_params$true_xi[3,2]),color="forestgreen") +
  ggtitle(paste0(scenario, " scenario parameter estimation for ξ under \nstratified sampling with unequal probabilities, \ndisplayed across 100 samples"))

### Plot Phi over 100 iterations

Phi_plot_data <- as.data.frame(rbind(pnorm(t(apply(metrics_Strat_s$xi_all, 1, c))),  
                                     pnorm(t(apply(metrics_Strat_ws$xi_all, 1, c))), 
                                     pnorm(t(apply(metrics_Strat_unsup$xi_all, 1, c)))))
colnames(Phi_plot_data) <- c("(S=1,C=1)", "(S=1,C=2)", "(S=1,C=3)", 
                             "(S=2,C=1)", "(S=2,C=2)", "(S=2,C=3)")
Phi_plot_data$Model <- c(rep("sOFMM", times=L), rep("wsOFMM", times=L), 
                         rep("wOFMM", times=L))
Phi_plot_data <- Phi_plot_data %>% gather("Phi_component", "value", -Model)
true_Phi <- true_params$true_Phi_mat
ggplot(Phi_plot_data, aes(x=Phi_component, y=value, fill=Model)) +
  theme_bw() +
  geom_boxplot() +
  ylab("P(Y=1|S,C)") +
  geom_segment(mapping=aes(x=0.5, xend=1.5, y=true_Phi[1,1], 
                           yend=true_Phi[1,1]),color="forestgreen") +
  geom_segment(mapping=aes(x=1.5, xend=2.5, y=true_Phi[2,1], 
                           yend=true_Phi[2,1]),color="forestgreen") +
  geom_segment(mapping=aes(x=2.5, xend=3.5, y=true_Phi[3,1], 
                           yend=true_Phi[3,1]),color="forestgreen") +
  geom_segment(mapping=aes(x=3.5, xend=4.5, y=true_Phi[1,2], 
                           yend=true_Phi[1,2]),color="forestgreen") +
  geom_segment(mapping=aes(x=4.5, xend=5.5, y=true_Phi[2,2], 
                           yend=true_Phi[2,2]),color="forestgreen") +
  geom_segment(mapping=aes(x=5.5, xend=6.5, y=true_Phi[3,2], 
                           yend=true_Phi[3,2]),color="forestgreen") +
  ggtitle(paste0(scenario, " scenario parameter estimation for P(Y=1|S,C) under \nstratified sampling with unequal probabilities, \ndisplayed across 100 samples"))


### Plot theta over 100 iteraitons

theta_mode_plot <- function(theta_plot_data, x_label) {
  p <- dim(theta_plot_data)[1]
  K <- dim(theta_plot_data)[2]
  Item <- factor(as.character(1:p), levels = as.character(p:1))
  theta_plot_data <- data.frame(theta_plot_data, Item)
  colnames(theta_plot_data) <- c(1:3, "Item")
  theta_plot_data <- theta_plot_data %>% gather("Class", "Level", 1:3) 
  patterns <- ggplot(theta_plot_data, aes(x=Class, y=Item, fill=Level)) + 
    theme_classic() +
    xlab(x_label) +
    geom_tile(color="gray") + 
    geom_text(aes(label = round(Level,2)), col="white", cex=2.5) +
    scale_fill_gradient(trans = "reverse")
  return(patterns)
}

p_true <- theta_mode_plot(sim_pop$true_global_patterns, "True Classes")
p_sOFMM <- theta_mode_plot(metrics_Strat_s$theta_mode, "sOFMM Classes")
p_wsOFMM <- theta_mode_plot(metrics_Strat_ws$theta_mode, "wsOFMM Classes")
p_wOFMM <- theta_mode_plot(metrics_Strat_unsup$theta_mode, "wOFMM Classes")
p_comb <- ggarrange(p_true, p_sOFMM + theme(axis.title.y = element_blank()), 
                    p_wsOFMM + theme(axis.title.y = element_blank()),
                    p_wOFMM + theme(axis.title.y = element_blank()), 
                    nrow = 1, common.legend = TRUE, legend = "right")
annotate_figure(p_comb, 
                top = text_grob(paste0(scen_samp, " scenario pattern elicitation under stratified sampling with unequal probabilities, \naveraged across 100 samples")))




#================ TROUBLESHOOTING NUMBER OF CLUSTERS, K ========================

wrong_K <- which(K_dist == 1)
model <- "wsOFMM"
for (i in 1:length(wrong_K)) {
  samp_n <- wrong_K[i]
  # wsOFMM model includes a variance adjustment
  sim_res_path <- paste0(wd, res_dir, model, "_results_adjRcpp_scen", scen_samp, 
                         "_samp", samp_n, ".RData")
  load(sim_res_path)
  analysis <- res$analysis_adj
  names(analysis) <- str_replace_all(names(analysis), "_adj", "")
  print(paste0("i: ", i, ", K: ", length(analysis$pi_med)))
  
  MCMC_path <- paste0(wd, res_dir, model, "_results_scen", scen_samp, 
                     "_samp", samp_n, ".RData")  # Output file
  load(MCMC_path)
  K_med <- res_MCMC$post_MCMC_out$K_med
  print(paste0("K_med: ", K_med))
}

samp_n <- 24
model <- "sOFMM"
MCMC_path <- paste0(wd, res_dir, model, "_results_scen", scen_samp, 
                    "_samp", samp_n, ".RData")  # Output file
load(MCMC_path)
# sOFMM
MCMC_out <- res$MCMC_out
# Get median number of classes with >= 5% of individuals, over all iterations
M <- dim(MCMC_out$pi_MCMC)[1]  # Number of stored MCMC iterations
K_med <- round(median(rowSums(MCMC_out$pi_MCMC >= 0.05)))
print(paste0("sOFMM K_fixed: ", dim(MCMC_out$pi_MCMC)[2]))
print(paste0("sOFMM K_med: ", K_med))

# Cluster individuals into reduced number of classes using agglomerative clustering
# Calculate pairwise distance matrix using Hamming distance: proportion of 
# iterations where two individuals have differing class assignments
distMat_s <- hamming.distance(t(MCMC_out$c_all_MCMC))
dendrogram <- hclust(as.dist(distMat_s), method = "complete") # Hierarchical clustering dendrogram
dend_plot <- as.dendrogram(dendrogram)
dend_plot %>% set("branches_k_color", k = 3) %>% plot()
red_c_all <- cutree(dendrogram, k = K_med)                  # Group individuals into K_med classes
table(red_c_all)

# Reduce and reorder parameter estimates using new classes
p <- 30
q <- 2
d <- 4
pi <- matrix(NA, nrow = M, ncol = K_med)
theta <- array(NA, dim = c(M, p, K_med, d))
xi <- array(NA, dim = c(M, K_med, q))
for (m in 1:M) {
  iter_order <- relabel_red_classes[m, ]
  pi_temp <- MCMC_out$pi_MCMC[m, iter_order]
  pi[m, ] <- pi_temp / sum(pi_temp)
  theta[m, , , ] <- MCMC_out$theta_MCMC[m, , iter_order, ]
  xi[m, , ] <- MCMC_out$xi_MCMC[m, iter_order, ]
}
post_MCMC_out <- list(K_med = K_med, pi = pi, theta = theta, xi = xi)
plot(pi[,1], type = "l", ylab = "Pi")
plot(pi[,2], type = "l", ylab = "Pi")
plot(pi[,3], type = "l", ylab = "Pi")
plot(pi[,4], type = "l", ylab = "Pi")

get_mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
# For each iteration, relabel new classes using the most common old class assignment
relabel_red_classes <- matrix(NA, nrow = M, ncol = K_med)   # Apply new classes to each iteration
for (k in 1:K_med) {
  relabel_red_classes[, k] <- apply(as.matrix(MCMC_out$c_all_MCMC[, red_c_all == k]), 
                                    1, get_mode)
}
# Posterior median estimate for theta across iterations
theta_med_temp <- apply(theta, c(2, 3, 4), median)
# Posterior modal exposure categories for each exposure item and reduced class
theta_modes <- apply(theta_med_temp, c(1, 2), which.max)
# Identify unique classes
unique_classes <- which(!duplicated(theta_modes, MARGIN = 2))
print(theta_modes)


# wsOFMM
model <- "wsOFMM"
# wsOFMM model includes a variance adjustment
sim_res_path <- paste0(wd, res_dir, model, "_results_adjRcpp_scen", scen_samp, 
                       "_samp", samp_n, ".RData")
load(sim_res_path)
MCMC_out <- res_MCMC$MCMC_out
post_MCMC_out <- res_MCMC$post_MCMC_out
print(paste0("wsOFMM K_fixed: ", dim(MCMC_out$pi_MCMC)[2]))
print(paste0("wsOFMM K_med: ", post_MCMC_out$K_med))

distMat_ws <- hamming.distance(t(MCMC_out$c_all_MCMC))
dendrogram <- hclust(as.dist(distMat_ws), method = "complete") # Hierarchical clustering dendrogram
dend_plot <- as.dendrogram(dendrogram)
dend_plot %>% set("branches_k_color", k = 3) %>% plot()
red_c_all <- cutree(dendrogram, k = K_med)                  # Group individuals into K_med classes
table(red_c_all)

plot(res_MCMC$post_MCMC_out$pi[,1], type = "l", ylab = "Pi")
plot(res_MCMC$post_MCMC_out$pi[,2], type = "l", ylab = "Pi")
plot(res_MCMC$post_MCMC_out$pi[,3], type = "l", ylab = "Pi")

# Posterior median estimate for theta across iterations
theta_med_temp <- apply(post_MCMC_out$theta, c(2, 3, 4), median)
# Posterior modal exposure categories for each exposure item and reduced class
theta_modes <- apply(theta_med_temp, c(1, 2), which.max)
# Identify unique classes
unique_classes <- which(!duplicated(theta_modes, MARGIN = 2))
print(theta_modes)

theta_iter <- MCMC_out$theta_MCMC[1000,,,]
theta_iter_modes <- apply(theta_iter, c(1, 2), which.max)
theta_iter_modes

leave <- c(24, 59, 80, 84, 94)
samp_n_seq <- setdiff(1:100, leave)

for (k in 1:3) {
  for (s in 1:2) {
    print(paste0("k: ", k, ", s: ", s))
    print(sum(sim_pop$Y_data==1 & sim_pop$true_Si==s & 
                                sim_pop$true_Ci==k))
    print(sum(sim_pop$true_Si==s & sim_pop$true_Ci==k))
    print("sample")
    print(sum(sim_samp$Y_data==1 & sim_samp$true_Si==s & 
                sim_samp$true_Ci==k))
    print(sum(sim_samp$true_Si==s & sim_samp$true_Ci==k))
  }
}
plot(MCMC_out$pi_MCMC[,1], type = "l", ylim=c(0.1, 0.7))
lines(MCMC_out$pi_MCMC[,2], col = "red")
lines(MCMC_out$pi_MCMC[,3], col = "blue")


#================ Marginal effect of C on Y ====================================
# Load simulated population data
pop_data_path <- paste0(wd, data_dir, "simdata_scen", scen_pop, "_iter", 
                        iter_pop, ".RData") 
load(pop_data_path)

as.matrix(get_marg_eff(sim_pop$true_xi, sim_pop$true_pi_s, sim_pop$true_pi, 
                       sim_pop$N_s, sim_pop$N), ncol = 1)

