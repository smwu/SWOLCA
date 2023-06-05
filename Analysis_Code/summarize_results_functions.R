#==================================================
# Helper functions for Plotting Simulation Results
# Programmer: SM Wu   
# Last Updated: 2023/05/20
#==================================================

#============== Get true population values =====================================

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

#============== Calculate distance =============================================

# `get_dist` returns the element-wise distance between 'par1' and 'par2' 
# according to the distance metric specified in 'dist_type'
# Inputs:
#   par1: First object. Can be vector, matrix, array
#   par2: Second object. Can be vector, matrix, array
#   dist_type: String specifying the distance type. Default is "mean_abs". Can
# also be "sum_sq" or "mean_sq"
# Output: Element-wise distance
get_dist <- function(par1, par2, dist_type = "mean_abs") {
  if (dist_type == "mean_abs") {  # Mean absolute error
    dist <- mean(abs(par1 - par2))
  } else if (dist_type == "sum_sq") {  # Frobenius norm / squared Euclidean norm
    dist <- sum((par1 - par2)^2)
  } else if (dist_type == "mean_sq") {  # MSE
    dist <- mean((par1 - par2)^2)
  } else {
    stop("Error: dist_type must be 'mean_abs', 'sum_sq', or 'mean_sq' ")
  }
  return(dist)
}

# `get_subset_dist` returns the mean absolute distance between parameter values
# after subsetting large_par to be the same dimension as small_par
# Inputs:
#   large_par: Parameter or estimate with larger number of classes
#   small_par: Parameter or estimate with smaller number of classes
#   param_name: Name of parameter. Either "theta", "pi", or "xi"
#   dist_type: String specifying the distance type. Must be "mean_abs", sum_sq", 
# or "mean_sq"
# Output: list containing: 
#   par_dist: Minimum mean absolute distance between the parameter and estimate 
# after subsetting
#   order_large: Optimal permutation ordering for large_par
get_subset_dist <- function(large_par, small_par, param_name, dist_type) {
  if (param_name == "theta") {
    large_K <- dim(large_par)[2]   ## Change to handle 0's
    sum(large_par[1, , 1] != 0)
    small_K <- dim(small_par)[2]
  } else if (param_name == "pi") {
    large_K <- length(large_par)
    small_K <- length(small_par)
  } else if (param_name == "xi") {
    large_K <- dim(large_par)[1] 
    small_K <- dim(small_par)[1]
  } else {
    stop("Error: 'param_name' must be either 'theta', 'pi', or 'xi'")
  }
  # Find all subsets of large_K with size equal to small_K
  sub_perms <- gtools::permutations(n = large_K, r = small_K)
  # Obtain dist (Frobenius norm) between large_par and small_par per permutation
  dist_sub_perms <- numeric(nrow(sub_perms))
  for (i in 1:nrow(sub_perms)) {
    if (param_name == "theta") {
      large_par_sub <- large_par[ , sub_perms[i, ], ]
    } else if (param_name == "pi") {
      large_par_sub <- large_par[sub_perms[i, ]]
    } else if (param_name == "xi") {
      large_par_sub <- large_par[sub_perms[i, ], ]
    }
    dist_sub_perms[i] <- get_dist(small_par, large_par_sub, "mean_abs")
  }
  # Lowest dist out of all permutations
  par_dist <- min(dist_sub_perms)
  # Ordering corresponding to lowest dist
  order_large <- sub_perms[which.min(dist_sub_perms), ]
  
  return(list(par_dist = par_dist, order_large = order_large))
}

# `get_theta_dist` returns the mean absolute distance between estimated and true 
# theta values, chosen over all permutations, with an option to subset the 
# number of classes
# Inputs:
#   est_theta: Estimated theta 3D array. pxKxd
#   true_theta: True theta 3D array. pxKxd
#   est_K: Estimated number of classes
#   true_K: True number of classes
#   subset: Boolean specifying whether or not to subset to the smaller number of 
# classes when calculating the distance. Default is TRUE, which applies less 
# penalty to having the wrong number of classes
#   dist_type: String specifying the distance type. Must be "mean_abs", sum_sq", 
# or "mean_sq"
# Outputs: a list with the following elements:
#   theta_dist: Mean absolute distance between estimated and true thetas using 
# optimal ordering
#   order: Optimal ordering to match estimated and true thetas. Kx1
#   est_theta_perm: Re-ordered estimated theta to match order of true theta. pxKxd
#   order_sub_est: Optimal subset ordering for estimated classes. (K_min)x1,
# where K_min is the minimum of est_K and true_K
#   order_sub_true: Optimal subset ordering for true classes. (K_min)x1
get_theta_dist <- function(est_theta, true_theta, est_K, true_K, subset = TRUE,
                           dist_type) {
  
  ### First get minimum distance using full vectors (with additional 0's)
  ### to get optimal ordering
  # Find all permutations of est_theta and true_theta with filler 0's
  all_perms <- gtools::permutations(n = dim(est_theta)[2], r = dim(true_theta)[2])
  # Obtain vector of mean absolute distance between est and true theta, 
  # calculated for each permutation
  dist_all_perms <- numeric(nrow(all_perms))
  for (i in 1:nrow(all_perms)) {
    est_theta_perm <- est_theta[,all_perms[i, ],]
    dist_all_perms[i] <- get_dist(est_theta_perm, true_theta, 
                                  dist_type = dist_type) 
  }
  # Obtain optimal ordering of classes
  order <- all_perms[which.min(dist_all_perms), ]
  est_theta_perm <- est_theta[ , order, ]
  
  # Initialize subset ordering of classes
  order_sub_est <- order_sub_true <- NULL
  
  ### Option to use this minimum distance, or calculate minimum distance after
  ### subsetting (subset == TRUE)
  if (!subset) {  # No subsetting
    # Lowest dist out of all permutations
    theta_dist <- min(dist_all_perms)
  } else {  # Calculate distance after subsetting
    if (est_K < true_K) {  # If missing a true class
      theta_sub <- get_subset_dist(large_par = true_theta[, 1:true_K, ], 
                                    small_par = est_theta[, 1:est_K, ], 
                                    param_name = "theta", dist_type = dist_type)
      theta_dist <- theta_sub$par_dist         # Distance
      order_sub_true <- theta_sub$order_large  # Subset order for true_theta
      order_sub_est <- 1:est_K                 # Subset order for est_theta
    } else if (true_K < est_K) {  # If extra class
      theta_sub <- get_subset_dist(large_par = est_theta[, 1:est_K, ], 
                                    small_par = true_theta[, 1:true_K, ],
                                    param_name = "theta", dist_type = dist_type)
      theta_dist <- theta_sub$par_dist
      order_sub_est <- theta_sub$order_large
      order_sub_true <- 1:true_K
    } else {  # true_K == est_K
      # Lowest dist out of all permutations
      theta_dist <- min(dist_all_perms)
      order_sub_est <- order
      order_sub_true <- 1:true_K
    }
  }
  
  ### Return dist, ordering, reordered estimate, and subsetted orderings
  return(list("theta_dist" = theta_dist, "order" = order, 
              "est_theta_perm" = est_theta_perm, 
              "order_sub_est" = order_sub_est,
              "order_sub_true" = order_sub_true))
}

# `get_pi_dist` returns the mean absolute distance between estimated and true pi
# values, with an option to subset the number of classes
# Inputs:
#   est_pi: Estimated pi vector. Kx1
#   true_pi: True pi vector. Kx1
#   order: Vector of optimal ordering of estimated pi elements. Kx1 
#   est_K: Estimated number of classes
#   true_K: True number of classes
#   subset: Boolean specifying whether or not to subset to the smaller number of 
# classes when calculating the distance. Default is TRUE, which applies less 
# penalty to having the wrong number of classes
#   order_sub_est: Optimal subset ordering for estimated classes. (K_min)x1,
# where K_min is the minimum of est_K and true_K
#   order_sub_true: Optimal subset ordering for true classes. (K_min)x1
#   dist_type: String specifying the distance type. Must be "mean_abs", sum_sq", 
# or "mean_sq"
# Outputs: a list with the following elements:
#   pi_dist: Mean absolute distance between est and true pis using optimal order
#   order: Optimal ordering to match estimated and true pis. Kx1
#   est_pi_perm: Re-ordered estimated pi to match order of true pi. Kx1
get_pi_dist <- function(est_pi, true_pi, order, est_K, true_K, subset = TRUE,
                        order_sub_est, order_sub_true, dist_type) {
  ### Use input optimal ordering
  ### Option to subset when calculating minimum distance (subset == TRUE)
  if (!subset) {  # No subsetting
    pi_dist <- get_dist(est_pi[order], true_pi, dist_type = dist_type)
  } else {  # Calculate distance after subsetting
    pi_dist <- get_dist(est_pi[order_sub_est], true_pi[order_sub_true], 
                        dist_type = dist_type)
  }
  
  ### Return dist, ordering, and reordered estimate
  return(list("pi_dist" = pi_dist, "order" = order, "est_pi_perm" = est_pi[order]))
}


# `get_xi_dist` returns the mean absolute distance between estimated and true xi
# values, with an option to subset the number of classes
# Inputs:
#   est_xi: Estimated xi matrix. Kxq
#   true_xi: True xi matrix. Kxq
#   order: Vector of optimal ordering of estimated xi elements. Kx1
#   est_K: Estimated number of classes
#   true_K: True number of classes
#   subset: Boolean specifying whether or not to subset to the smaller number of 
# classes when calculating the distance. Default is TRUE, which applies less 
# penalty to having the wrong number of classes
#   order_sub_est: Optimal subset ordering for estimated classes. (K_min)x1,
# where K_min is the minimum of est_K and true_K
#   order_sub_true: Optimal subset ordering for true classes. (K_min)x1
#   dist_type: String specifying the distance type. Must be "mean_abs", sum_sq", 
# or "mean_sq"
# Outputs: a list with the following elements:
#   xi_dist: Mean absolute distance between est and true xis using optimal order
#   order: Optimal ordering to match estimated and true xis. Kx1
#   est_xi_perm: Re-ordered estimated xi to match order of true xi. Kxq
get_xi_dist <- function(est_xi, true_xi, order, est_K, true_K, subset = TRUE,
                        order_sub_est, order_sub_true, dist_type) {
  
  ### Use input optimal ordering
  ### Option to subset when calculating minimum distance (subset == TRUE)
  if (!subset) {  # No subsetting
    xi_dist <- get_dist(est_xi[order], true_xi, dist_type = dist_type)
  } else {  # Calculate distance after subsetting
    xi_dist <- get_dist(est_xi[order_sub_est], true_xi[order_sub_true], 
                        dist_type = dist_type)
  }
  
  ### Return dist, ordering, and reordered estimate
  return(list("xi_dist" = xi_dist, "order" = order, "est_xi_perm" = est_xi[order]))
}

# `get_Phi_dist` returns mean absolute distance between estimated and true Phi
# Inputs:
#   est_Phi: Estimated Phi vector; nx1
#   true_Phi: True Phi vector; nx1
#   dist_type: String specifying the distance type. Must be "mean_abs", sum_sq", 
# or "mean_sq"
# Outputs: a list with the following elements:
#   Phi_dist: Mean absolute distance between estimated and true Phis, calculated
# as a mean across all individual
get_Phi_dist <- function(est_Phi, true_Phi, dist_type) {
  Phi_dist <- get_dist(est_Phi, true_Phi, dist_type = dist_type)
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
#   pop_data_path: String path for population dataset
#   sim_data_path: String with main part of sample dataset path, ending in "_samp"
#   sim_res_path: String with main part of sample results path, ending in "_samp"
#   model: String indicating model. Must be one of "wsOFMM", "sOFMM", or "wOFMM"
#   samp_n_seq: Vector of sample iterations to summarize over
#   plot: Boolean indicating if output for plots is needed. Default is TRUE
#   marg: Boolean indicating if marginal outcome probabilities should be used. 
# Default is FALSE
#   subset: Boolean specifying whether or not to subset to the smaller number of 
# classes when calculating the distance. Default is TRUE, which applies less 
# penalty to having the wrong number of classes
#   dist_type: String specifying the distance type. Must be "mean_abs", sum_sq", 
# or "mean_sq". Default is "mean_abs"
# Outputs: list with the following components:
#   K_bias: bias for number of clusters K
#   pi_bias: bias for cluster probabilities pi
#   theta_bias: bias for item response probabilities theta
#   xi_bias: bias for probit coefficients xi
#   posthoc_pi_bias: if posthoc = TRUE, bias for pi applying weights posthoc
#   pi_all: if plot_pi = TRUE, estimated pi values across samples
#   theta_mode: if plot_theta = TRUE, estimate modal theta values across samples
#   mode_mis: if plot_theta = TRUE, number of modal mismatches compared to true theta
#   xi_all: if plot_xi = TRUE, estimated xi values across samples
#   pi_cover_avg: if coverage = TRUE, coverage of pi estimate
#   theta_cover_avg: if coverage = TRUE, coverage of theta estimate
#   xi_cover_avg: if coverage = TRUE, coverage of xi estimate
#   runtime_avg: Average runtime across iterations
get_metrics <- function(pop_data_path, sim_data_path, sim_res_path, model,
                        samp_n_seq, plot = TRUE, marg = FALSE, 
                        subset = TRUE, dist_type = "mean_abs") {
  
  #============== Load data and initialize variables ===========================
  # Load simulated population data
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
  pi_cover <- matrix(0, nrow=L, ncol=length(true_params$true_pi))
  theta_cover <- array(0, c(L, dim(true_params$true_theta)[c(1,2)]))
  xi_cover <- array(0, c(L, dim(true_params$true_xi)))
  # MSE for all iterations
  pi_mse_all <- theta_mse_all <- xi_mse_all <- rep(NA, L)
  
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
    
    # Obtain true observed population parameters
    # Need to reload so that filler dimensions do not keep adding over iterations
    true_params <- get_true_params(sim_pop = sim_pop)
    true_K <- as.vector(sim_pop$true_K)
    
    # Read in sample data. If file does not exist, move on to next iteration
    sim_data_path_full <- paste0(sim_data_path, samp_n, ".RData")
    if (!file.exists(sim_data_path_full)) {
      print(paste0("File does not exist: ", sim_data_path_full))
      next
    }
    load(sim_data_path_full)
    sim_samp <- sim_data
          # sim_data_path <- paste0(wd, data_dir, "simdata_scen", scen_samp, "_iter", 
          #                         iter_pop, "_samp", samp_n, ".mat") 
          # sim_samp <- readMat(sim_data_path)$sim.data
          # names(sim_samp) <- str_replace_all(dimnames(sim_samp)[[1]], "[.]", "_")
    
    # Read in results data
    sim_res_path_full <- paste0(sim_res_path, samp_n, ".RData")
    if (!file.exists(sim_res_path_full)) {
      print(paste0("File does not exist: ", sim_res_path_full))
      next
    } 
    if (model == "wsOFMM") {
      # wsOFMM model includes a variance adjustment
      load(sim_res_path_full)
      analysis <- res$analysis_adj
      names(analysis) <- str_replace_all(names(analysis), "_adj", "")
      runtime_all[l] <- res$runtime
    } else if (model == "sOFMM") { 
      load(sim_res_path_full)
      analysis <- res$analysis
      runtime_all[l] <- res$runtime
    } else if (model == "wOFMM") {
      load(sim_res_path_full)
      analysis <- res$analysis
      runtime_all[l] <- res$runtime
      n <- length(sim_samp$true_Si)
      ##### CHECK THIS!!!!!!!!!!!!!!
      Phi_med <- numeric(n)                    # Initialize individual outcome probabilities
      # Calculate posterior class membership, p(c_i=k|-), for each class k
      for (i in 1:n) {
        if (marg) {
          Phi_med[i] <- analysis$xi_med[analysis$c_all[i], ]
        }
        Phi_med[i] <- analysis$xi_med[analysis$c_all[i], sim_data$true_Si[i]]
      }
      analysis$Phi_med <- Phi_med
    } else {
      stop("Error: model must be one of 'wsOFMM', 'sOFMM', or 'wOFMM'")
    }
    
    S <- length(unique(sim_samp$true_Si))  # Number of strata
    M <- dim(analysis$theta_red)[1]        # Number of MCMC iterations
    p <- dim(analysis$theta_red)[2]        # Number of exposure items
    d <- dim(analysis$theta_red)[4]        # Number of exposure levels
    K <- length(analysis$pi_med)           # Number of classes
    Q <- dim(analysis$xi_med)[2]           # Number of additional covariates
    
    # For the effect modifier scenario, change true_xi to the marginal effect
    if (marg) {
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
    
    #============== Calculated mean absolute distance (abs bias) ===============
    ##### Number of classes, K
    K_dist[l] <- get_dist(K, true_K, dist_type = dist_type)
    
    ##### theta: get dist (Eucl norm) and optimal ordering
    theta_perm <- get_theta_dist(est_theta = analysis$theta_med, 
                                 true_theta = true_params$true_theta, 
                                 est_K = K, true_K = true_K, subset = subset,
                                 dist_type = dist_type)
    theta_dist[l] <- theta_perm$theta_dist
    order <- theta_perm$order
    if (subset) {
      order_sub_est <- theta_perm$order_sub_est
      order_sub_true <- theta_perm$order_sub_true
      K_min <- length(order_sub_est)
    } else {
      order_sub_est <- order
      order_sub_true <- 1:length(order)
      K_min <- K
    }
    
    ##### pi 
    pi_perm <- get_pi_dist(est_pi = analysis$pi_med, 
                           true_pi = true_params$true_pi, order = order, 
                           est_K = K, true_K = true_K, subset = subset,
                           order_sub_est = order_sub_est, 
                           order_sub_true = order_sub_true,
                           dist_type = dist_type)
    pi_dist[l] <- pi_perm$pi_dist
    
    ##### xi
    xi_perm <- get_xi_dist(est_xi = analysis$xi_med, 
                           true_xi = true_params$true_xi, order = order, 
                           est_K = K, true_K = true_K, subset = subset,
                           order_sub_est = order_sub_est, 
                           order_sub_true = order_sub_true,
                           dist_type = dist_type)
    xi_dist[l] <- xi_perm$xi_dist
    
    ##### Phi
    Phi_dist[l] <- get_Phi_dist(est_Phi = analysis$Phi_med,
                                true_Phi = c(sim_samp$true_Phi),
                                dist_type = "mean_sq")$Phi_dist
    
    #============== Calculate coverage and CI widths ===========================
    ##### pi
    # Obtain credible intervals for each of the K true clusters
    pi_CI <- apply(analysis$pi_red[, order_sub_est], 2, 
                   function(x) quantile(x, c(0.025, 0.975)))
    # Assign 1 if interval covers true value, 0 if not
    # If a class is missing, defaults to 0 (not covered)
    pi_cover[l, order_sub_true] <- ifelse(
      (true_params$true_pi[order_sub_true] >= pi_CI[1,]) & 
      (true_params$true_pi[order_sub_true] <= pi_CI[2,]), 1, 0)
    # CI width averaged over the components
    pi_var_all[l] <- mean(apply(pi_CI, 2, diff))
    # MSE
    pi_mse_all[l] <- mean(apply(analysis$pi_red, 1, function(x) 
      get_dist(x[order_sub_est], true_params$true_pi[order_sub_true], 
               "mean_sq")))
    
    
    ##### theta
    # Theta mode consumption levels for each item and class (pxK)
    est_modes <- apply(analysis$theta_med[, order_sub_est, ], c(1,2), which.max)
    true_modes <- apply(true_params$true_theta[, order_sub_true, ], c(1,2), 
                        which.max)
    # True modal probabilities for each item and class (pxK)
    true_theta_modal <- apply(true_params$true_theta[ , order_sub_true, ], 
                              c(1,2), max) 
    theta_var_temp <- numeric(K_min)
    for (k in 1:K_min) {
      # Subset theta for cluster k
      est_theta_k <- analysis$theta_red[, , order_sub_est[k], ]
      # Each row provides the indices for one row of modal probabilities
      modal_idx <- cbind(rep(1:M, each = p), rep(1:p, times = M), 
                         rep(est_modes[, k], times = M))
      # estimated probabilities for the mode for cluster k (Mxp)
      est_theta_k_modal <- matrix(est_theta_k[modal_idx], ncol = p, byrow = TRUE)
      # Obtain credible intervals for each item 
      # Margins of apply are the dimensions that should be preserved
      theta_CI <- apply(est_theta_k_modal, 2, 
                        function(x) quantile(x, c(0.025, 0.975)))
      theta_cover[l, , order_sub_true[k]] <- ifelse(
        (true_theta_modal[, k] >= theta_CI[1, ]) &
        (true_theta_modal[, k] <= theta_CI[2, ]), 1, 0)
      # CI width measures variation in estimating the modes for each k,
      # averaged over the items
      theta_var_temp[k] <- mean(apply(theta_CI, 2, diff))
    }
    # CI width averaged over the classes
    theta_var_all[l] <- mean(theta_var_temp)
    # MSE
    theta_mse_all[l] <- mean(apply(analysis$theta_red, 1, function(x) 
      get_dist(x[, order_sub_est, ], true_params$true_theta[, order_sub_true, ], 
               "mean_sq")))
    
    ##### xi
    if (model == "wOFMM") {
      # For two-step model, 'coefCI(fitglm)' is used in the Matlab code to 
      # extract CI from the regression model
      xi_CI <- array(NA, dim = c(2, dim(analysis$xi_med[order_sub_est, , 
                                                        drop = FALSE])))
      xi_CI[1, , ] <- analysis$xi_med_lb[order_sub_est, , drop = FALSE]
      xi_CI[2, , ] <- analysis$xi_med_ub[order_sub_est, , drop = FALSE]
      # Get MSE as bias^2 + var, recovering var of the mean from CI and s.e.
      # Then take mean over components of xi
      xi_se <- apply(xi_CI, c(2,3), function(x) (x[2] - x[1])/2 / qnorm(0.975))
      xi_mse_all[l] <- mean( (analysis$xi_med[order_sub_est, ] - 
                      true_params$true_xi[order_sub_true, ])^2 + xi_se^2 )
    } else {
      # Obtain credible intervals for each component
      # Be careful not to drop size-1 dimension
      xi_CI <- apply(analysis$xi_red[, order_sub_est, , drop = FALSE], c(2, 3), 
                     function(x) quantile(x, c(0.025, 0.975)))
      # MSE
      xi_mse_all[l] <- mean(apply(analysis$xi_red, 1, function(x) 
        get_dist(x[order_sub_est, ], true_params$true_xi[order_sub_true, ], 
                 "mean_sq")))
    }
    # Assign 1 if interval covers true value, 0 if not
    xi_cover[l, order_sub_true, ] <- ifelse(
      (true_params$true_xi[order_sub_true, ] >= xi_CI[1, , ]) & 
      (true_params$true_xi[order_sub_true, ] <= xi_CI[2, , ]), 1, 0)
    # CI width averaged over the components
    xi_var_all[l] <- mean(apply(xi_CI, c(2, 3), diff))
    
    
    #============== Parameter estimate plot outputs ============================
    ##### theta
    if (plot) {
      # Theta mode consumption levels for each item and class (pxK)
      est_modes <- apply(analysis$theta_med[, order, ], c(1,2), which.max)
      true_modes <- apply(true_params$true_theta[, order, ], c(1,2), 
                          which.max)
      if (K != true_K) { # mismatch of number of classes
        # Expand estimated array size if necessary
        if (dim(theta_mode_all)[3] < K) {
          filler <- array(NA, dim=c(L, dim(analysis$theta_med)[1], extra))
          theta_mode_all <- abind(theta_mode_all, filler, along = 3)
        }
      } 
      theta_mode_all[l, , 1:length(order)] <- est_modes
      # Mode mismatches
      mode_mis_all[l] <- sum(abs(est_modes[, 1:true_K] - 
                                   sim_pop$true_global_patterns))
      
      ##### pi
      if (K != true_K) { # mismatch of number of classes
        # Expand estimated matrix size if necessary
        if (dim(pi_all)[2] < K) {
          filler <- array(NA, dim=c(L, extra))
          pi_all <- abind(pi_all, filler, along = 2)
        }
      }
      pi_all[l, 1:length(order)] <- analysis$pi_med[order]
      
      ##### xi
      if (K != true_K) { # mismatch of number of classes
        # Expand estimated array size if necessary
        if (dim(xi_all)[2] < K) {
          filler <- array(NA, dim=c(L, extra, Q))
          xi_all <- abind(xi_all, filler, along = 2)
        }
      }
      xi_all[l, 1:length(order), ] <- analysis$xi_med[order, ] 
    }
  }
  
  #============== Calculate bias^2 averaged over sample iterations =============
  K_bias <- mean(K_dist, na.rm = TRUE)
  pi_bias <- mean(pi_dist, na.rm = TRUE)
  theta_bias <- mean(theta_dist, na.rm = TRUE)
  xi_bias <- mean(xi_dist, na.rm = TRUE)
  Phi_bias <- mean(Phi_dist, na.rm = TRUE)
  
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
  ret_list <- list(K_bias = K_bias, pi_bias = pi_bias, pi_var = pi_var, 
                   theta_bias = theta_bias, theta_var = theta_var, 
                   xi_bias = xi_bias, xi_var = xi_var, 
                   Phi_bias = Phi_bias, pi_cover_avg = pi_cover_avg,
                   theta_cover_avg = theta_cover_avg, xi_cover_avg = xi_cover_avg,
                   runtime_avg = runtime_avg, K_dist = K_dist, pi_dist = pi_dist,
                   theta_dist = theta_dist, xi_dist = xi_dist, 
                   pi_mse_all = pi_mse_all, theta_mse_all = theta_mse_all,
                   xi_mse_all = xi_mse_all)
  
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


# `save_scen_metrics` obtains and saves summary metrics for each scenario
# Inputs:
#   scen_pop: Four-digit population scenario. First digit specifies strength of 
# pattern, second digit specifies diet-outcome relationship, third digit 
# specifies selection (stratum) variable type, fourth digit specifies clustering
#   scen_samp: Six-digit sample scenario. First digit specifies strength of 
# pattern, second digit specifies diet-outcome relationship, third digit 
# specifies selection (stratum) variable type, fourth digit specifies clustering,
# fifth digit specifies sample size, sixth digit specifies sampling design
#   WSOLCA: Boolean indicating if WSOLCA results should be obtained
#   SOLCA: Boolean indicating if SOLCA results should be obtained
#   WOLCA: Boolean indicating if WOLCA results should be obtained
#   WSOLCA_name: String with file name identifier for WSOLCA results
#   SOLCA_name: String with file name identifier for SOLCA results
#   WOLCA_name: String with file name identifier for WOLCA results
#   marg: Boolean indicating if marginal outcome probabilities should be used. 
# Default is FALSE
#   save_name: String with main file name for saving all metrics
# Outputs: Saves 'metrics_all' list with all metrics in an .RData file
save_scen_metrics <- function(scen_pop, scen_samp, WSOLCA = TRUE, SOLCA = TRUE,
                              WOLCA = TRUE, WSOLCA_name, SOLCA_name, WOLCA_name,
                              marg = FALSE, save_name) {
  # Set parameters and paths
  wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/"  # Working directory
  data_dir <- "Data/"               # Simulated data directory
  res_dir <- "Results/"             # Model results directory
  analysis_dir <- "Analysis_Code/"  # Analysis directory where metrics will be saved
  iter_pop <- 1                     # Population iteration
  samp_n_seq <- 1:100               # Vector of sample iterations to summarize over
  plot <- TRUE                      # Whether output for plots is needed
  subset <- TRUE                    # Whether to subset to classes when calculating distance
  dist_type <- "mean_abs"           # Bias distance type
  # Population simulated data path
  pop_data_path <- paste0(wd, data_dir, "simdata_scen", scen_pop, "_iter",
                          iter_pop, ".RData")
  # Main portion of sample simulated data path
  sim_data_path <- paste0(wd, data_dir, "simdata_scen", scen_samp, "_iter",
                          iter_pop, "_samp")
  
  # Get metrics for models
  metrics_all <- list()
  if (WSOLCA) {
    print("Getting WSOLCA results...")
    model <- "wsOFMM"
    # Main portion of model results path
    sim_res_path <- paste0(wd, res_dir, model, WSOLCA_name, scen_samp, "_samp")
    metrics_ws <- get_metrics(pop_data_path = pop_data_path,
                              sim_data_path = sim_data_path, 
                              sim_res_path = sim_res_path, model = model,
                              samp_n_seq = samp_n_seq, plot = plot, marg = marg, 
                              subset = subset, dist_type = dist_type)
    metrics_all$metrics_ws <- metrics_ws
  } 
  if (SOLCA) {
    print("Getting SOLCA results...")
    model <- "sOFMM"
    # Main portion of model results path
    sim_res_path <- paste0(wd, res_dir, model, SOLCA_name, scen_samp, "_samp")
    metrics_s <- get_metrics(pop_data_path = pop_data_path,
                             sim_data_path = sim_data_path,
                             sim_res_path = sim_res_path, model = model,
                             samp_n_seq = samp_n_seq, plot = plot, marg = marg, 
                             subset = subset, dist_type = dist_type)
    metrics_all$metrics_s <- metrics_s
  }
  if (WOLCA) {
    print("Getting WOLCA results...")
    model <- "wOFMM"
    # Main portion of model results path
    sim_res_path <- paste0(wd, res_dir, model, WOLCA_name, scen_samp, "_samp")
    metrics_unsup <- get_metrics(pop_data_path = pop_data_path,
                                 sim_data_path = sim_data_path,
                                 sim_res_path = sim_res_path, model = model,
                                 samp_n_seq = samp_n_seq, plot = plot, marg = marg, 
                                 subset = subset, dist_type = dist_type)
    metrics_all$metrics_unsup <- metrics_unsup
  }
  
  # Save summary metrics
  save(metrics_all, 
       file = paste0(wd, analysis_dir, save_name, scen_samp, ".RData"))
  
}

#============== Create Table 1 =================================================

# `create_table1` uses 'gt()' to create a table of summary metrics comparing
# sampling schemes
# Inputs:
#   wd: String specifying working directory
#   analysis_dir: String specifying analysis directory where summaries are saved
#   format: String specifying kable table format. Must be "latex", "html", 
# "pipe", "simple", or "rst". Default is "latex"
# Output: Formatted table with absolute bias, CI width, and coverage
create_table1 <- function(wd, analysis_dir, format) {
  metrics_summ <- as.data.frame(matrix(NA, nrow = 9, ncol = 12))
  colnames(metrics_summ) <- c("Sampling Scheme", "Model", 
                              "$K$ Abs Bias", "$\\pi$ Abs Bias",  
                              "$\\theta$ Abs Bias", "$\\xi$ Abs Bias", 
                              "$\\pi$ CI Width", "$\\theta$ CI Width", "$\\xi$ CI Width", 
                              "$\\pi$ Coverage","$\\theta$ Coverage", "$\\xi$ Coverage")
  metrics_summ[, 1] <- c(rep("SRS", 3), rep("Stratified", 3), rep("Stratified Cluster", 3))
  metrics_summ[, 2] <- rep(c("SOLCA", "WOLCA", "WSOLCA"), 3)  
  # output_inds <- 1:7
  output_inds <- c(1, 2, 4, 6, 3, 5, 7)
  # SRS results
  load(paste0(wd, analysis_dir, "metrics_scen", 111213, ".RData"))
  metrics_summ[1, -c(1,2)] <- c(metrics_all$metrics_s[output_inds], 
                                mean(metrics_all$metrics_s$pi_cover_avg), 
                                mean(metrics_all$metrics_s$theta_cover_avg),
                                mean(metrics_all$metrics_s$xi_cover_avg))
  metrics_summ[2, -c(1,2)] <- c(metrics_all$metrics_unsup[output_inds], 
                                mean(metrics_all$metrics_unsup$pi_cover_avg), 
                                mean(metrics_all$metrics_unsup$theta_cover_avg),
                                mean(metrics_all$metrics_unsup$xi_cover_avg))
  metrics_summ[3, -c(1,2)] <- c(metrics_all$metrics_ws[output_inds], 
                                mean(metrics_all$metrics_ws$pi_cover_avg), 
                                mean(metrics_all$metrics_ws$theta_cover_avg),
                                mean(metrics_all$metrics_ws$xi_cover_avg))

  # Stratified results
  load(paste0(wd, analysis_dir, "metrics_scen", 111211, ".RData"))
  metrics_summ[4, -c(1,2)] <- c(metrics_all$metrics_s[output_inds], 
                                mean(metrics_all$metrics_s$pi_cover_avg), 
                                mean(metrics_all$metrics_s$theta_cover_avg),
                                mean(metrics_all$metrics_s$xi_cover_avg))
  metrics_summ[5, -c(1,2)] <- c(metrics_all$metrics_unsup[output_inds], 
                                mean(metrics_all$metrics_unsup$pi_cover_avg), 
                                mean(metrics_all$metrics_unsup$theta_cover_avg),
                                mean(metrics_all$metrics_unsup$xi_cover_avg))
  metrics_summ[6, -c(1,2)] <- c(metrics_all$metrics_ws[output_inds], 
                                mean(metrics_all$metrics_ws$pi_cover_avg), 
                                mean(metrics_all$metrics_ws$theta_cover_avg),
                                mean(metrics_all$metrics_ws$xi_cover_avg))
  
  # Stratified cluster results
  load(paste0(wd, analysis_dir, "metrics_scen", 111212, ".RData"))
  metrics_summ[7, -c(1,2)] <- c(metrics_all$metrics_s[output_inds], 
                                mean(metrics_all$metrics_s$pi_cover_avg), 
                                mean(metrics_all$metrics_s$theta_cover_avg),
                                mean(metrics_all$metrics_s$xi_cover_avg))
  metrics_summ[8, -c(1,2)] <- c(metrics_all$metrics_unsup[output_inds], 
                                mean(metrics_all$metrics_unsup$pi_cover_avg), 
                                mean(metrics_all$metrics_unsup$theta_cover_avg),
                                mean(metrics_all$metrics_unsup$xi_cover_avg))
  metrics_summ[9, -c(1,2)] <- c(metrics_all$metrics_ws[output_inds], 
                                mean(metrics_all$metrics_ws$pi_cover_avg), 
                                mean(metrics_all$metrics_ws$theta_cover_avg),
                                mean(metrics_all$metrics_ws$xi_cover_avg))
  
  metrics_summ %>% 
    kbl(digits = 4, align = "rrrrrrrrrrrr", booktabs = TRUE, format = format,
        caption = "Summary of mean absolute bias, 95% credible interval width, and coverage for simulations based on posterior samples.") %>%
    kable_classic() %>%
    kable_styling(full_width = FALSE)
}


#================ Create Figure 1 ==============================================
# `plot_mse_boxplot` plots MSE as a grouped boxplot
# Inputs:
#   wd: String specifying working directory
#   analysis_dir: String specifying analysis directory where summaries are saved
#   save_names: String vector with main file names for saving all metrics
#   scenarios: Numeric vector specifying sampling scenarios to include
#   scen_names: String vector specifying labels for sampling scenarios
#   overall_name: String specifying overall grouping name for sampling scenarios
#   param: String specifying parameter to plot. Must be "pi", "theta", or "xi"
#   lower_lim: Lower limit for y-axis. Default is 0
#   upper_lim: Upper limit for y_axis. Default is 1
#   xlab: x-axis label
#   ylab: y-axis label
# Output: ggplot plot object with grouped boxplot of MSEs for the three models
plot_mse_boxplot <- function(wd, analysis_dir, save_names, scenarios, 
                             scen_names, overall_name, param,
                             lower_lim = 0, upper_lim = 1, xlab, ylab) {
  if (param == "pi") {
    param_dist <- "pi_mse_all"
  } else if (param == "theta") {
    param_dist <- "theta_mse_all"
  } else if (param == "xi") {
    param_dist <- "xi_mse_all"
  } else {
    stop("Error: 'param' must be 'pi', 'theta', or 'xi' ")
  }
  
  # Combine data from various scenarios into one dataframe
  L <- 100
  num_scen <- length(scenarios)
  plot_df <- as.data.frame(matrix(NA, nrow = 3*L, ncol = num_scen))
  for (i in 1:num_scen) {
    scen_samp <- scenarios[i]
    save_name <- save_names[i]
    load(paste0(wd, analysis_dir, save_name, scen_samp, ".RData"))
    plot_df[, i] <- c(metrics_all$metrics_s[[param_dist]], 
                      metrics_all$metrics_ws[[param_dist]],
                      metrics_all$metrics_unsup[[param_dist]])
  }
  colnames(plot_df) <- scen_names
  # Add column indicating model
  plot_df$Model <- factor(c(rep("SOLCA", times=L), rep("WSOLCA", times=L), 
                            rep("WOLCA", times=L)), 
                          levels = c("SOLCA", "WSOLCA", "WOLCA"))
  
  # Create a single Sampling column by gathering sampling columns together
  # Dimensions (9*L)x3 with columns Model, overall_name, MSE
  plot_df <- plot_df %>% 
    gather({{overall_name}}, "MSE", -Model)
  plot_df[[overall_name]] <- factor(plot_df[[overall_name]], 
                                    levels = scen_names)
  
  # Create grouped boxplot of mse values for the scenarios
  p <- plot_df %>% 
    ggplot(aes_string(x = overall_name, y = "MSE", fill = "Model")) + 
    theme_bw() + 
    geom_boxplot() + 
    ylim(lower_lim, upper_lim) + 
    xlab({{xlab}}) + ylab({{ylab}})
  
  return(p)
}

# `plot_rmse_boxplot` functions the same way as `plot_mse_boxplot` but uses 
# RMSE rather than MSE
plot_rmse_boxplot <- function(wd, analysis_dir, save_names, scenarios, 
                              scen_names, overall_name, param,
                              lower_lim = 0, upper_lim = 1, xlab, ylab) {
  if (param == "pi") {
    param_dist <- "pi_mse_all"
  } else if (param == "theta") {
    param_dist <- "theta_mse_all"
  } else if (param == "xi") {
    param_dist <- "xi_mse_all"
  } else {
    stop("Error: 'param' must be 'pi', 'theta', or 'xi' ")
  }
  
  # Combine data from various scenarios into one dataframe
  L <- 100
  num_scen <- length(scenarios)
  plot_df <- as.data.frame(matrix(NA, nrow = 3*L, ncol = num_scen))
  for (i in 1:num_scen) {
    scen_samp <- scenarios[i]
    save_name <- save_names[i]
    load(paste0(wd, analysis_dir, save_name, scen_samp, ".RData"))
    plot_df[, i] <- c(sqrt(metrics_all$metrics_s[[param_dist]]), 
                      sqrt(metrics_all$metrics_unsup[[param_dist]]),
                      sqrt(metrics_all$metrics_ws[[param_dist]]))
  }
  colnames(plot_df) <- scen_names
  # Add column indicating model
  plot_df$Model <- factor(c(rep("SOLCA", times=L), rep("WOLCA", times=L), 
                            rep("WSOLCA", times=L)), 
                          levels = c("SOLCA", "WOLCA", "WSOLCA"))
  
  # Create a single Sampling column by gathering sampling columns together
  # Dimensions (9*L)x3 with columns Model, overall_name, MSE
  plot_df <- plot_df %>% 
    gather({{overall_name}}, "RMSE", -Model)
  plot_df[[overall_name]] <- factor(plot_df[[overall_name]], 
                                    levels = scen_names)
  
  # Create grouped boxplot of mse values for the scenarios
  p <- plot_df %>% 
    ggplot(aes_string(x = overall_name, y = "RMSE", fill = "Model")) + 
    theme_bw() + scale_fill_brewer(palette="Set2") + 
    geom_boxplot() + 
    ylim(lower_lim, upper_lim) + 
    xlab({{xlab}}) + ylab({{ylab}})  + 
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12))
  
  return(p)
}

#================== Plot parameters across iterations ==========================
# `theta_mode_plot` plots the modal theta patterns for each class across all 
# simulation iterations, for a given model
# Inputs: 
#   theta_plot_data: Modal patterns
#   x_label: x-axis label
# Output: Returns a ggplot object that plots the modal pattern with text 
# displaying the average mode across iterations
theta_mode_plot <- function(theta_plot_data, x_label) {
  p <- dim(theta_plot_data)[1]
  K <- dim(theta_plot_data)[2]
  Item <- factor(as.character(1:p), levels = as.character(p:1))
  theta_plot <- data.frame(theta_plot_data, Item)
  colnames(theta_plot) <- c(1:K, "Item")
  theta_plot <- theta_plot %>% gather("Class", "Level", 1:K) 
  patterns <- ggplot(theta_plot, aes(x=Class, y=Item, fill=Level)) + 
    theme_classic() +
    xlab(x_label) +
    geom_tile(color="gray") + 
    geom_text(aes(label = round(Level,2)), col="white", cex=2.5) +
    scale_fill_gradient(trans = "reverse")
  return(patterns)
}

# `plot_theta_patterns` plots the modal theta patterns for each class across
# all simulation iterations, comparing all models to the true patterns
# Inputs: 
#   wd: String specifying working directory
#   data_dir: String specifying directory where data are saved
#   analysis_dir: String specifying analysis directory where summaries are saved
#   scen_pop: Four-digit population scenario
#   scen_samp: Six-digit sample scenario
#   scen_name: String with main file name with saved metrics
# Output: Plot with the modal patterns for all models beside each other,
# with text displaying the average mode across iterations
plot_theta_patterns <- function(wd, data_dir, analysis_dir, 
                                scen_pop, iter_pop, scen_samp, scen_name) {
  # Load true population data
  load(paste0(wd, data_dir, "simdata_scen", scen_pop,"_iter", iter_pop, ".RData"))
  # Obtain true observed population parameters
  true_params <- get_true_params(sim_pop = sim_pop)  
  
  # Load simulated sample data
  load(paste0(wd, analysis_dir, scen_name, scen_samp, ".RData"))
  
  p_true <- theta_mode_plot(sim_pop$true_global_patterns, "True Classes") + 
    guides(fill = guide_legend(reverse = FALSE)) +
    labs(fill = "Modal θ Level")
  p_SOLCA <- theta_mode_plot(metrics_all$metrics_s$theta_mode, "SOLCA Classes")
  p_WOLCA <- theta_mode_plot(metrics_all$metrics_unsup$theta_mode, "WOLCA Classes")
  p_WSOLCA <- theta_mode_plot(metrics_all$metrics_ws$theta_mode, "WSOLCA Classes")
  p_comb <- ggarrange(p_true, 
                      p_SOLCA + theme(axis.title.y = element_blank()), 
                      p_WOLCA + theme(axis.title.y = element_blank()), 
                      p_WSOLCA + theme(axis.title.y = element_blank()),
                      nrow = 1, common.legend = TRUE, legend = "top")
  return(p_comb)
}

# plot pi
plot_pi_patterns <- function(wd, data_dir, analysis_dir, y_lim = c(0,1),
                             scen_pop, iter_pop, scen_samp, scen_name) {
  # Load true population data
  load(paste0(wd, data_dir, "simdata_scen", scen_pop,"_iter", iter_pop, ".RData"))
  # Obtain true observed population parameters
  true_params <- get_true_params(sim_pop = sim_pop)  
  
  # Load simulated sample data
  load(paste0(wd, analysis_dir, scen_name, scen_samp, ".RData"))
  L <- length(metrics_all$metrics_ws$K_dist)
  
  pi_plot_data <- as.data.frame(rbind(metrics_all$metrics_s$pi_all,
                                      metrics_all$metrics_unsup$pi_all,
                                      metrics_all$metrics_ws$pi_all))
  colnames(pi_plot_data) <- 1:ncol(pi_plot_data)
  pi_plot_data$Model <- c(rep("SOLCA", times=L), rep("WOLCA", times=L),
                          rep("WSOLCA", times=L))
  pi_plot_data <- pi_plot_data %>% gather("pi_component", "value", -Model)
  ggplot(pi_plot_data, aes(x=pi_component, y=value, fill=Model)) +
    theme_bw() + scale_fill_brewer(palette="RdYlBu") + 
    geom_boxplot() +
    geom_segment(mapping=aes(x=0.5, xend=1.5, y=true_params$true_pi[1],
                             yend=true_params$true_pi[1]),color="#d7191c") +
    geom_segment(mapping=aes(x=1.5, xend=2.5, y=true_params$true_pi[2],
                             yend=true_params$true_pi[2]),color="#d7191c") +
    geom_segment(mapping=aes(x=2.5, xend=3.5, y=true_params$true_pi[3],
                             yend=true_params$true_pi[3]),color="#d7191c") +
    # ggtitle(paste0("Parameter estimation for π across 100 samples")) + 
    xlab("Latent Class") + ylab("π Value") + ylim(y_lim[1],y_lim[2])
}

# plot xi
plot_xi_patterns <- function(wd, data_dir, analysis_dir, y_lim = c(0,1),
                             scen_pop, iter_pop, scen_samp, scen_name) {
  # Load true population data
  load(paste0(wd, data_dir, "simdata_scen", scen_pop,"_iter", iter_pop, ".RData"))
  # Obtain true observed population parameters
  true_params <- get_true_params(sim_pop = sim_pop)  
  
  # Load simulated sample data
  load(paste0(wd, analysis_dir, scen_name, scen_samp, ".RData"))
  L <- length(metrics_all$metrics_ws$K_dist)
  
  xi_plot_data <- as.data.frame(rbind(t(apply(metrics_all$metrics_s$xi_all, 1, c)),
                                      t(apply(metrics_all$metrics_unsup$xi_all, 1, c)),
                                      t(apply(metrics_all$metrics_ws$xi_all, 1, c))))
  colnames(xi_plot_data) <- c("S=1\nC=1", "S=1\nC=2", "S=1\nC=3",
                              "S=2\nC=1", "S=2\nC=2", "S=2\nC=3")
  xi_plot_data$Model <- c(rep("SOLCA", times=L), rep("WOLCA", times=L),
                          rep("WSOLCA", times=L))
  xi_plot_data <- xi_plot_data %>% gather("xi_component", "value", -Model)
  ggplot(xi_plot_data, aes(x=xi_component, y=value, fill=Model)) +
    theme_bw() + scale_fill_brewer(palette="RdYlBu") + 
    geom_boxplot() +
    geom_segment(mapping=aes(x=0.5, xend=1.5, y=true_params$true_xi[1,1],
                             yend=true_params$true_xi[1,1]),color="#d7191c") +
    geom_segment(mapping=aes(x=1.5, xend=2.5, y=true_params$true_xi[2,1],
                             yend=true_params$true_xi[2,1]),color="#d7191c") +
    geom_segment(mapping=aes(x=2.5, xend=3.5, y=true_params$true_xi[3,1],
                             yend=true_params$true_xi[3,1]),color="#d7191c") +
    geom_segment(mapping=aes(x=3.5, xend=4.5, y=true_params$true_xi[1,2],
                             yend=true_params$true_xi[1,2]),color="#d7191c") +
    geom_segment(mapping=aes(x=4.5, xend=5.5, y=true_params$true_xi[2,2],
                             yend=true_params$true_xi[2,2]),color="#d7191c") +
    geom_segment(mapping=aes(x=5.5, xend=6.5, y=true_params$true_xi[3,2],
                             yend=true_params$true_xi[3,2]),color="#d7191c") +
    # ggtitle(paste0("Parameter estimation for ξ across 100 samples")) + 
    xlab("Covariate Levels") + ylab("ξ Value") + ylim(y_lim[1],y_lim[2])
}

# Plot Phi
plot_Phi_patterns <- function(wd, data_dir, analysis_dir, y_lim = c(0,1),
                             scen_pop, iter_pop, scen_samp, scen_name) {
  # Load true population data
  load(paste0(wd, data_dir, "simdata_scen", scen_pop,"_iter", iter_pop, ".RData"))
  # Obtain true observed population parameters
  true_params <- get_true_params(sim_pop = sim_pop)  
  
  # Load simulated sample data
  load(paste0(wd, analysis_dir, scen_name, scen_samp, ".RData"))
  L <- length(metrics_all$metrics_ws$K_dist)
  
  Phi_plot_data <- as.data.frame(rbind(pnorm(t(apply(metrics_all$metrics_s$xi_all, 1, c))),
                                      pnorm(t(apply(metrics_all$metrics_unsup$xi_all, 1, c))),
                                      pnorm(t(apply(metrics_all$metrics_ws$xi_all, 1, c)))))
  colnames(Phi_plot_data) <- c("S=1\nC=1", "S=1\nC=2", "S=1\nC=3",
                               "S=2\nC=1", "S=2\nC=2", "S=2\nC=3")
  Phi_plot_data$Model <- c(rep("SOLCA", times=L), rep("WOLCA", times=L),
                          rep("WSOLCA", times=L))
  Phi_plot_data <- Phi_plot_data %>% gather("Phi_component", "value", -Model)
  true_Phi <- true_params$true_Phi_mat
  ggplot(Phi_plot_data, aes(x=Phi_component, y=value, fill=Model)) +
    theme_bw() + scale_fill_brewer(palette="RdYlBu") + 
    geom_boxplot() +
    ylab("P(Y=1|S,C)") +
    geom_segment(mapping=aes(x=0.5, xend=1.5, y=true_Phi[1,1],
                             yend=true_Phi[1,1]),color="#d7191c") +
    geom_segment(mapping=aes(x=1.5, xend=2.5, y=true_Phi[2,1],
                             yend=true_Phi[2,1]),color="#d7191c") +
    geom_segment(mapping=aes(x=2.5, xend=3.5, y=true_Phi[3,1],
                             yend=true_Phi[3,1]),color="#d7191c") +
    geom_segment(mapping=aes(x=3.5, xend=4.5, y=true_Phi[1,2],
                             yend=true_Phi[1,2]),color="#d7191c") +
    geom_segment(mapping=aes(x=4.5, xend=5.5, y=true_Phi[2,2],
                             yend=true_Phi[2,2]),color="#d7191c") +
    geom_segment(mapping=aes(x=5.5, xend=6.5, y=true_Phi[3,2],
                             yend=true_Phi[3,2]),color="#d7191c") +
    # ggtitle(paste0("Parameter estimation for outcome probabilities across 100 samples")) + 
      xlab("Covariate Levels") + ylab("Outcome Probability") + ylim(y_lim[1],y_lim[2])
}

#============== Plot parameters, by latent class, across sampling designs ======
# `plot_pi_samp` plots pi values across iterations for the three latent classes
# as a grouped boxplot for the three models, and facetted by sampling design
# Inputs:
#   wd: String specifying working directory
#   analysis_dir: String specifying analysis directory where summaries are saved
#   scen_pop: 4-digit number specifying population scenario
#   iter_pop: Number specifying population iteration. Default is 1
# Output: ggplot plot object with grouped boxplot of pis for the three models
plot_pi_samp <- function(wd, analysis_dir, scen_pop, iter_pop = 1) {
  
  # Load true population data
  load(paste0(wd, data_dir, "simdata_scen", scen_pop,"_iter", iter_pop, ".RData"))
  # Obtain true observed population parameters
  true_params <- get_true_params(sim_pop = sim_pop)
  
  # Load saved metrics
  load(paste0(wd, analysis_dir, "metrics_scen", 111213, ".RData"))
  metrics_SRS_s <- metrics_all$metrics_s
  metrics_SRS_ws <- metrics_all$metrics_ws
  metrics_SRS_unsup <- metrics_all$metrics_unsup
  load(paste0(wd, analysis_dir, "metrics_scen", 111211, ".RData"))
  metrics_Strat_s <- metrics_all$metrics_s
  metrics_Strat_ws <- metrics_all$metrics_ws
  metrics_Strat_unsup <- metrics_all$metrics_unsup
  load(paste0(wd, analysis_dir, "metrics_scen", 111212, ".RData"))
  metrics_Clus_s <- metrics_all$metrics_s
  metrics_Clus_ws <- metrics_all$metrics_ws
  metrics_Clus_unsup <- metrics_all$metrics_unsup
  
  # Create plot data
  L <- length(metrics_SRS_ws$K_dist)
  pi_plot_data <- as.data.frame(rbind(metrics_SRS_s$pi_all,
                                      metrics_Strat_s$pi_all,
                                      metrics_Clus_s$pi_all,
                                      metrics_SRS_unsup$pi_all,
                                      metrics_Strat_unsup$pi_all,
                                      metrics_Clus_unsup$pi_all,
                                      metrics_SRS_ws$pi_all,
                                      metrics_Strat_ws$pi_all,
                                      metrics_Clus_ws$pi_all))
  colnames(pi_plot_data) <- 1:ncol(pi_plot_data)
  pi_plot_data$Model <- c(rep("SOLCA", times=3*L), 
                          rep("WOLCA", times=3*L),
                          rep("WSOLCA", times=3*L))
  pi_plot_data$Sampling <- factor(rep(c("SRS", "Stratified", "Stratified Cluster"), 
                                      each = L, times = 3),
                                  levels = c("SRS", "Stratified", 
                                             "Stratified Cluster"))
  pi_plot_data <- pi_plot_data %>% 
    gather("pi_component", "value", -c(Model, Sampling))
  
  # Create and return plot
  gg <- ggplot(pi_plot_data, aes(x=pi_component, y=value, fill=Model)) +
    theme_bw() + scale_fill_brewer(palette="Set2") + 
    geom_boxplot() +
    facet_grid(~Sampling) + 
    geom_segment(mapping=aes(x=0.5, xend=1.5, y=true_params$true_pi[1],
                             yend=true_params$true_pi[1]),color="#d7191c") +
    geom_segment(mapping=aes(x=1.5, xend=2.5, y=true_params$true_pi[2],
                             yend=true_params$true_pi[2]),color="#d7191c") +
    geom_segment(mapping=aes(x=2.5, xend=3.5, y=true_params$true_pi[3],
                             yend=true_params$true_pi[3]),color="#d7191c") +
    # ggtitle(paste0("Parameter estimation for π across 100 samples")) + 
    xlab("Latent Class") + ylab("π Value") + 
    theme(legend.position="top")  + 
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12))
  return(gg)
}

plot_xi_samp <- function(wd, analysis_dir, scen_pop, iter_pop = 1) {
  # Load true population data
  load(paste0(wd, data_dir, "simdata_scen", scen_pop,"_iter", iter_pop, ".RData"))
  # Obtain true observed population parameters
  true_params <- get_true_params(sim_pop = sim_pop)
  
  # Load saved metrics
  load(paste0(wd, analysis_dir, "metrics_scen", 111213, ".RData"))
  metrics_SRS_s <- metrics_all$metrics_s
  metrics_SRS_ws <- metrics_all$metrics_ws
  metrics_SRS_unsup <- metrics_all$metrics_unsup
  load(paste0(wd, analysis_dir, "metrics_scen", 111211, ".RData"))
  metrics_Strat_s <- metrics_all$metrics_s
  metrics_Strat_ws <- metrics_all$metrics_ws
  metrics_Strat_unsup <- metrics_all$metrics_unsup
  load(paste0(wd, analysis_dir, "metrics_scen", 111212, ".RData"))
  metrics_Clus_s <- metrics_all$metrics_s
  metrics_Clus_ws <- metrics_all$metrics_ws
  metrics_Clus_unsup <- metrics_all$metrics_unsup
  
  # Create plot data
  L <- length(metrics_SRS_ws$K_dist)
  xi_plot_data <- as.data.frame(rbind(t(apply(metrics_SRS_s$xi_all, 1, c)),
                                      t(apply(metrics_Strat_s$xi_all, 1, c)),
                                      t(apply(metrics_Clus_s$xi_all, 1, c)),
                                      t(apply(metrics_SRS_unsup$xi_all, 1, c)),
                                      t(apply(metrics_Strat_unsup$xi_all, 1, c)),
                                      t(apply(metrics_Clus_unsup$xi_all, 1, c)),
                                      t(apply(metrics_SRS_ws$xi_all, 1, c)),
                                      t(apply(metrics_Strat_ws$xi_all, 1, c)),
                                      t(apply(metrics_Clus_ws$xi_all, 1, c))))
  colnames(xi_plot_data) <- c("S=1\nC=1", "S=1\nC=2", "S=1\nC=3",
                              "S=2\nC=1", "S=2\nC=2", "S=2\nC=3")
  xi_plot_data$Model <- c(rep("SOLCA (Unweighted)", times=3*L), 
                          rep("WOLCA (Two-Step)", times=3*L),
                          rep("WSOLCA", times=3*L))
  xi_plot_data$Sampling <- factor(rep(c("SRS", "Stratified", "Stratified Cluster"), 
                                      each = L, times = 3),
                                  levels = c("SRS", "Stratified", 
                                             "Stratified Cluster"))
  xi_plot_data <- xi_plot_data %>% 
    gather("xi_component", "value", -c(Model, Sampling))
  
  # Create and return plot
  gg <- ggplot(xi_plot_data, aes(x=xi_component, y=value, fill=Model)) +
    theme_bw() + scale_fill_brewer(palette="Set2") + 
    geom_boxplot() +
    facet_grid(~Sampling) + 
    geom_segment(mapping=aes(x=0.5, xend=1.5, y=true_params$true_xi[1,1],
                             yend=true_params$true_xi[1,1]),color="#d7191c") +
    geom_segment(mapping=aes(x=1.5, xend=2.5, y=true_params$true_xi[2,1],
                             yend=true_params$true_xi[2,1]),color="#d7191c") +
    geom_segment(mapping=aes(x=2.5, xend=3.5, y=true_params$true_xi[3,1],
                             yend=true_params$true_xi[3,1]),color="#d7191c") +
    geom_segment(mapping=aes(x=3.5, xend=4.5, y=true_params$true_xi[1,2],
                             yend=true_params$true_xi[1,2]),color="#d7191c") +
    geom_segment(mapping=aes(x=4.5, xend=5.5, y=true_params$true_xi[2,2],
                             yend=true_params$true_xi[2,2]),color="#d7191c") +
    geom_segment(mapping=aes(x=5.5, xend=6.5, y=true_params$true_xi[3,2],
                             yend=true_params$true_xi[3,2]),color="#d7191c") +
    # ggtitle(paste0("Parameter estimation for ξ across 100 samples")) + 
    xlab("Covariate Levels") + ylab("ξ Value") + 
    theme(legend.position="top")  + 
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12))
  return(gg)
}

#============== Create Appendix Tables =========================================

# `create_app_tables` uses 'gt()' to create a table of summary metrics comparing
# all groups of simulation scenarios
# Inputs:
#   wd: String specifying working directory
#   analysis_dir: String specifying analysis directory where summaries are saved
#   save_names: String vector with main file names for saving all metrics
#   scenarios: Numeric vector specifying sampling scenarios to include
#   scen_names: String vector specifying labels for sampling scenarios
#   overall_name: String specifying overall grouping name for sampling scenarios
#   format: String specifying kable table format. Must be "latex", "html", 
# "pipe", "simple", or "rst". Default is "latex"
# Output: Formatted table with absolute bias, CI width, and coverage
create_app_tables <- function(wd, analysis_dir, save_names, scenarios, 
                              scen_names, overall_name, format = "latex") {
  num_scen <- length(scenarios)
  
  metrics_summ <- as.data.frame(matrix(NA, nrow = 3*length(scenarios), 
                                       ncol = 12))
  colnames(metrics_summ) <- c(overall_name, "Model", 
                              "$K$ Abs Bias", "$\\pi$ Abs Bias",  
                              "$\\theta$ Abs Bias", "$\\xi$ Abs Bias", 
                              "$\\pi$ CI Width", "$\\theta$ CI Width", "$\\xi$ CI Width", 
                              "$\\pi$ Coverage","$\\theta$ Coverage", "$\\xi$ Coverage")
  metrics_summ[, 1] <- rep(scen_names, each = 3)
  metrics_summ[, 2] <- rep(c("SOLCA", "WOLCA", "WSOLCA"), num_scen)  
  # output_inds <- 1:7
  output_inds <- c(1, 2, 4, 6, 3, 5, 7)
  for (i in 1:num_scen) {
    scen_samp <- scenarios[i]
    save_name <- save_names[i]
    load(paste0(wd, analysis_dir, save_name, scen_samp, ".RData"))
    row_ind <- 3 * (i-1)
    metrics_summ[row_ind + 1, -c(1,2)] <- c(metrics_all$metrics_s[output_inds], 
                                  mean(metrics_all$metrics_s$pi_cover_avg), 
                                  mean(metrics_all$metrics_s$theta_cover_avg),
                                  mean(metrics_all$metrics_s$xi_cover_avg))
    metrics_summ[row_ind + 2, -c(1,2)] <- c(metrics_all$metrics_unsup[output_inds], 
                                            mean(metrics_all$metrics_unsup$pi_cover_avg), 
                                            mean(metrics_all$metrics_unsup$theta_cover_avg),
                                            mean(metrics_all$metrics_unsup$xi_cover_avg))
    metrics_summ[row_ind + 3, -c(1,2)] <- c(metrics_all$metrics_ws[output_inds], 
                                  mean(metrics_all$metrics_ws$pi_cover_avg), 
                                  mean(metrics_all$metrics_ws$theta_cover_avg),
                                  mean(metrics_all$metrics_ws$xi_cover_avg))
  }
  
  metrics_summ %>% 
    kbl(digits = 4, align = "rrrrrrrrrrrr", booktabs = TRUE, format = format,
        caption = "Summary of mean absolute bias, 95% credible interval width, and coverage for simulations based on posterior samples.") %>%
    kable_classic() %>%
    kable_styling(full_width = FALSE)
}

