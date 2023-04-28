#===================================================================
# Simulate data with stratum-specific class-membership probabilities
# Informative sampling: strata variable influences class and outcome
# Programmer: SM Wu
# Date updated: 2023/04/10
#===================================================================

library(LaplacesDemon)
library(SimCorMultRes)
library(dplyr)

# `convert_to_reference` converts a probit coefficient matrix to a vector of 
# coefficients using reference cell coding
# Inputs:
#   xi: Matrix of probit coefficients in factor variable coding; Kxq
# Outputs:
#   beta: Vector of probit coefficients in reference cell coding with interactions; (K*q)x1
convert_to_reference <- function(xi) {
  # Convert to b1 + b2*I(C=2) + b3*I(C=3) + b4*I(S=2) + b5*I(C=2,S=2) + b6*I(C=3,S=2)
  xi_vec <- c(xi)
  beta <- numeric(length(xi_vec))
  beta[1] <- xi_vec[1]
  beta[2] <- xi_vec[2] - xi_vec[1]
  beta[3] <- xi_vec[3] - xi_vec[1]
  beta[4] <- xi_vec[4] - xi_vec[1]
  beta[5] <- xi_vec[5] - xi_vec[4] - beta[2]
  beta[6] <- xi_vec[6] - xi_vec[4] - beta[3]
  return(beta)
}

# `create_pop` creates and saves simulated population data
# Inputs:
#   scenario: Four-digit population scenario. First digit specifies strength of 
# pattern, second digit specifies diet-outcome relationship, third digit 
# specifies selection (stratum) variable type, fourth digit specifies clustering
#   iter_pop: Population iteration. Default is 1.
#   pop_data_path: Path for saving simulated population data
# Outputs: 'sim_pop' list with the following elements:
#   N: Population size
#   N_s: Vector of stratum population sizes. Sx1
#   p: Number of exposure items
#   d: Number of exposure categories
#   S: Number of strata
#   true_K: True number of latent classes, K
#   true_pi_s: Matrix of true pis within each stratum. SxK  
#   true_pi: Vector of true pis overall. Kx1
#   true_global_patterns: Matrix of true global exposure patterns defined by
# modal category. pxK
#   true_global_thetas: Array of true thetas. pxKxd 
#   true_Si = true_Si: Vector of true stratum indicators. Nx1 
#   true_Ci = true_Ci: Vector of true latent class indicators. Nx1 
#   cluster_id: Vector of true cluster indicators. Nx1 
#   cluster_size: Cluster size, same for all clusters
#   X_data = X_data: Matrix of categorical exposure for all individuals. Nxp
#   Y_data = Y_data: Vector of binary outcome for all individuals. Nx1
#   true_xi: Matrix of true xis. Kxq
#   true_Phi_mat: Matrix of true outcome probabilities, pnorm(true_xi). Kxq
#   true_Phi: Vector of true outcome probabilities for all individuals. Nx1
create_pop <- function(scenario, iter_pop = 1, pop_data_path) {
  # Set seed
  set.seed(iter_pop)
  # Extract scenario indicators
  scenario_vec <- as.numeric(unlist(strsplit(as.character(scenario), "")))
  
  # Dimensions
  p <- 30                  # Number of food items
  d <- 4                   # Number of response levels (assumed constant across items)
  S <- 2                   # Number of subpops
  K <- 3                   # Number of global classes
  N_s <- c(20000, 60000)   # Stratum sizes
  N <- sum(N_s)            # Population size    
  
  #================ Specify strength of dietary patterns =======================
  if (scenario_vec[1] == 1) {
    # Strength of dietary pattern is 85:5:5:5
    clust_mode <- 0.85       # Probability of true consumption level occurring
    non_mode <- 0.05         # Probability of other consumption levels occurring
  } else if (scenario_vec[1] == 2) {
    # Strength of dietary pattern is 55:15:15:15
    clust_mode <- 0.55
    non_mode <- 0.15
  } else {
    stop("Error: Not a valid scenario.")
  }
  
  #================ Specify diet-outcome relationship ==========================
  global_1 <- c(rep(1, times = 0.5 * p), rep(3, times = 0.5 * p))
  global_2 <- c(rep(4, times = 0.2 * p), rep(2, times = 0.8 * p))
  if (scenario_vec[2] == 1) {
    # Diet-outcome relationship set to Baseline
    global_3 <- c(rep(3, times = 0.3 * p), rep(4, times = 0.4 * p), 
                  rep(1, times = 0.3 * p))
  } else if (scenario_vec[2] == 2) {
    # Diet-outcome relationship set to Supervised
    global_3 <- c(rep(4, times = 0.1 * p), rep(1, times = 0.4 * p), 
                  rep(3, times = 0.5 * p))
  } else {
    stop("Error: Not a valid scenario.")
  } 
  #================ True modal patterns for each class =========================
  true_global_patterns <- cbind(global_1, global_2, global_3)
  
  #================ Specify selection (stratum) variable type ==================
  if (scenario_vec[3] == 1) {
    # Selection variable S is a confounder associated with X and Y
    true_pi_s <- matrix(c(0.2, 0.4, 0.4,
                          0.7, 0.2, 0.1), nrow = 2, byrow = TRUE)
  } else if (scenario_vec[3] == 2) {
    # Selection variable S is a precision variable associated with Y
    true_pi_s <- matrix(c(0.2, 0.4, 0.4,
                          0.2, 0.4, 0.4), nrow = 2, byrow = TRUE)
  } else {
    stop("Error: Not a valid scenario.")
  }  
  #================ True parameter values for class membership probs, pi =======
  true_pi <- c(t(true_pi_s) %*% as.matrix(N_s/N))
  
  #================ True class and stratum variables for all individuals =======
  pop_inds <- 1:N
  true_Si <- c(rep(1, times = N_s[1]), rep(2, times = N_s[2]))
  true_Ci <- rep(1, N)
  for (s in 1:S) {
    pop_s <- pop_inds[true_Si == s]
    for (k in 2:K) {
      inds_k <- sample(pop_s, round(N_s[s] * true_pi_s[s, k]))
      true_Ci[inds_k] <- k
      pop_s <- setdiff(pop_s, inds_k)
    }
  }
  
  #================ True parameter values for item consumption probs, theta ====
  true_global_thetas <- array(non_mode, dim = c(p, K, d))
  for (j in 1:p) {
    for (k in 1:K) {
      true_global_thetas[j, k, true_global_patterns[j, k]] <- clust_mode
    }
  }
  
  #================ Categorical exposure =======================================
  X_data <- matrix(0, nrow = N, ncol = p)
  for (i in 1:N) {
    for (j in 1:p) {
      X_data[i, j] <- rcat(n = 1, p = true_global_thetas[j, true_Ci[i], ]) 
    }
  }
  
  #================ Binary outcome and true outcome probabilities ==============
  if (scenario_vec[4] == 1) {
    # No clustering in the data
    cluster_id <- pop_inds
    cluster_size <- 1
    xi_under <- matrix(c(1, 0.3, -0.5, 0.5, -0.7, -1.3), nrow = 3, byrow = FALSE)
    # Matrix of true outcome probabilities
    true_Phi_mat <- pnorm(xi_under) 
    # True outcome probabilities for each individual
    true_Phi <- numeric(N)   
    # Outcome data
    Y_data <- numeric(N) 
    for (s in 1:S) {
      for (k in 1:K) {
        # Individuals in stratum s class k
        N_s_k <- pop_inds[true_Si == s & true_Ci == k]
        with_outcome <- sample(N_s_k, round(length(N_s_k) * true_Phi_mat[k, s]))
        Y_data[with_outcome] <- 1
        true_Phi[N_s_k] <- true_Phi_mat[k, s]
      }
    }
  } else if (scenario_vec[4] == 2) {
    # Clustered data with 1000 clusters of size 80
    # Simulate correlated binary outcomes 'SimCorMultRes' package
    cluster_size <- 80  # Cluster size
    # Assume exchangeable correlation matrix
    latent_correlation_matrix <- toeplitz(c(1, rep(0.5, cluster_size - 1)))
    # Underlying true parameter values for probit model coefficients, xi
    xi_under <- matrix(c(1, 0.3, -0.5, 0.5, -0.7, -1.3), nrow = 3, byrow = FALSE)
    # Convert xi to cell reference coding beta coefficients
    beta_coefs <- convert_to_reference(xi_under)
    intercepts <- beta_coefs[1]
    betas <- beta_coefs[-1]
    # Create probit regression design matrix with interactions
    C2 <- ifelse(true_Ci == 2, 1, 0)
    C3 <- ifelse(true_Ci == 3, 1, 0)
    S2 <- ifelse(true_Si == 2, 1, 0)
    V_design <- data.frame(C2 = C2, C3 = C3, S2 = S2, C2S2 = C2*S2, C3S2 = C3*S2)
    # Simulate correlated binary outcomes
    sim_binary <- rbin(clsize = cluster_size, intercepts = intercepts,
                       betas = betas, xformula = ~C2+C3+S2+C2S2+C3S2, 
                       xdata = V_design,
                       cor.matrix = latent_correlation_matrix, link = "probit")
    # Cluster indicator for all individuals. Stratum 1 includes clusters 1-250.
    # Stratum 2 includes clusters 251-1000
    cluster_id <- sim_binary$simdata$id  
    # Outcome data
    Y_data <- sim_binary$simdata$y       
    # Matrix of true outcome probabilities
    true_Phi_mat <- matrix(NA, nrow=K, ncol=S)  
    # True outcome probabilities for each individual
    true_Phi <- numeric(N)               
    for (s in 1:S) {
      for (k in 1:K) {
        # Individuals in stratum s class k
        N_s_k <- pop_inds[true_Si == s & true_Ci == k]  
        true_Phi_mat[k, s] <- sum(Y_data[N_s_k]) / length(N_s_k)
        true_Phi[N_s_k] <- true_Phi_mat[k, s]
      }
    }
  } else {
    stop("Error: Not a valid scenario.")
  }
  
  #================ True parameter values for probit model coefficients, xi ====
  # True xi in the finite population
  true_xi <- qnorm(true_Phi_mat)   
  
  #================ Save and return data =======================================
  sim_pop <- list(N = N, N_s = N_s, p = p, d = d, S = S, true_K = K, 
                  true_pi_s = true_pi_s, true_pi = true_pi, 
                  true_global_patterns = true_global_patterns,
                  true_global_thetas = true_global_thetas, true_Si = true_Si, 
                  true_Ci = true_Ci, cluster_id = cluster_id, 
                  cluster_size = cluster_size, X_data = X_data, Y_data = Y_data, 
                  true_xi = true_xi, true_Phi_mat = true_Phi_mat, true_Phi = true_Phi)
  save(sim_pop, file = pop_data_path)
  return(sim_pop)
}

# `create_samp` creates and saves simulated sample data
# Inputs:
#   sim_pop: Simulated population output from `create_pop` function
#   scenario: Six-digit sample scenario. First digit specifies strength of 
# pattern, second digit specifies diet-outcome relationship, third digit 
# specifies selection (stratum) variable type, fourth digit specifies clustering,
# fifth digit specifies sample size, sixth digit specifies sampling design
#   samp_n: Sample iteration
#   samp_data_path: Path for saving simulated sample data
# Outputs: 'sim_data' list with the following elements:
#   samp_ind: Vector of population indices for sampled individuals. nx1
#   sample_wt: Vector of sampling weights for sampled individuals. nx1
#   true_Si = true_Si: Vector of stratum indicators. Nx1 
#   cluster_id: Vector of true cluster indicators. Nx1 
#   X_data = X_data: Matrix of categorical exposure for all individuals. Nxp
#   Y_data = Y_data: Vector of binary outcome for all individuals. Nx1
#   N: Population size
#   N_s: Vector of stratum population sizes. Sx1
#   p: Number of exposure items
#   d: Number of exposure categories
#   S: Number of strata
#   true_K: True number of latent classes, K
#   true_pi_s: Matrix of true pis within each stratum. SxK  
#   true_pi: Vector of true pis overall. Kx1
#   true_global_patterns: Matrix of true global exposure patterns defined by
# modal category. pxK
#   true_global_thetas: Array of true thetas. pxKxd 
#   true_xi: Matrix of true xis. Kxq
#   true_Ci = true_Ci: Vector of true latent class indicators. Nx1 
#   true_Phi_mat: Matrix of true outcome probabilities, pnorm(true_xi). Kxq
#   true_Phi: Vector of true outcome probabilities for all individuals. Nx1
create_samp <- function(sim_pop, scenario, samp_n, samp_data_path) {
  set.seed(100 + samp_n)
  # Extract scenario indicators
  scenario_vec <- as.numeric(unlist(strsplit(as.character(scenario), "")))
  
  #================ Specify sample size ========================================
  
  if (scenario_vec[5] == 1) {
    # Sample size is 5%
    n <- ceiling(0.05 * sim_pop$N)
  } else if (scenario_vec[5] == 2) {
    # Sample size is 10%
    n <- ceiling(0.1 * sim_pop$N)
  } else if (scenario_vec[5] == 3) {
    # Sample size is 1%
    n <- ceiling(0.01 * sim_pop$N)
  } else {
    stop("Error: Not a valid scenario.")
  } 
  
  #================ Specify sampling design ====================================
  sample_wt_temp <- numeric(sim_pop$N)  # Temp weights for population
  samp_ind_temp <- numeric(sim_pop$N)   # Temp sample indicator for population
  pop_inds <- 1:sim_pop$N
  clus_inds <- 1:max(sim_pop$cluster_id)
  if (scenario_vec[6] == 1) {
    # Survey design is stratified sampling with unequal probabilities
    n_s <- c(n, n) / sim_pop$S  # Sample strata are equal sizes by design
    for (s in 1:sim_pop$S) {
      # Individuals in stratum s
      pop_s <- pop_inds[sim_pop$true_Si == s]  
      # Sample from stratum s
      samp_ind_s <- sample(pop_s, n_s[s])   
      samp_ind_temp[samp_ind_s] <- 1          
      sample_wt_temp[samp_ind_s] <- sim_pop$N_s[s] / n_s[s]
    }
  } else if (scenario_vec[6] == 2) {
    # Survey design is SRS
    # Sampled individuals
    samp_ind_SRS <- sample(pop_inds, n)  
    samp_ind_temp[samp_ind_SRS] <- 1
    sample_wt_temp[samp_ind_SRS] <- sim_pop$N / n
  } else if (scenario_vec[6] == 3) {
    # Survey design is stratified cluster sampling
    n_s <- c(n/2, n/2)  # Sample strata are equal sizes by design
    n_clus_s <- n_s / sim_pop$cluster_size  # Number of clusters to sample per stratum
    for (s in 1:sim_pop$S) {
      # Clusters in stratum s
      clus_s <- unique(sim_pop$cluster_id[sim_pop$true_Si == s])  
      # Sample clusters from stratum s
      samp_clus_s <- sample(clus_s, n_clus_s[s])    
      # Sampled individuals
      samp_ind_s <- pop_inds[sim_pop$cluster_id %in% samp_clus_s]
      samp_ind_temp[samp_ind_s] <- 1          
      sample_wt_temp[samp_ind_s] <- sim_pop$N_s[s] / n_s[s]
    }
  } else {
    stop("Error: Not a valid scenario.")
  }
  
  #================ Subset to sampled individuals ==============================
  samp_ind <- pop_inds[samp_ind_temp > 0]
  sample_wt <- sample_wt_temp[samp_ind]
  X_data <- sim_pop$X_data[samp_ind, ]
  Y_data <- sim_pop$Y_data[samp_ind]
  cluster_id <- sim_pop$cluster_id[samp_ind]
  true_Si <- sim_pop$true_Si[samp_ind]
  true_Ci <- sim_pop$true_Ci[samp_ind]
  
  #================ Save and return output =====================================
  sim_data <- list(samp_ind = samp_ind, sample_wt = sample_wt, true_Si = true_Si, 
                   cluster_id = cluster_id, X_data = X_data, Y_data = Y_data, 
                   N = sim_pop$N, N_s = sim_pop$N_s, p = sim_pop$p, 
                   d = sim_pop$d, S = sim_pop$S, true_K = sim_pop$true_K, 
                   true_pi_s = sim_pop$true_pi_s, true_pi = sim_pop$true_pi, 
                   true_global_patterns = sim_pop$true_global_patterns,
                   true_global_thetas = sim_pop$true_global_thetas, 
                   true_xi = sim_pop$true_xi, true_Ci = true_Ci, 
                   true_Phi_mat = sim_pop$true_Phi_mat, 
                   true_Phi = sim_pop$true_Phi)
  save(sim_data, file = samp_data_path)
  return(sim_data)
}


# Define directories
wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/"
# wd <- "~/Documents/Harvard/Research/Briana/supRPC/wsOFMM/"
data_dir <- "Data/"

#==================== Create population scenarios ==============================
scenarios <- 1211
scenarios <- c(1111, 2111, 1121)
iter_pop <- 1
scenarios <- c(1111, 2111, 1211, 1121, 1112)
for (scenario in scenarios) {
  pop_data_path <- paste0(wd, data_dir, "simdata_scen", scenario, "_iter", 
                          iter_pop, ".RData") 
  sim_pop <- create_pop(scenario = scenario, iter_pop = iter_pop, 
                        pop_data_path = pop_data_path)
}

# Sanity checks
prop.table(table(sim_pop$true_Ci[sim_pop$true_Si == 1]))
prop.table(table(sim_pop$true_Ci[sim_pop$true_Si == 2]))
sim_pop$true_Phi_mat

#==================== Create sampling scenarios ================================
scenarios <- 121111
scenarios <- c(211111, 112111, 111121, 111131)
scenarios <- 112111
samp_n_seq <- 1
scenarios <- c(111111, 211111, 121111, 112111, 111211, 111121, 111131, 111112, 
               111113)
samp_n_seq <- 1:100
for (scenario in scenarios) {
  # Get population scenario
  scen_pop <- substring(scenario, 1, 4)
  sim_pop_path <- paste0(wd, data_dir, "simdata_scen", scen_pop, "_iter", 
                         iter_pop, ".RData")
  load(sim_pop_path)
  for (samp_n in samp_n_seq) {
    samp_data_path <- paste0(wd, data_dir, "simdata_scen", scenario, "_iter", 
                             iter_pop, "_samp", samp_n, ".RData") 
    sim_data <- create_samp(sim_pop = sim_pop, scenario = scenario, 
                            samp_n = samp_n, samp_data_path = samp_data_path)
  }
}

#==================== Sanity checks ============================================
# table(true_Si)
# prop.table(table(true_Ci))
# all(apply(true_global_thetas, c(1, 2), sum) == 1)
# apply(X_data[true_Ci == 1, ], 2, median)
# apply(X_data[true_Ci == 1, ], 2, mean)
# for (k in 1:K) {print(mean(Y_data[true_Ci == k]))}
# for (s in 1:S) {print(mean(Y_data[true_Si == s]))}
# sum(beta[c(1,2)]) == xi_vec[2]
# sum(beta[c(1,2,4,5)]) == xi_vec[5]
# sum(beta[c(1,3,4,6)]) == xi_vec[6]
# mean(sim_binary$simdata[sim_binary$simdata$C2S2 == 1, "y"])
# mean(sim_binary$simdata[sim_binary$simdata$C3S2 == 1, "y"])
# mean(sim_binary$simdata[sim_binary$simdata$S2 == 1, "y"])
# table(sim_binary$simdata$time)
# # Check that no cluster is in both strata
# intersect(unique(cluster_id[true_Si == 1]), unique(cluster_id[true_Si == 2]))
# table(sim_data$sample_wt)
# prop.table(table(sim_data$true_Si))
# prop.table(table(sim_data$true_Ci))
