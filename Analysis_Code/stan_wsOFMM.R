# Running WSOLCA using Stan
# Author: Stephanie Wu
# Date updated: 2023/03/08

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

# Create Stan model
mod_stan <- stan_model(paste0(analysis_dir, "mixture_model.stan"))

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
# xi_mat <- matrix(analysis$xi_med, byrow = FALSE, nrow = K)  # KxS matrix
xi_mat <- matrix(colMeans(analysis$xi_red), byrow = FALSE, nrow = K)  # KxS matrix

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
# xi_med_adj <- apply(xi_red_adj, c(2,3), median)
xi_med_adj <- apply(xi_red_adj, c(2,3), mean)