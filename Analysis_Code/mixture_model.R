library(csSampling)
library(rstan)
library(survey)
library(tidyverse)
library(R.matlab)
set.seed(11152022)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

######### Helper functions ###################
## 'grad_par()' helper function nested withi 'withReplicates()' to obtain gradient.
## Stan will pass warnings from calling 0 chains, but will still create an 
## out_stan object for the 'grad_log_prob()' method
grad_par <- function(pwts, svydata, stanmod, standata, par_stan, par_hat) {
  standata$weights <- pwts
  out_stan <- sampling(object = stanmod, data = standata, pars = par_stan,
                       chains = 0, warmup = 0)
  gradpar <- grad_log_prob(out_stan, par_hat)
  return(gradpar)
}

## 'DEadj()' helper function to apply matrix rotation
DEadj <- function(par, par_hat, R2R1) {
  par_adj <- (par - par_hat) %*% R2R1 + par_hat
  return(par_adj)
}

setwd("~/Documents/Harvard/Research/Briana/supRPC/wsOFMM")
data_dir <- "Data/"
res_dir <- "Results/"
analysis_dir <- "Analysis_Code/"
model <- "wsOFMM"
scen_samp <- 101
iter_pop <- 1
samp_n <- 1

### ADD IN FOR LOOP ACROSS ITERATIONS
### NEED SAMPLED DATA AS AN OUTPUT

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
} else {
  # Load simulated sample data for the iteration
  sim_samp <- readMat(sim_samp_path)$sim.data
  names(sim_samp) <- str_replace_all(dimnames(sim_samp)[[1]], "[.]", "_")
  
  # Load model output and extract analysis portion
  output <- readMat(sim_res_path)
  analysis <- output$analysis
  names(analysis) <- str_replace_all(dimnames(analysis)[[1]], "[.]", "_")
} 

# Posterior estimates from MCMC sampler
K <- c(analysis$k_red)
p <- dim(analysis$theta_med)[1]
d <- dim(analysis$theta_med)[3]
n <- length(analysis$c_i)
q <- length(analysis$xi_med)
X <- sim_samp$X_data
y <- c(sim_samp$Y_data)
S <- 2

# Change xi to reference cell coding 
# No need to adjust for label switching
xi_ref <- analysis$xi_med[1]
for (i in 2:q) {
  if (i < (K+S)) {
    xi_ref[i] <- analysis$xi_med[i] - analysis$xi_med[1]
  } else {
    xi_ref[i] <- analysis$xi_med[i] - analysis$xi_med[K+S-1] - xi_ref[i - K]
  }
}

# Converting theta array to vector order: dim j -> dim k -> dim r
# First 1:p for k=1 and d=1, then 1:p for k=2 and d=1, then 1:p for k=3 and d=1, 
# then 1:p for k=1 and d=2, then 1:p for k=2 and d=2, then 1:p for k=3 and d=2,... 
names_array <- array(1:(p*K*d), dim =c(p, K, d))
par_hat <- c(analysis$pi_med, c(analysis$theta_med), xi_ref)
names(par_hat) <- c(paste0("pi", 1:K), 
                    paste0("theta", 1:p, "_", rep(1:K, each=p), "_", rep(1:d, each=p*K)),
                    paste0("xi", 1:q))
## Can I run grad_log_prob on a stan model with 0 chains? Yes, see grad_par function
## How do I plug in the posterior estimates from the MCMC sampler? Use par_hat


# Create Stan model
mod_stan <- stan_model(paste0(analysis_dir, "mixture_model.stan"))

# Stan parameters of interest
par_stan <- c('pi', 'theta', 'xi')  # subset of parameters interested in

# Create Stan data

# Define probit design matrix
probit_data <- data.frame(s = sim_samp$true_Si,
                          c = factor(analysis$c_i, levels=1:K))
V <- model.matrix(~ c * s, data = probit_data)
# Create array with assigned classes
V_k <- array(NA, dim=c(n, q, K))
for (k in 1:K) {
  temp_data <- data.frame(s = sim_samp$true_Si, 
                          c = factor(rep(k, n), levels=1:K))
  V_k[, ,k] <- model.matrix(~ c * s, data = temp_data) 
}

data_stan <- list(K = K, p = p, d = d, n = n, q = q, X = X, y = y, V = V, V_k = V_k)


# Create survey design
svy_data <- data.frame(s = sim_samp$true_Si, 
                       x = sim_samp$X_data,
                       y = sim_samp$Y_data, 
                       wts = sim_samp$sample_wt)
svydes <- svydesign(ids = ~1, strata = ~s, weights = ~wts, data = svy_data)
# create svrepdesign
svyrep <- as.svrepdesign(design = svydes, type = ctrl_rep$type, 
                         replicates = ctrl_rep$replicates)

# Run Stan model
# Stan will pass warnings from calling 0 chains, but will still create an 
# out_stan object for the 'grad_log_prob()' method
out_stan <- sampling(object = mod_stan, data = data_stan, pars = par_stan,
                     chains = 0, warmup = 0)
# 368 unconstrained parameters
# adjust_transform = FALSE
temp <- stan(file = paste0(analysis_dir, "mixture_model.stan"),
             data = data_stan, 
             iter = 100, chains = 4)

upars <- unconstrain_pars(out_stan, list("pi" = c(analysis$pi_med),
                                "theta" = analysis$theta_med,
                                "xi" = xi_ref))
grad_log_prob(out_stan, upars)

# Estimate Hessian
Hhat <- -1*optimHess(par_hat, 
                     gr = function(x){grad_log_prob(out_stan, x, 
                                                    adjust_transform = FALSE)})

# Estimate Jhat = Var(gradient)
print('gradient evaluation')

rep_temp <- withReplicates(design = svyrep, theta = grad_par, 
                           stanmod = mod_stan, standata = data_stan, 
                           par_stan = par_stan, par_hat = par_hat)
Jhat <- vcov(rep_temp)

# compute adjustment
Hi <- solve(Hhat)
V1 <- Hi %*% Jhat %*% Hi
R1 <- chol(V1)
R2i <- chol(Hi)
R2 <- solve(R2i)
R2R1 <- R2 %*% R1

# Adjust samples
par_adj <- aaply(par_samps, 1, DEadj, par_hat = par_hat, R2R1 = R2R1, 
                 .drop = FALSE)

output <- list(stan_fit = out_stan, sampled_params = par_samps, 
               adjusted_params = par_adj)


  
  


