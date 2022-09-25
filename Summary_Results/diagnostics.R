library(R.matlab)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(dplyr)
library(plyr)
library(tidyr)
library(rstan)
library(bayesplot)
library(coda)

setwd("C:/Users/Lang/Documents/Harvard/Research/Briana/supRPC/wsOFMM/Toy_Example")
# setwd("/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/Summary_Results")

#============================== Code for Graphics ======================================
# Read in data
scen <- 14
iter <- 2
samp_n <- 1

# Simulated data
pop_data <- readMat(paste0("simdata_scen2_iter", iter, ".mat"))$sim.data
samp_data <- readMat(paste0("simdata_scen14_iter", iter, "_samp", samp_n, ".mat"))$sim.data
names(pop_data) <- c("true_pi", "true_xi", "true_global_patterns", "true_global_thetas", "samp_ind", 
                     "sample_wt", "norm_const", "true_Si", "true_Ci", "X_data", "Y_data",
                     "true_Phi", "true_K")
names(samp_data) <- c(names(pop_data), "true_Li")

# Model output
data_ws <- readMat(paste0("wsOFMM_latent_results_scen", scen, "_iter", iter, 
                        "_samp", samp_n, ".mat"))
data_s <- readMat(paste0("sOFMM_latent_results_scen", scen, "_iter", iter, 
                       "_samp", samp_n, ".mat"))
data_ws <- data_ws$analysis
data_s <- data_s$analysis
names(data_ws) <- c("pi_med", "theta_med", "xi_med", "k_red", "pi_red", "theta_red", 
                    "xi_red", "c_i", "Phi_med", "z_i", "loglik_med", "dic6", "aebic")
names(data_s) <- c(names(data_ws), "posthoc_pi")


# Create dataframe to analyze MCMC chains

df_ws <- data.frame(data_ws$pi_red, data_ws$theta_red, data_ws$xi_red)
df_s <- data.frame(data_s$pi_red, data_s$theta_red, data_s$xi_red)

# Assuming both models have same number of components
pi_dim <- ncol(data_ws$pi_red)
theta_dim <- dim(data_ws$theta_red)
theta_names <- paste0(paste0(paste0("theta_", 1:theta_dim[2]), 
                             "_", rep(1:theta_dim[3], each=theta_dim[2])), 
                      "_", rep(1:theta_dim[4], each=(theta_dim[2]*theta_dim[3])))
xi_dim <- ncol(data_ws$xi_red)
colnames(df_ws) <- colnames(df_s) <- c(paste0("pi_", 1:pi_dim),
                                       theta_names,
                                       paste0("xi_", 1:xi_dim))
combined_df <- rbind(df_s, df_ws)
combined_df$Model <- c(rep("sOFMM", times=nrow(df_s)), rep("wsOFMM", times=nrow(df_ws)))
reshape_df <- combined_df %>% gather("pi", "pi_value", 1:pi_dim) %>%
  gather("theta", "theta_value", 1:length(theta_names)) %>%
  gather("xi", "xi_value", 1:xi_dim)

posthoc_df <- data.frame(paste0("pi_", 1:pi_dim), data_s$posthoc_pi)
colnames(posthoc_df) <- c("pi", "pi_value")
wsOFMM_pi_med <- data.frame(paste0("pi_", 1:pi_dim), t(data_ws$pi_med))
colnames(wsOFMM_pi_med) <- c("pi", "pi_value")
point_df <- rbind(posthoc_df, wsOFMM_pi_med)
point_df$Model <- c(rep("sOFMM", times=nrow(posthoc_df)), rep("wsOFMM", times=nrow(wsOFMM_pi_med)))

# Plot comparing pi for wsOFMM and sOFMM, with true values in green and sample values in red

pi_plot <- reshape_df %>% ggplot(aes(x=pi, y=pi_value, fill=Model)) + 
  geom_boxplot() +
  ggtitle(paste0("Parameter comparison for pi for iteration ", iter)) +
  theme(legend.position="top", plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept= (sum(pop_data[[9]]==1)/nrow(pop_data[[9]])), linetype="dashed", color="forestgreen") + 
  geom_hline(yintercept= (sum(pop_data[[9]]==2)/nrow(pop_data[[9]])), linetype="dashed", color="forestgreen") + 
  geom_hline(yintercept= (sum(samp_data[[9]]==1)/nrow(samp_data[[9]])), linetype="dashed", color="red") +
  geom_hline(yintercept= (sum(samp_data[[9]]==2)/nrow(samp_data[[9]])), linetype="dashed", color="red")
# pi_plot + geom_point(data=posthoc_df, aes(x=pi, y=pi_value), inherit.aes = FALSE, color="orange", size=3) +
#   geom_point(data=wsOFMM_pi_med, aes(x=pi, y=pi_value), inherit.aes = FALSE, color="blue", size=3)

pi_plot + geom_point(data=point_df, aes(x=pi, y=pi_value, color=Model), inherit.aes = FALSE, size=5)


# Plot comparing theta for wsOFMM and sOFMM, with true values in green

theta_labels <- paste0(paste0(paste0(1:theta_dim[2]), 
                             "_", rep(1:theta_dim[3], each=theta_dim[2])), 
                      "_", rep(1:theta_dim[4], each=(theta_dim[2]*theta_dim[3])))
reshape_df %>% ggplot(aes(x=theta, y=theta_value, fill=Model)) + 
  geom_boxplot() +
  ggtitle(paste0("Parameter comparison for theta for iteration 1")) +
  theme(legend.position="top", plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept= 0.9, linetype="dashed", color="forestgreen") + 
  geom_hline(yintercept= 0.1, linetype="dashed", color="forestgreen") +
  scale_x_discrete(labels=theta_labels) + 
  theme(axis.text.x=element_text(angle=40)) + 
  xlab("theta(j,k,r)")

# Plot comparing median xi for wsOFMM and sOFMM, with true values in green

reshape_df %>% ggplot(aes(xi, y=xi_value, fill=Model)) + geom_boxplot() +
  ggtitle(paste0("Parameter comparison for xi for iteration 1")) +
  theme(legend.position="top", plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept=pop_data[[2]][1], linetype="dashed", color="forestgreen") + 
  geom_hline(yintercept=pop_data[[2]][2], linetype="dashed", color="forestgreen") +
  geom_hline(yintercept=pop_data[[2]][3], linetype="dashed", color="forestgreen") +
  geom_hline(yintercept=pop_data[[2]][4], linetype="dashed", color="forestgreen") 

X_data <- samp_data[[10]]
S_data <- samp_data[[8]]
design_mat <- cbind(S_data, X_data)
X_data <- as.factor(X_data)
X_dummy <- dummy_cols(X_data)

reshape_df_ws <- df_ws %>% gather("pi", "pi_value", 1:pi_dim) %>%
  gather("theta", "theta_value", 1:length(theta_names)) %>%
  gather("xi", "xi_value", 1:xi_dim)
pi_summary <- ddply(reshape_df_ws, "pi", summarise, pi_median=median(pi_value), 
                    pi_lower=quantile(pi_value, 0.025), 
                    pi_upper=quantile(pi_value, 0.975))
theta_summary <- ddply(reshape_df_ws, "theta", summarise, theta_median=median(theta_value),
                       theta_lower=quantile(theta_value, 0.025),
                       theta_upper=quantile(theta_value, 0.975))
xi_summary <- ddply(reshape_df_ws, "xi", summarise, xi_median=median(xi_value),
                    xi_lower = quantile(xi_value, 0.025),
                    xi_upper = quantile(xi_value, 0.975))


#============================ Code for MCMC Diagnostics ========================

mcmc_hist(df_ws, pars=c("pi_1", "pi_2"))
mcmc_trace(df_ws, pars=paste0("pi_", 1:pi_dim))
mcmc_trace(df_ws, pars=theta_names)
mcmc_trace(df_ws) + ggtitle(paste0("Trace Plots for wsOFMM Iteration ", iter))

mcmc_trace(df_s) + ggtitle(paste0("Trace Plots for sOFMM Iteration ", iter))

ggplot(reshape_df_ws, aes(x=pi_value)) +
  geom_histogram(aes(y=..density..))+ 
  facet_wrap(~pi) +
  geom_vline(aes(xintercept=pi_median), data=pi_summary, color="red")

#============================ MSE Summaries across Iterations ==================

# Get summary statistics
m <- nrow(data_ws$pi_red)
summaries <- as.data.frame(matrix(NA, nrow=m, ncol=4))
colnames(summaries) <- c("MSE_pi", "MSE_theta", "MSE_xi", "MSE_Phi")
K_match <- ncol(data_ws$pi_red) == pop_data$true_K


agg_results <- function(scenarios, samp_n) {
  len <- length(scenarios)
  output <- NULL
  Model <- NULL
  Scenario <- NULL
  
  for (scen in scenarios) {
    if (scen > 4) {
      results_s <- readMat(paste0("summary_sOFMM_scen", scen, "_samp", samp_n, ".mat"))
      results_ws <- readMat(paste0("summary_wsOFMM_scen", scen, "_samp", samp_n, ".mat"))
    } else {
      results_s <- readMat(paste0("summary_sOFMM_scen", scen, "_pop", ".mat"))
      results_ws <- readMat(paste0("summary_wsOFMM_scen", scen, "_pop", ".mat"))
    }
    res_s <- results_s$all
    res_ws <- results_ws$all
    res_s <- do.call(cbind.data.frame, results_s$all)
    res_ws <- do.call(cbind.data.frame, results_ws$all)
    colnames(res_s) <- colnames(res_ws) <- c("K_red", "K_match", "Sens", "Spec", "DIC6", "AEBIC", "Phi_MSE", 
                                             "Theta_MSE", "Corr(S,C)", "Corr(S,Y)", "Corr(C,Y)")
    output <- rbind(output, rbind(res_s, res_ws))
    Model <- c(Model, c(rep("sOFMM", times=nrow(res_s)), rep("wsOFMM", times=nrow(res_ws))))
    Scenario <- c(Scenario, rep(scen, times=nrow(res_s)+nrow(res_ws)))
  }
  
  output$Model <- Model
  output$Scenario <- Scenario
  
  return(output)
}

scenarios <- seq(5, 16, by=1)
output <- agg_results(scenarios, 1)


pdf(file="plots_scen_5_16.pdf")

output %>% ggplot(aes(x=Model, y=K_red, fill=Model)) + 
  geom_boxplot() +
  ggtitle(paste0("Model Performance Comparison for K_red Across Scenarios 5-16")) +
  theme(legend.position="top", plot.title = element_text(hjust = 0.5)) +
  stat_summary(aes(label=round(..y..,2)), fun = median, geom="text", 
               position=position_nudge(y=2), size=4) + 
  facet_wrap(~Scenario)
  
output %>% ggplot(aes(x=Model, y=Sens, fill=Model)) + 
  geom_boxplot() +
  ggtitle(paste0("Model Performance Comparison for Sensitivity Across Scenarios 5-16")) +
  theme(legend.position="top", plot.title = element_text(hjust = 0.5)) +
  stat_summary(aes(label=round(..y..,5)), fun = median, geom="text", 
               position=position_nudge(y=0.1), size=4) + 
  facet_wrap(~Scenario)

output %>% ggplot(aes(x=Model, y=Spec, fill=Model)) + 
  geom_boxplot() +
  ggtitle(paste0("Model Performance Comparison for Specificity Across Scenarios 5-16")) +
  theme(legend.position="top", plot.title = element_text(hjust = 0.5)) +
  stat_summary(aes(label=round(..y..,5)), fun = median, geom="text", 
               position=position_nudge(y=0.01), size=4) + 
  facet_wrap(~Scenario)

output %>% ggplot(aes(x=Model, y=Phi_MSE, fill=Model)) + 
  geom_boxplot() +
  ggtitle(paste0("Model Performance Comparison for Phi_MSE Across Scenarios 5-16")) +
  theme(legend.position="top", plot.title = element_text(hjust = 0.5)) +
  stat_summary(aes(label=round(..y..,5)), fun = median, geom="text", 
               position=position_nudge(y=0.03), size=4) + 
  facet_wrap(~Scenario)

output %>% ggplot(aes(x=Model, y=Theta_MSE, fill=Model)) + 
  geom_boxplot() +
  ggtitle(paste0("Model Performance Comparison for Theta_MSE Across Scenarios 5-16")) +
  theme(legend.position="top", plot.title = element_text(hjust = 0.5)) +
  stat_summary(aes(label=round(..y..,5)), fun = median, geom="text", 
               position=position_nudge(y=0.06), size=4) + 
  facet_wrap(~Scenario)

output %>% ggplot(aes(x=Model, y=DIC6, fill=Model)) + 
  geom_boxplot() +
  ggtitle(paste0("Model Performance Comparison for DIC6 Across Scenarios 5-16")) +
  theme(legend.position="top", plot.title = element_text(hjust = 0.5)) +
  stat_summary(aes(label=round(..y..,0)), fun = median, geom="text", 
               position=position_nudge(y=7000), size=4) + 
  facet_wrap(~Scenario)

output %>% ggplot(aes(x=Model, y=AEBIC, fill=Model)) + 
  geom_boxplot() +
  ggtitle(paste0("Model Performance Comparison for AEBIC Across Scenarios 5-16")) +
  theme(legend.position="top", plot.title = element_text(hjust = 0.5)) +
  stat_summary(aes(label=round(..y..,0)), fun = median, geom="text", 
               position=position_nudge(y=5000), size=4) + 
  facet_wrap(~Scenario)

dev.off()

# output %>% ggviolin(x="Model", y="Phi_MSE", fill="Model", add="boxplot") + 
#   facet_wrap(~Scenario)


# output %>% ggviolin(x="Scenario", y="Phi_MSE", fill="Model", add="boxplot")

# mcmc_scen3_iter68 <- readMat("wsOFMM_MCMC_scen3_iter68.mat")
# data <- mcmc_scen3_iter68$MCMC.out
# names(data) <- c("pi", "theta", "c_i", "xi", "z_i", "loglik", "runtime")
# data_mcmc <- mcmcUpgrade(mcmc(data))
# coda::traceplot(data_mcmc)
