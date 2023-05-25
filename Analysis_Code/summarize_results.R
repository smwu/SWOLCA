#==================================
## Summarize and Plot Model Results
## Programmer: SM Wu   
## Last Updated: 2023/05/20
#==================================

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

source("Analysis_Code/summarize_results_functions.R")

#================ Summarize and save results ===================================

scen_pop <- 1112
scen_samp <- 111211
WSOLCA_name <- "_results_2mod20000_adjRcpp_scen"
SOLCA_name <- "_results_2mod20000_scen"
WOLCA_name <- "_results_wt_2mod20000_scen"

# Baseline: stratified, 5% sample size, 85/5/5/5, non-supervised, confounder
save_scen_metrics(scen_pop = 1112, scen_samp = 111211, WSOLCA = TRUE, 
                  SOLCA = TRUE, WOLCA = TRUE, WSOLCA_name = WSOLCA_name, 
                  SOLCA_name = SOLCA_name, WOLCA_name = WOLCA_name,
                  save_name = "metrics_scen")
# Sampling schemes: stratified cluster
save_scen_metrics(scen_pop = 1112, scen_samp = 111212, WSOLCA = TRUE, 
                  SOLCA = TRUE, WOLCA = TRUE, WSOLCA_name = WSOLCA_name, 
                  SOLCA_name = SOLCA_name, WOLCA_name = WOLCA_name,
                  save_name = "metrics_scen")
# Sampling schemes: SRS
save_scen_metrics(scen_pop = 1112, scen_samp = 111213, WSOLCA = TRUE, 
                  SOLCA = TRUE, WOLCA = TRUE, WSOLCA_name = WSOLCA_name, 
                  SOLCA_name = SOLCA_name, WOLCA_name = WOLCA_name,
                  save_name = "metrics_scen")
# Sample size: 10%
save_scen_metrics(scen_pop = 1112, scen_samp = 111221, WSOLCA = TRUE, 
                  SOLCA = TRUE, WOLCA = TRUE, WSOLCA_name = WSOLCA_name, 
                  SOLCA_name = SOLCA_name, WOLCA_name = WOLCA_name,
                  save_name = "metrics_scen")
# Sample size: 1%
save_scen_metrics(scen_pop = 1112, scen_samp = 111231, WSOLCA = TRUE, 
                  SOLCA = TRUE, WOLCA = TRUE, WSOLCA_name = WSOLCA_name, 
                  SOLCA_name = SOLCA_name, WOLCA_name = WOLCA_name,
                  save_name = "metrics_scen")
# Weak patterns
save_scen_metrics(scen_pop = 2112, scen_samp = 211211, WSOLCA = TRUE, 
                  SOLCA = TRUE, WOLCA = TRUE, WSOLCA_name = WSOLCA_name, 
                  SOLCA_name = SOLCA_name, WOLCA_name = WOLCA_name,
                  save_name = "metrics_scen")
# Supervised
save_scen_metrics(scen_pop = 1212, scen_samp = 121211, WSOLCA = TRUE, 
                  SOLCA = TRUE, WOLCA = TRUE, WSOLCA_name = WSOLCA_name, 
                  SOLCA_name = SOLCA_name, WOLCA_name = WOLCA_name,
                  save_name = "metrics_scen")
# Effect modifier
save_scen_metrics(scen_pop = 1122, scen_samp = 112211, WSOLCA = TRUE, 
                  SOLCA = TRUE, WOLCA = TRUE, WSOLCA_name = WSOLCA_name, 
                  SOLCA_name = SOLCA_name, WOLCA_name = WOLCA_name,
                  save_name = "metrics_scen")
# Effect modifier marginal
save_scen_metrics(scen_pop = 1122, scen_samp = 112211, WSOLCA = TRUE, 
                  SOLCA = TRUE, WOLCA = TRUE, 
                  WSOLCA_name = "_results_effmod_adjRcpp_scen", 
                  SOLCA_name = "_results_effmod_scen", 
                  WOLCA_name = "_results_wt_effmod_scen",
                  marg = TRUE, 
                  save_name = "metrics_marg_scen")

# Load simulated population data
load(paste0(wd, data_dir, "simdata_scen", scen_pop,"_iter", iter_pop, ".RData"))

# Obtain true observed population parameters
true_params <- get_true_params(sim_pop = sim_pop)  

#================ TABLE METRICS SUMMARY ========================================

wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/" # Working directory
data_dir <- "Data/"               # Simulated data directory
res_dir <- "Results/"             # Model results directory
analysis_dir <- "Analysis_Code/"  # Analysis directory where metrics are saved

create_table1(wd = wd, analysis_dir = analysis_dir)


#================ PLOT MSE SUMMARY =========================================

# Sampling scenarios
scenarios <- c(111213, 111211, 111212)
scen_names <- c("SRS", "Stratified", "Stratified Cluster")
save_names <- rep("metrics_scen", 3)
samp_pi <- plot_rmse_boxplot(wd = wd, analysis_dir = analysis_dir, 
                          save_names = save_names, scenarios = scenarios,
                          scen_names = scen_names, overall_name = "Sampling", 
                          param = "pi", upper_lim = 0.75, xlab = "Sampling Scheme", 
                          ylab = expression("RMSE for "~pi))
samp_theta <- plot_rmse_boxplot(wd = wd, analysis_dir = analysis_dir, 
                             save_names = save_names, scenarios = scenarios,
                             scen_names = scen_names, overall_name = "Sampling", 
                             param = "theta", upper_lim = 0.75, xlab = "Sampling Scheme", 
                             ylab = expression("RMSE for "~theta))
samp_xi <- plot_rmse_boxplot(wd = wd, analysis_dir = analysis_dir, 
                          save_names = save_names, scenarios = scenarios,
                          scen_names = scen_names, overall_name = "Sampling", 
                          param = "xi", upper_lim = 0.75, xlab = "Sampling Scheme", 
                          ylab = expression("RMSE for "~xi))
ggarrange(samp_pi, samp_theta, samp_xi, nrow = 3, ncol = 1, common.legend = TRUE)
ggarrange(samp_pi, samp_theta, samp_xi, nrow = 1, ncol = 3, common.legend = TRUE)

# Pattern strength and separation
scenarios <- c(111211, 211211, 121211)
scen_names <- c("Mode 85%", "Mode 55%", "Overlapping")
pattern_pi <- plot_mse_boxplot(wd = wd, analysis_dir = analysis_dir, 
                            save_names = save_names, scenarios = scenarios,
                            scen_names = scen_names, overall_name = "Pattern", 
                            param = "pi", upper_lim = 0.1, 
                            xlab = "Pattern Strength and Separation", 
                            ylab = expression("MSE for "~pi))
pattern_theta <- plot_mse_boxplot(wd = wd, analysis_dir = analysis_dir, 
                                  save_names = save_names, scenarios = scenarios,
                               scen_names = scen_names, overall_name = "Pattern", 
                               param = "theta", upper_lim = 0.1, 
                               xlab = "Pattern Strength and Separation", 
                               ylab = expression("MSE for "~theta))
pattern_xi <- plot_mse_boxplot(wd = wd, analysis_dir = analysis_dir, 
                               save_names = save_names, scenarios = scenarios,
                            scen_names = scen_names, overall_name = "Pattern", 
                            param = "xi", upper_lim = 0.15, 
                            xlab = "Pattern Strength and Separation", 
                            ylab = expression("MSE for "~xi))
ggarrange(pattern_pi, pattern_theta, pattern_xi, nrow = 1, ncol = 3, common.legend = TRUE)

# Sample size
scenarios <- c(111221, 111211, 111231)
scen_names <- c("10%", "5%", "1%")
ss_pi <- plot_mse_boxplot(wd = wd, analysis_dir = analysis_dir, 
                          save_names = save_names, scenarios = scenarios,
                               scen_names = scen_names, overall_name = "SS", 
                               param = "pi", upper_lim = 0.1, 
                               xlab = "Sample Size", 
                               ylab = expression("MSE for "~pi))
ss_theta <- plot_mse_boxplot(wd = wd, analysis_dir = analysis_dir, 
                             save_names = save_names, scenarios = scenarios,
                                  scen_names = scen_names, overall_name = "SS", 
                                  param = "theta", upper_lim = 0.1, 
                                  xlab = "Sample Size", 
                                  ylab = expression("MSE for "~theta))
ss_xi <- plot_mse_boxplot(wd = wd, analysis_dir = analysis_dir, 
                          save_names = save_names, scenarios = scenarios,
                               scen_names = scen_names, overall_name = "SS", 
                               param = "xi", upper_lim = 0.8, 
                               xlab = "Sample Size", 
                               ylab = expression("MSE for "~xi))
ggarrange(ss_pi, ss_theta, ss_xi, nrow = 1, ncol = 3, common.legend = TRUE)


# Selection bias
scenarios <- c(111211, 112211, 112211)
scen_names <- c("Confounder", "Precision \nMeasured", "Precision \nUnmeasured")
save_names <- c("metrics_scen", "metrics_scen", "metrics_marg_scen")
selection_pi <- plot_mse_boxplot(wd = wd, analysis_dir = analysis_dir, 
                          save_names = save_names, scenarios = scenarios,
                          scen_names = scen_names, overall_name = "selection", 
                          param = "pi", upper_lim = 0.025, 
                          xlab = "Selection Bias", 
                          ylab = expression("MSE for "~pi))
selection_theta <- plot_mse_boxplot(wd = wd, analysis_dir = analysis_dir, 
                             save_names = save_names, scenarios = scenarios,
                             scen_names = scen_names, overall_name = "selection", 
                             param = "theta", upper_lim = 0.05, 
                             xlab = "Selection Bias", 
                             ylab = expression("MSE for "~theta))
selection_xi <- plot_mse_boxplot(wd = wd, analysis_dir = analysis_dir, 
                          save_names = save_names, scenarios = scenarios,
                          scen_names = scen_names, overall_name = "selection", 
                          param = "xi", upper_lim = 0.2, 
                          xlab = "Selection Bias", 
                          ylab = expression("MSE for "~xi))
ggarrange(selection_pi, selection_theta, selection_xi, nrow = 1, ncol = 3, 
          common.legend = TRUE)

#================ PLOT BIAS SUMMARY =========================================

scenarios <- c(111213, 111211, 111212)
scen_names <- c("SRS", "Stratified", "Stratified Cluster")
p_pi <- plot_bias_boxplot(wd = wd, analysis_dir = analysis_dir, 
                  save_names = save_names, scenarios = scenarios,
                  scen_names = scen_names, overall_name = "Sampling", 
                  param = "pi", upper_lim = 0.15, xlab = "Sampling Scheme", 
                  ylab = expression("Mean Absolute Bias for "~pi))
p_theta <- plot_bias_boxplot(wd = wd, analysis_dir = analysis_dir, 
                  save_names = save_names, scenarios = scenarios,
                  scen_names = scen_names, overall_name = "Sampling", 
                  param = "theta", upper_lim = 0.1, xlab = "Sampling Scheme", 
                  ylab = expression("Mean Absolute Bias for "~theta))
p_xi <- plot_bias_boxplot(wd = wd, analysis_dir = analysis_dir, 
                  save_names = save_names, scenarios = scenarios,
                  scen_names = scen_names, overall_name = "Sampling", 
                  param = "xi", upper_lim = 0.4, xlab = "Sampling Scheme", 
                  ylab = expression("Mean Absolute Bias for "~xi))
ggarrange(p_pi, p_theta, p_xi, nrow = 3, ncol = 1, common.legend = TRUE)
ggarrange(p_pi, p_theta, p_xi, nrow = 1, ncol = 3, common.legend = TRUE)

### Plot bias as grouped boxplot 
plot_bias_boxplot <- function(wd, analysis_dir, save_names, scenarios, 
                              scen_names, overall_name, param,
                              lower_lim = 0, upper_lim = 1, xlab, ylab) {
  if (param == "pi") {
    param_dist <- "pi_dist"
  } else if (param == "theta") {
    param_dist <- "theta_dist"
  } else if (param == "xi") {
    param_dist <- "xi_dist"
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
  # Dimensions (9*L)x3 with columns Model, Sampling Scheme, Mean Absolute Bias
  plot_df <- plot_df %>% 
    gather({{overall_name}}, "Bias", -Model)
  
  # Create grouped boxplot of pi values for the sampling schemes
  p <- plot_df %>% 
    ggplot(aes_string(x = overall_name, y = "Bias", fill = "Model")) + 
    theme_bw() + 
    geom_boxplot() + 
    ylim(lower_lim, upper_lim) + 
    xlab({{xlab}}) + ylab({{ylab}})
  
  return(p)
}


#============== Create Appendix Tables =========================================
wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/" # Working directory
analysis_dir <- "Analysis_Code/"  # Analysis directory where metrics are saved

# Sampling scenarios
scenarios <- c(111213, 111211, 111212)
scen_names <- c("SRS", "Stratified", "Stratified Cluster")
save_names <- rep("metrics_scen", 3)
create_app_tables(wd, analysis_dir, save_names, scenarios, scen_names, 
                  overall_name = "Sampling Scheme", format = "html")

# Pattern strength and separation
scenarios <- c(111211, 211211, 121211)
scen_names <- c("Mode 85%", "Mode 55%", "Overlapping")
save_names <- rep("metrics_scen", 3)
create_app_tables(wd, analysis_dir, save_names, scenarios, scen_names, 
                  overall_name = "Pattern", format = "html")
  
# Sample size
scenarios <- c(111221, 111211, 111231)
scen_names <- c("10%", "5%", "1%")
save_names <- rep("metrics_scen", 3)
create_app_tables(wd, analysis_dir, save_names, scenarios, scen_names, 
                  overall_name = "Sample Size", format = "html")

# Selection bias
scenarios <- c(111211, 112211, 112211)
scen_names <- c("Confounder", "Precision \nMeasured", "Precision \nUnmeasured")
save_names <- c("metrics_scen", "metrics_scen", "metrics_marg_scen")
create_app_tables(wd, analysis_dir, save_names, scenarios, scen_names, 
                  overall_name = "Selection Bias", format = "html")

#================== OLD CODE ===================================================
load(paste0(wd, analysis_dir, "metrics_scen", 111211, ".RData"))


create_table(metrics_s = metrics_all$metrics_s, 
             metrics_ws = metrics_all$metrics_ws, 
             metrics_unsup = metrics_all$metrics_unsup,
             scen_samp = scen_samp)

### Create table of metrics with bias and variance
# Inputs:
#   metrics_s: Summary metrics for SOLCA
#   metrics_ws: Summary metrics for WSOLCA
#   metrics_unsup: Summary metrics for WOLCA
#   scen_samp: Sampling scenario
# Output: Formatted table with absolute bias, CI width, and coverage
create_table <- function(metrics_s, metrics_ws, metrics_unsup, scen_samp) {
  metrics_summ <- as.data.frame(matrix(NA, nrow=3, ncol=12))
  colnames(metrics_summ) <- c("Sampling Scheme", "Model", 
                              "K Bias^2", "$\\pi$ Bias^2", "$\\pi$ CI width", 
                              "$\\theta$ Bias^2", "$\\theta$ CI width", 
                              "$\\xi$ Bias^2", "$\\xi$ CI width", 
                              "$\\pi$ Coverage","$\\theta$ Coverage", "$\\xi$ Coverage")
  metrics_summ[, 1] <- rep("Stratified", 3)
  metrics_summ[, 2] <- rep(c("Unwtd(sOFMM)", "Wtd(wsOFMM)", "Unsup(wOFMM)"), 1)  ## latent versions
  output_inds <- 1:7
  metrics_summ[1, -c(1,2)] <- c(metrics_s[output_inds], 
                                mean(metrics_s$pi_cover_avg), 
                                mean(metrics_s$theta_cover_avg),
                                mean(metrics_s$xi_cover_avg))
  metrics_summ[2, -c(1,2)] <- c(metrics_ws[output_inds], 
                                mean(metrics_ws$pi_cover_avg), 
                                mean(metrics_ws$theta_cover_avg),
                                mean(metrics_ws$xi_cover_avg))
  metrics_summ[3, -c(1,2)] <- c(metrics_unsup[output_inds], 
                                mean(metrics_unsup$pi_cover_avg), 
                                mean(metrics_unsup$theta_cover_avg),
                                mean(metrics_unsup$xi_cover_avg))
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
}



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


