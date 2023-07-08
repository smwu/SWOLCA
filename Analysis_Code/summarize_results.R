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
WSOLCA_name <- "_results_comb_adjRcpp_scen"
SOLCA_name <- "_results_scen"
WOLCA_name <- "_results_wt_scen"

# Baseline: stratified, 5% sample size, 85/5/5/5, non-supervised, confounder
save_scen_metrics(scen_pop = 1112, scen_samp = 111211, WSOLCA = TRUE, 
                  SOLCA = TRUE, WOLCA = TRUE, WSOLCA_name = WSOLCA_name, 
                  SOLCA_name = SOLCA_name, WOLCA_name = WOLCA_name,
                  covs = "true_Si", save_name = "metrics_scen")
# Sampling schemes: stratified cluster
save_scen_metrics(scen_pop = 1112, scen_samp = 111212, WSOLCA = TRUE, 
                  SOLCA = TRUE, WOLCA = TRUE, WSOLCA_name = WSOLCA_name, 
                  SOLCA_name = SOLCA_name, WOLCA_name = WOLCA_name,
                  covs = "true_Si", save_name = "metrics_scen")
# Sampling schemes: SRS
save_scen_metrics(scen_pop = 1112, scen_samp = 111213, WSOLCA = TRUE, 
                  SOLCA = TRUE, WOLCA = TRUE, WSOLCA_name = WSOLCA_name, 
                  SOLCA_name = SOLCA_name, WOLCA_name = WOLCA_name,
                  covs = "true_Si", save_name = "metrics_scen")
# Sample size: 10%
save_scen_metrics(scen_pop = 1112, scen_samp = 111221, WSOLCA = TRUE, 
                  SOLCA = TRUE, WOLCA = TRUE, WSOLCA_name = WSOLCA_name, 
                  SOLCA_name = SOLCA_name, WOLCA_name = WOLCA_name,
                  covs = "true_Si", save_name = "metrics_scen")
# Sample size: 1%
save_scen_metrics(scen_pop = 1112, scen_samp = 111231, WSOLCA = TRUE, 
                  SOLCA = TRUE, WOLCA = TRUE, WSOLCA_name = WSOLCA_name, 
                  SOLCA_name = SOLCA_name, WOLCA_name = WOLCA_name,
                  covs = "true_Si", save_name = "metrics_scen")
# Weak patterns
save_scen_metrics(scen_pop = 2112, scen_samp = 211211, WSOLCA = TRUE, 
                  SOLCA = TRUE, WOLCA = TRUE, WSOLCA_name = WSOLCA_name, 
                  SOLCA_name = SOLCA_name, WOLCA_name = WOLCA_name,
                  covs = "true_Si", save_name = "metrics_scen")
# Supervised
save_scen_metrics(scen_pop = 1212, scen_samp = 121211, WSOLCA = TRUE, 
                  SOLCA = TRUE, WOLCA = TRUE, WSOLCA_name = WSOLCA_name, 
                  SOLCA_name = SOLCA_name, WOLCA_name = WOLCA_name,
                  covs = "true_Si", save_name = "metrics_scen")
# Effect modifier
save_scen_metrics(scen_pop = 1122, scen_samp = 112211, WSOLCA = TRUE, 
                  SOLCA = TRUE, WOLCA = TRUE, WSOLCA_name = WSOLCA_name, 
                  SOLCA_name = SOLCA_name, WOLCA_name = WOLCA_name,
                  covs = "true_Si", save_name = "metrics_scen")
# Effect modifier marginal
save_scen_metrics(scen_pop = 1122, scen_samp = 112211, WSOLCA = TRUE, 
                  SOLCA = TRUE, WOLCA = TRUE, 
                  WSOLCA_name = WSOLCA_name, SOLCA_name = SOLCA_name, 
                  WOLCA_name = WOLCA_name,
                  # WSOLCA_name = "_results_effmod_adjRcpp_scen", 
                  # SOLCA_name = "_results_effmod_scen", 
                  # WOLCA_name = "_results_wt_effmod_scen",
                  covs = NULL, save_name = "metrics_marg_scen")
# Confounder marginal
save_scen_metrics(scen_pop = 1112, scen_samp = 111211, WSOLCA = TRUE, 
                  SOLCA = TRUE, WOLCA = TRUE, 
                  WSOLCA_name = "_results_effmod_adjRcpp_scen", 
                  SOLCA_name = "_results_effmod_scen", 
                  WOLCA_name = "_results_wt_effmod_scen",
                  covs = NULL, save_name = "metrics_marg_scen")
# Additional confounders
save_scen_metrics(scen_pop = 1132, scen_samp = 113211, WSOLCA = TRUE, 
                  SOLCA = TRUE, WOLCA = TRUE, WSOLCA_name = WSOLCA_name, 
                  SOLCA_name = SOLCA_name, WOLCA_name = WOLCA_name,
                  covs = "additional", save_name = "metrics_scen")

#================ TABLE METRICS SUMMARY ========================================

wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/" # Working directory
data_dir <- "Data/June22/"               # Simulated data directory
res_dir <- "Results/July3/"             # Model results directory
analysis_dir <- "Analysis_Code/"  # Analysis directory where metrics are saved

create_table1(wd = wd, analysis_dir = analysis_dir, format = "latex")


#================ PLOT RMSE SUMMARY =========================================

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
pattern_pi <- plot_rmse_boxplot(wd = wd, analysis_dir = analysis_dir, 
                            save_names = save_names, scenarios = scenarios,
                            scen_names = scen_names, overall_name = "Pattern", 
                            param = "pi", upper_lim = 0.4, 
                            xlab = "Pattern Strength and Separation", 
                            ylab = expression("RMSE for "~pi))
pattern_theta <- plot_rmse_boxplot(wd = wd, analysis_dir = analysis_dir, 
                                  save_names = save_names, scenarios = scenarios,
                               scen_names = scen_names, overall_name = "Pattern", 
                               param = "theta", upper_lim = 0.4, 
                               xlab = "Pattern Strength and Separation", 
                               ylab = expression("RMSE for "~theta))
pattern_xi <- plot_rmse_boxplot(wd = wd, analysis_dir = analysis_dir, 
                               save_names = save_names, scenarios = scenarios,
                            scen_names = scen_names, overall_name = "Pattern", 
                            param = "xi", upper_lim = 0.75, 
                            xlab = "Pattern Strength and Separation", 
                            ylab = expression("RMSE for "~xi))
ggarrange(pattern_pi, pattern_theta, pattern_xi, nrow = 1, ncol = 3, common.legend = TRUE)

# Sample size
scenarios <- c(111221, 111211, 111231)
scen_names <- c("10%", "5%", "1%")
ss_pi <- plot_rmse_boxplot(wd = wd, analysis_dir = analysis_dir, 
                          save_names = save_names, scenarios = scenarios,
                               scen_names = scen_names, overall_name = "SS", 
                               param = "pi", upper_lim = 0.75, 
                               xlab = "Sample Size", 
                               ylab = expression("RMSE for "~pi))
ss_theta <- plot_rmse_boxplot(wd = wd, analysis_dir = analysis_dir, 
                             save_names = save_names, scenarios = scenarios,
                                  scen_names = scen_names, overall_name = "SS", 
                                  param = "theta", upper_lim = 0.75, 
                                  xlab = "Sample Size", 
                                  ylab = expression("RMSE for "~theta))
ss_xi <- plot_rmse_boxplot(wd = wd, analysis_dir = analysis_dir, 
                          save_names = save_names, scenarios = scenarios,
                               scen_names = scen_names, overall_name = "SS", 
                               param = "xi", upper_lim = 0.75, 
                               xlab = "Sample Size", 
                               ylab = expression("RMSE for "~xi))
ggarrange(ss_pi, ss_theta, ss_xi, nrow = 1, ncol = 3, common.legend = TRUE)


# Selection bias
scenarios <- c(111211, 112211, 113211)
scen_names <- c("Confounder \nMeasured", "Precision \nUnmeasured", 
                "Additional \nConfounders")
save_names <- c("metrics_scen", "metrics_marg_scen", "metrics_scen")
selection_pi <- plot_rmse_boxplot(wd = wd, analysis_dir = analysis_dir, 
                          save_names = save_names, scenarios = scenarios,
                          scen_names = scen_names, overall_name = "selection", 
                          param = "pi", upper_lim = 0.5, 
                          xlab = "Selection Bias", 
                          ylab = expression("RMSE for "~pi))
selection_theta <- plot_rmse_boxplot(wd = wd, analysis_dir = analysis_dir, 
                             save_names = save_names, scenarios = scenarios,
                             scen_names = scen_names, overall_name = "selection", 
                             param = "theta", upper_lim = 0.5, 
                             xlab = "Selection Bias", 
                             ylab = expression("RMSE for "~theta))
selection_xi <- plot_rmse_boxplot(wd = wd, analysis_dir = analysis_dir, 
                          save_names = save_names, scenarios = scenarios,
                          scen_names = scen_names, overall_name = "selection", 
                          param = "xi", upper_lim = 0.75, 
                          xlab = "Selection Bias", 
                          ylab = expression("RMSE for "~xi))
ggarrange(selection_pi, selection_theta, selection_xi, nrow = 1, ncol = 3, 
          common.legend = TRUE)

# Examining xi across scenarios
ggarrange(samp_xi, pattern_xi, ss_xi, selection_xi, nrow = 2, ncol = 2,
         common.legend = TRUE, legend = "top")

samp_xi


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
scen_names <- c("10% (n=8000)", "5% (n=4000)", "1% (n=800)")
save_names <- rep("metrics_scen", 3)
create_app_tables(wd, analysis_dir, save_names, scenarios, scen_names, 
                  overall_name = "Sample Size", format = "html")

# Selection bias
scenarios <- c(111211, 112211, 113211)
scen_names <- c("Confounder \nMeasured", "Precision \nUnmeasured", 
                "Additional \nConfounders")
save_names <- c("metrics_scen", "metrics_marg_scen", "metrics_scen")
create_app_tables(wd, analysis_dir, save_names, scenarios, scen_names, 
                  overall_name = "Selection Bias", format = "html")


#================= Plot theta patterns =========================================
# Load simulated population data

# Stratified sampling
p_theta_strat <- plot_theta_patterns(wd = wd, data_dir = data_dir, 
                                     analysis_dir = analysis_dir,
                    scen_pop = 1112, iter_pop = 1, scen_samp = 111211,
                    scen_name = "metrics_scen") 
annotate_figure(p_theta_strat, 
  top = text_grob("Modal pattern per latent class, averaged across 100 iterations"))

# Stratified cluster sampling
p_theta_clus <- plot_theta_patterns(wd = wd, data_dir = data_dir, 
                                    analysis_dir = analysis_dir,
                    scen_pop = 1112, iter_pop = 1, scen_samp = 111212,
                    scen_name = "metrics_scen")

# Investigate wrong pattern for WOLCA
# Load summary result
load(paste0(wd, analysis_dir, "metrics_scen", 111211, ".RData"))
K_wrong <- which(metrics_all$metrics_unsup$K_dist != 0)
sim_res_path_full <- paste0(wd, res_dir, "wOFMM", "_results_wt_2mod20000_scen", 
                            111211, "_samp", K_wrong, ".RData")
load(sim_res_path_full)
dim(res$analysis$theta_med)  # Only 2 classes
theta_plot_data <- apply(res$analysis$theta_med, c(1,2), which.max)
theta_mode_plot(theta_plot_data, "WOLCA Classes")



#========= Plot pi over all iterations, for sampling designs ===================

p_pi_samp <- plot_pi_samp(wd = wd, analysis_dir = analysis_dir, scen_pop = 1112)
plot_xi_samp(wd = wd, analysis_dir = analysis_dir, scen_pop = 1112)
ggarrange(p_pi_samp, 
          p_theta_strat,
          samp_xi + theme(legend.position = "top"), 
          ncol = 3, widths = c(1,1,0.8))

#================ COVERAGE PLOTS ===============================================
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
plot_cov <- data.frame(
  coverage = c(mean(metrics_SRS_s$pi_cover_avg), 
               mean(metrics_SRS_unsup$pi_cover_avg),
               mean(metrics_SRS_ws$pi_cover_avg),
               mean(metrics_Strat_s$pi_cover_avg), 
               mean(metrics_Strat_unsup$pi_cover_avg),
               mean(metrics_Strat_ws$pi_cover_avg),
               mean(metrics_Clus_s$pi_cover_avg), 
               mean(metrics_Clus_unsup$xi_cover_avg),
               mean(metrics_Clus_ws$pi_cover_avg),
               mean(metrics_SRS_s$theta_cover_avg), 
               mean(metrics_SRS_unsup$theta_cover_avg),
               mean(metrics_SRS_ws$theta_cover_avg),
               mean(metrics_Strat_s$theta_cover_avg), 
               mean(metrics_Strat_unsup$theta_cover_avg),
               mean(metrics_Strat_ws$theta_cover_avg),
               mean(metrics_Clus_s$theta_cover_avg), 
               mean(metrics_Clus_unsup$theta_cover_avg),
               mean(metrics_Clus_ws$theta_cover_avg),
               mean(metrics_SRS_s$xi_cover_avg), 
               mean(metrics_SRS_unsup$xi_cover_avg),
               mean(metrics_SRS_ws$xi_cover_avg),
               mean(metrics_Strat_s$xi_cover_avg), 
               mean(metrics_Strat_unsup$xi_cover_avg),
               mean(metrics_Strat_ws$xi_cover_avg),
               mean(metrics_Clus_s$xi_cover_avg), 
               mean(metrics_Clus_unsup$xi_cover_avg),
               mean(metrics_Clus_ws$xi_cover_avg)),
  param = factor(c(rep("π", 9), rep("θ", 9), rep("ξ", 9)), 
                 levels = c("π","θ","ξ")),
  model = rep(c("SOLCA", "WOLCA", "WSOLCA"), times=9),
  sampling = rep(rep(c("SRS", "Stratified", "Stratified\n Cluster"), each = 3), 
                 times=3)
)
plot_cov %>%  ggplot(aes(x=param, y=coverage, color=model, shape=sampling)) +
  theme_bw() + scale_color_brewer(palette="Set2") + 
  labs(x = "Parameter", y = "Coverage", color = "Model", shape = "Sampling \nScheme") +
  geom_point(size=3, position = position_dodge(0.6)) +
  geom_hline(aes(yintercept=0.95)) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))

plot_cov %>%  ggplot(aes(x=param, y=coverage, fill=model)) +
  theme_bw() + scale_fill_brewer(palette="RdYlBu") + 
  labs(x = "Parameter", y = "Coverage", fill = "Model") +
  geom_point(size=2.5, position = position_dodge(0.3), shape=21) +
  geom_hline(aes(yintercept=0.95)) + ylim(0.5,1) + 
  facet_grid(~sampling)


#================= Plot pi over 100 iterations ==============================

# Stratified sampling
p1 <- plot_pi_patterns(wd = wd, data_dir = data_dir, analysis_dir = analysis_dir,
                       scen_pop = 1112, iter_pop = 1, scen_samp = 111211,
                       scen_name = "metrics_scen", y_lim = c(0,1))
p2 <- plot_xi_patterns(wd = wd, data_dir = data_dir, analysis_dir = analysis_dir,
                       scen_pop = 1112, iter_pop = 1, scen_samp = 111211,
                       scen_name = "metrics_scen", y_lim = c(-2,2))
p3 <- plot_Phi_patterns(wd = wd, data_dir = data_dir, analysis_dir = analysis_dir,
                        scen_pop = 1112, iter_pop = 1, scen_samp = 111211,
                        scen_name = "metrics_scen", y_lim = c(0,1))
ggarrange(p1, p2, p3, ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")

# Stratified cluster sampling
p4 <- plot_pi_patterns(wd = wd, data_dir = data_dir, analysis_dir = analysis_dir,
                       scen_pop = 1112, iter_pop = 1, scen_samp = 111212,
                       scen_name = "metrics_scen")
p5 <- plot_xi_patterns(wd = wd, data_dir = data_dir, analysis_dir = analysis_dir,
                       scen_pop = 1112, iter_pop = 1, scen_samp = 111212,
                       scen_name = "metrics_scen", y_lim = c(-2,2))
p6 <- plot_Phi_patterns(wd = wd, data_dir = data_dir, analysis_dir = analysis_dir,
                        scen_pop = 1112, iter_pop = 1, scen_samp = 111212,
                        scen_name = "metrics_scen")
ggarrange(p4, p5, p6, ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")

ggarrange(p1, p3, p4, p6, ncol = 2, nrow = 2, common.legend = TRUE, 
          legend = "right")

# SRS
p7 <- plot_pi_patterns(wd = wd, data_dir = data_dir, analysis_dir = analysis_dir,
                       scen_pop = 1112, iter_pop = 1, scen_samp = 111213,
                       scen_name = "metrics_scen", y_lim = c(0,1))
p8 <- plot_xi_patterns(wd = wd, data_dir = data_dir, analysis_dir = analysis_dir,
                       scen_pop = 1112, iter_pop = 1, scen_samp = 111213,
                       scen_name = "metrics_scen", y_lim = c(-2,2))
p9 <- plot_Phi_patterns(wd = wd, data_dir = data_dir, analysis_dir = analysis_dir,
                        scen_pop = 1112, iter_pop = 1, scen_samp = 111213,
                        scen_name = "metrics_scen", y_lim = c(0,1))
# Plot pi for the three scenarios
ggarrange(p1, p4, p7, nrow = 1, common.legend = TRUE)


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


#================== OLD CODE ===================================================
# plot_pi_cov <- data.frame(
#   coverage = c(metrics_SRS_s$xi_cover_avg, metrics_SRS_unsup$xi_cover_avg, 
#                metrics_SRS_ws$xi_cover_avg,
#                metrics_Strat_s$xi_cover_avg, metrics_Strat_unsup$xi_cover_avg,
#                metrics_Strat_ws$xi_cover_avg, 
#                metrics_Clus_s$xi_cover_avg, metrics_Clus_unsup$xi_cover_avg,
#                metrics_Clus_ws$xi_cover_avg),
#   pi = rep(c("1", "2", "3"), times=9),
#   model = rep(rep(c("SOLCA", "WOLCA", "WSOLCA"), each=3), times=3),
#   sampling = rep(c("SRS", "Stratified", "Stratified\n Cluster"), each=9)
# )
# pi_plot <- plot_pi_cov %>%  ggplot(aes(x=pi, y=coverage, fill=model)) +
#   theme_bw() + scale_fill_brewer(palette="RdYlBu") + 
#   labs(x = expression(pi), y = "Coverage", fill = "Model") +
#   geom_point(size=2.5, position = position_dodge(0.3), pch = 21) +
#   geom_hline(aes(yintercept=0.95)) +
#   scale_y_continuous(limits=c(0,1), breaks=sort(c(seq(0, 1, length.out=5), 0.95))) +
#   facet_grid(~sampling)
# 
# plot_theta_cov <- data.frame(
#   coverage = c(metrics_SRS_s$theta_cover_avg, metrics_SRS_ws$theta_cover_avg,
#                metrics_SRS_unsup$theta_cover_avg, metrics_Strat_s$theta_cover_avg,
#                metrics_Strat_ws$theta_cover_avg, metrics_Strat_unsup$theta_cover_avg),
#   theta = rep(c("theta_1", "theta_2", "theta_3"), times=6),
#   model = rep(rep(c("SOLCA (Unweighted)", "WSOLCA (Weighted)",
#                     "WOLCA (Unsupervised)"), each=3), times=2),
#   sampling = rep(c("SRS", "Stratified"), each=9)
# )
# theta_plot <- plot_theta_cov %>%  ggplot(aes(x=theta, y=coverage, color=model)) +
#   theme_bw() +
#   labs(x = expression(theta), y = "Coverage", color = "Model") +
#   geom_point(size=2.5, position = position_dodge(0.3)) +
#   geom_hline(aes(yintercept=0.95)) +
#   scale_y_continuous(limits=c(0,1), breaks=sort(c(seq(0, 1, length.out=5), 0.95))) +
#   facet_grid(~sampling) +
#   ggtitle(paste0(scenario, " scenario coverage for θ over 100 samples"))
# 
# plot_xi_cov <- data.frame(
#   coverage = c(metrics_SRS_s$xi_cover_avg, metrics_SRS_ws$xi_cover_avg,
#                metrics_SRS_unsup$xi_cover_avg, metrics_Strat_s$xi_cover_avg,
#                metrics_Strat_ws$xi_cover_avg, metrics_Strat_unsup$xi_cover_avg),
#   xi = rep(c("xi_1", "xi_2", "xi_3", "xi_4", "xi_5", "xi_6"), times=6),
#   model = rep(rep(c("SOLCA (Unweighted)", "WSOLCA (Weighted)",
#                     "WOLCA (Unsupervised)"), each=6), times=2),
#   sampling = rep(c("SRS", "Stratified"), each=18)
# )
# xi_plot <- plot_xi_cov %>%  ggplot(aes(x=xi, y=coverage, color=model)) +
#   theme_bw() +
#   theme(legend.position="top") +
#   labs(x = expression(xi), y = "Coverage", color = "Model") +
#   geom_point(size=2.5, position = position_dodge(0.3)) +
#   geom_hline(aes(yintercept=0.95)) +
#   scale_y_continuous(limits=c(0,1), breaks=sort(c(seq(0, 1, length.out=5), 0.95))) +
#   facet_grid(~sampling) +
#   ggtitle(paste0(scenario, " scenario coverage for ξ over 100 samples"))
# 
# ggarrange(pi_plot, theta_plot, nrow = 1, common.legend = TRUE, legend = "top")
# xi_plot

# load(paste0(wd, analysis_dir, "metrics_scen", 111211, ".RData"))
# 
# 
# create_table(metrics_s = metrics_all$metrics_s, 
#              metrics_ws = metrics_all$metrics_ws, 
#              metrics_unsup = metrics_all$metrics_unsup,
#              scen_samp = scen_samp)
# 
# ### Create table of metrics with bias and variance
# # Inputs:
# #   metrics_s: Summary metrics for SOLCA
# #   metrics_ws: Summary metrics for WSOLCA
# #   metrics_unsup: Summary metrics for WOLCA
# #   scen_samp: Sampling scenario
# # Output: Formatted table with absolute bias, CI width, and coverage
# create_table <- function(metrics_s, metrics_ws, metrics_unsup, scen_samp) {
#   metrics_summ <- as.data.frame(matrix(NA, nrow=3, ncol=12))
#   colnames(metrics_summ) <- c("Sampling Scheme", "Model", 
#                               "K Bias^2", "$\\pi$ Bias^2", "$\\pi$ CI width", 
#                               "$\\theta$ Bias^2", "$\\theta$ CI width", 
#                               "$\\xi$ Bias^2", "$\\xi$ CI width", 
#                               "$\\pi$ Coverage","$\\theta$ Coverage", "$\\xi$ Coverage")
#   metrics_summ[, 1] <- rep("Stratified", 3)
#   metrics_summ[, 2] <- rep(c("Unwtd(sOFMM)", "Wtd(wsOFMM)", "Unsup(wOFMM)"), 1)  ## latent versions
#   output_inds <- 1:7
#   metrics_summ[1, -c(1,2)] <- c(metrics_s[output_inds], 
#                                 mean(metrics_s$xi_cover_avg), 
#                                 mean(metrics_s$theta_cover_avg),
#                                 mean(metrics_s$xi_cover_avg))
#   metrics_summ[2, -c(1,2)] <- c(metrics_ws[output_inds], 
#                                 mean(metrics_ws$xi_cover_avg), 
#                                 mean(metrics_ws$theta_cover_avg),
#                                 mean(metrics_ws$xi_cover_avg))
#   metrics_summ[3, -c(1,2)] <- c(metrics_unsup[output_inds], 
#                                 mean(metrics_unsup$xi_cover_avg), 
#                                 mean(metrics_unsup$theta_cover_avg),
#                                 mean(metrics_unsup$xi_cover_avg))
#   metrics_summ %>% 
#     gt(caption = paste0(scen_samp, " scenario metrics of posterior parameter estimates, averaged over 100 samples")) %>%
#     cols_label("$\\pi$ Bias^2" = "π Bias^2", "$\\theta$ Bias^2" = "θ Bias^2", 
#                "$\\xi$ Bias^2" = "ξ Bias^2",
#                "$\\pi$ CI width" = "π CI width", "$\\theta$ CI width" = "θ CI width", 
#                "$\\xi$ CI width" = "ξ CI width", "$\\pi$ Coverage" = "π Coverage",
#                "$\\theta$ Coverage" = "θ Coverage", "$\\xi$ Coverage" = "ξ Coverage") %>%
#     fmt_number(
#       columns = 3:12,
#       decimals = 4)
# }
# 
# 
# 
# #================ TROUBLESHOOTING NUMBER OF CLUSTERS, K ========================
# 
# wrong_K <- which(K_dist == 1)
# model <- "wsOFMM"
# for (i in 1:length(wrong_K)) {
#   samp_n <- wrong_K[i]
#   # wsOFMM model includes a variance adjustment
#   sim_res_path <- paste0(wd, res_dir, model, "_results_adjRcpp_scen", scen_samp, 
#                          "_samp", samp_n, ".RData")
#   load(sim_res_path)
#   analysis <- res$analysis_adj
#   names(analysis) <- str_replace_all(names(analysis), "_adj", "")
#   print(paste0("i: ", i, ", K: ", length(analysis$pi_med)))
#   
#   MCMC_path <- paste0(wd, res_dir, model, "_results_scen", scen_samp, 
#                       "_samp", samp_n, ".RData")  # Output file
#   load(MCMC_path)
#   K_med <- res_MCMC$post_MCMC_out$K_med
#   print(paste0("K_med: ", K_med))
# }
# 
# samp_n <- 24
# model <- "sOFMM"
# MCMC_path <- paste0(wd, res_dir, model, "_results_scen", scen_samp, 
#                     "_samp", samp_n, ".RData")  # Output file
# load(MCMC_path)
# # sOFMM
# MCMC_out <- res$MCMC_out
# # Get median number of classes with >= 5% of individuals, over all iterations
# M <- dim(MCMC_out$pi_MCMC)[1]  # Number of stored MCMC iterations
# K_med <- round(median(rowSums(MCMC_out$pi_MCMC >= 0.05)))
# print(paste0("sOFMM K_fixed: ", dim(MCMC_out$pi_MCMC)[2]))
# print(paste0("sOFMM K_med: ", K_med))
# 
# # Cluster individuals into reduced number of classes using agglomerative clustering
# # Calculate pairwise distance matrix using Hamming distance: proportion of 
# # iterations where two individuals have differing class assignments
# distMat_s <- hamming.distance(t(MCMC_out$c_all_MCMC))
# dendrogram <- hclust(as.dist(distMat_s), method = "complete") # Hierarchical clustering dendrogram
# dend_plot <- as.dendrogram(dendrogram)
# dend_plot %>% set("branches_k_color", k = 3) %>% plot()
# red_c_all <- cutree(dendrogram, k = K_med)                  # Group individuals into K_med classes
# table(red_c_all)
# 
# # Reduce and reorder parameter estimates using new classes
# p <- 30
# q <- 2
# d <- 4
# pi <- matrix(NA, nrow = M, ncol = K_med)
# theta <- array(NA, dim = c(M, p, K_med, d))
# xi <- array(NA, dim = c(M, K_med, q))
# for (m in 1:M) {
#   iter_order <- relabel_red_classes[m, ]
#   pi_temp <- MCMC_out$pi_MCMC[m, iter_order]
#   pi[m, ] <- pi_temp / sum(pi_temp)
#   theta[m, , , ] <- MCMC_out$theta_MCMC[m, , iter_order, ]
#   xi[m, , ] <- MCMC_out$xi_MCMC[m, iter_order, ]
# }
# post_MCMC_out <- list(K_med = K_med, pi = pi, theta = theta, xi = xi)
# plot(pi[,1], type = "l", ylab = "Pi")
# plot(pi[,2], type = "l", ylab = "Pi")
# plot(pi[,3], type = "l", ylab = "Pi")
# plot(pi[,4], type = "l", ylab = "Pi")
# 
# get_mode <- function(v) {
#   uniqv <- unique(v)
#   uniqv[which.max(tabulate(match(v, uniqv)))]
# }
# # For each iteration, relabel new classes using the most common old class assignment
# relabel_red_classes <- matrix(NA, nrow = M, ncol = K_med)   # Apply new classes to each iteration
# for (k in 1:K_med) {
#   relabel_red_classes[, k] <- apply(as.matrix(MCMC_out$c_all_MCMC[, red_c_all == k]), 
#                                     1, get_mode)
# }
# # Posterior median estimate for theta across iterations
# theta_med_temp <- apply(theta, c(2, 3, 4), median)
# # Posterior modal exposure categories for each exposure item and reduced class
# theta_modes <- apply(theta_med_temp, c(1, 2), which.max)
# # Identify unique classes
# unique_classes <- which(!duplicated(theta_modes, MARGIN = 2))
# print(theta_modes)
# 
# 
# # wsOFMM
# model <- "wsOFMM"
# # wsOFMM model includes a variance adjustment
# sim_res_path <- paste0(wd, res_dir, model, "_results_adjRcpp_scen", scen_samp, 
#                        "_samp", samp_n, ".RData")
# load(sim_res_path)
# MCMC_out <- res_MCMC$MCMC_out
# post_MCMC_out <- res_MCMC$post_MCMC_out
# print(paste0("wsOFMM K_fixed: ", dim(MCMC_out$pi_MCMC)[2]))
# print(paste0("wsOFMM K_med: ", post_MCMC_out$K_med))
# 
# distMat_ws <- hamming.distance(t(MCMC_out$c_all_MCMC))
# dendrogram <- hclust(as.dist(distMat_ws), method = "complete") # Hierarchical clustering dendrogram
# dend_plot <- as.dendrogram(dendrogram)
# dend_plot %>% set("branches_k_color", k = 3) %>% plot()
# red_c_all <- cutree(dendrogram, k = K_med)                  # Group individuals into K_med classes
# table(red_c_all)
# 
# plot(res_MCMC$post_MCMC_out$pi[,1], type = "l", ylab = "Pi")
# plot(res_MCMC$post_MCMC_out$pi[,2], type = "l", ylab = "Pi")
# plot(res_MCMC$post_MCMC_out$pi[,3], type = "l", ylab = "Pi")
# 
# # Posterior median estimate for theta across iterations
# theta_med_temp <- apply(post_MCMC_out$theta, c(2, 3, 4), median)
# # Posterior modal exposure categories for each exposure item and reduced class
# theta_modes <- apply(theta_med_temp, c(1, 2), which.max)
# # Identify unique classes
# unique_classes <- which(!duplicated(theta_modes, MARGIN = 2))
# print(theta_modes)
# 
# theta_iter <- MCMC_out$theta_MCMC[1000,,,]
# theta_iter_modes <- apply(theta_iter, c(1, 2), which.max)
# theta_iter_modes
# 
# leave <- c(24, 59, 80, 84, 94)
# samp_n_seq <- setdiff(1:100, leave)
# 
# for (k in 1:3) {
#   for (s in 1:2) {
#     print(paste0("k: ", k, ", s: ", s))
#     print(sum(sim_pop$Y_data==1 & sim_pop$true_Si==s & 
#                 sim_pop$true_Ci==k))
#     print(sum(sim_pop$true_Si==s & sim_pop$true_Ci==k))
#     print("sample")
#     print(sum(sim_samp$Y_data==1 & sim_samp$true_Si==s & 
#                 sim_samp$true_Ci==k))
#     print(sum(sim_samp$true_Si==s & sim_samp$true_Ci==k))
#   }
# }
# plot(MCMC_out$pi_MCMC[,1], type = "l", ylim=c(0.1, 0.7))
# lines(MCMC_out$pi_MCMC[,2], col = "red")
# lines(MCMC_out$pi_MCMC[,3], col = "blue")
# 
# 
