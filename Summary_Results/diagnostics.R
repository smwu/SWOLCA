library(R.matlab)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(dplyr)
library(rstan)
library(bayesplot)
library(coda)

setwd("C:/Users/Lang/Documents/Harvard/Research/Briana/supRPC/wsOFMM/Summary_Results")
setwd("/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/Summary_Results")
# scenarios <- seq(5, 16)
# 
# len <- length(scenarios)
# 
# output <- as.data.frame(matrix(NA, nrow=len*2, ncol=12))
# colnames(output) <- c("Scenario", "Model", "K_red", "Phi_MSE", "Theta_MSE", "Sens", 
#                       "Spec", "DIC6", "AEBIC", "Corr(S,C)", "Corr(S,Y)", "Corr(C,Y)")
# 
# for (scenario in scenarios) {
#   if (scenario > 4) {
#     results_s <- readMat(paste0("summary_sOFMM_scen", scenario, "_sample", ".mat"))
#     results_ws <- readMat(paste0("summary_wsOFMM_scen", scenario, "_sample", ".mat"))
#   } else {
#     results_s <- readMat(paste0("summary_sOFMM_scen", scenario, "_pop", ".mat"))
#     results_ws <- readMat(paste0("summary_wsOFMM_scen", scenario, "_pop", ".mat"))
#   }
#   res_s <- results_s$res
#   res_ws <- results_ws$res
#   start <- scenario-scenarios[1]+1
#   output[start*2-1, -2] <- c(scenario, res_s[[1]][1], res_s[[7]][1], res_s[[8]][1], 
#                             res_s[[3]][1], res_s[[4]][1], res_s[[5]][1], res_s[[6]][1], 
#                             res_s[[9]], res_s[[10]], res_s[[11]])
#   output[start*2-1, 2] <- "sOFMM"
#   output[start*2, -2] <- c(scenario, res_ws[[1]][1], res_ws[[7]][1], res_ws[[8]][1], 
#                             res_ws[[3]][1], res_ws[[4]][1], res_ws[[5]][1], res_ws[[6]][1], 
#                             res_ws[[9]], res_ws[[10]], res_ws[[11]])
#   output[start*2, 2] <- "wsOFMM"
# }
# 
# p1 <- output %>% 
#   ggplot(aes(x=Scenario, y=Phi_MSE, fill=Model)) +
#   geom_bar(stat = "identity", position="dodge") + 
#   ggtitle("Phi_MSE") + scale_x_continuous(breaks=5:16)
# 
# p2 <- output %>% 
#   ggplot(aes(x=Scenario, y=Theta_MSE, fill=Model)) +
#   geom_bar(stat = "identity", position="dodge") + 
#   ggtitle("Theta_MSE") + scale_x_continuous(breaks=5:16)
# 
# p3 <- output %>% 
#   ggplot(aes(x=Scenario, y=Sens, fill=Model)) +
#   geom_bar(stat = "identity", position="dodge") + 
#   ggtitle("Sens") + scale_x_continuous(breaks=5:16)
# 
# p4 <- output %>% 
#   ggplot(aes(x=Scenario, y=Spec, fill=Model)) +
#   geom_bar(stat = "identity", position="dodge") + 
#   ggtitle("Spec") + scale_x_continuous(breaks=5:16)
# 
# p5 <- output %>% 
#   ggplot(aes(x=Scenario, y=DIC6, fill=Model)) +
#   geom_bar(stat = "identity", position="dodge") + 
#   ggtitle("DIC-6") + scale_x_continuous(breaks=5:16)
# 
# p6 <- output %>% 
#   ggplot(aes(x=Scenario, y=AEBIC, fill=Model)) +
#   geom_bar(stat = "identity", position="dodge") + 
#   ggtitle("AEBIC") + scale_x_continuous(breaks=5:16)
# 
# ggarrange(p1, p2, p3, p4, p5, p6, ncol=2, nrow=3, common.legend=TRUE)
# # facet_grid(.~Scenario)




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
