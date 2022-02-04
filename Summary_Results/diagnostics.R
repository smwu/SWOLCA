library(R.matlab)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(dplyr)
library(rstan)
library(bayesplot)
library(coda)

setwd("C:/Users/Lang/Documents/Harvard/Research/Briana/supRPC/wsOFMM")
setwd("/n/home01/stephwu18/wsOFMM/")
scenarios <- seq(5, 16)

len <- length(scenarios)
# K_red <- numeric(len)
# Phi_mse <- numeric(len)
# theta_mse <- numeric(len)
# sens <- numeric(len)
# spec <- numeric(len)
# dic6 <- numeric(len)
# aebic <- numeric(len)
# corr_s_c <- numeric(len)
# corr_s_y <- numeric(len)
output <- as.data.frame(matrix(NA, nrow=len*2, ncol=12))
colnames(output) <- c("Scenario", "Model", "K_red", "Phi_MSE", "Theta_MSE", "Sens", 
                      "Spec", "DIC6", "AEBIC", "Corr(S,C)", "Corr(S,Y)", "Corr(C,Y)")

for (scenario in scenarios) {
  if (scenario > 4) {
    results_s <- readMat(paste0("summary_sOFMM_scen", scenario, "_sample", ".mat"))
    results_ws <- readMat(paste0("summary_wsOFMM_scen", scenario, "_sample", ".mat"))
  } else {
    results_s <- readMat(paste0("summary_sOFMM_scen", scenario, "_pop", ".mat"))
    results_ws <- readMat(paste0("summary_wsOFMM_scen", scenario, "_pop", ".mat"))
  }
  res_s <- results_s$res
  res_ws <- results_ws$res
  start <- scenario-scenarios[1]+1
  output[start*2-1, -2] <- c(scenario, res_s[[1]][1], res_s[[7]][1], res_s[[8]][1], 
                            res_s[[3]][1], res_s[[4]][1], res_s[[5]][1], res_s[[6]][1], 
                            res_s[[9]], res_s[[10]], res_s[[11]])
  output[start*2-1, 2] <- "sOFMM"
  output[start*2, -2] <- c(scenario, res_ws[[1]][1], res_ws[[7]][1], res_ws[[8]][1], 
                            res_ws[[3]][1], res_ws[[4]][1], res_ws[[5]][1], res_ws[[6]][1], 
                            res_ws[[9]], res_ws[[10]], res_ws[[11]])
  output[start*2, 2] <- "wsOFMM"
}

p1 <- output %>% 
  ggplot(aes(x=Scenario, y=Phi_MSE, fill=Model)) +
  geom_bar(stat = "identity", position="dodge") + 
  ggtitle("Phi_MSE") + scale_x_continuous(breaks=5:16)

p2 <- output %>% 
  ggplot(aes(x=Scenario, y=Theta_MSE, fill=Model)) +
  geom_bar(stat = "identity", position="dodge") + 
  ggtitle("Theta_MSE") + scale_x_continuous(breaks=5:16)

p3 <- output %>% 
  ggplot(aes(x=Scenario, y=Sens, fill=Model)) +
  geom_bar(stat = "identity", position="dodge") + 
  ggtitle("Sens") + scale_x_continuous(breaks=5:16)

p4 <- output %>% 
  ggplot(aes(x=Scenario, y=Spec, fill=Model)) +
  geom_bar(stat = "identity", position="dodge") + 
  ggtitle("Spec") + scale_x_continuous(breaks=5:16)

p5 <- output %>% 
  ggplot(aes(x=Scenario, y=DIC6, fill=Model)) +
  geom_bar(stat = "identity", position="dodge") + 
  ggtitle("DIC-6") + scale_x_continuous(breaks=5:16)

p6 <- output %>% 
  ggplot(aes(x=Scenario, y=AEBIC, fill=Model)) +
  geom_bar(stat = "identity", position="dodge") + 
  ggtitle("AEBIC") + scale_x_continuous(breaks=5:16)

ggarrange(p1, p2, p3, p4, p5, p6, ncol=2, nrow=3, common.legend=TRUE)
# facet_grid(.~Scenario)

# mcmc_scen3_iter68 <- readMat("wsOFMM_MCMC_scen3_iter68.mat")
# data <- mcmc_scen3_iter68$MCMC.out
# names(data) <- c("pi", "theta", "c_i", "xi", "z_i", "loglik", "runtime")
# data_mcmc <- mcmcUpgrade(mcmc(data))
# coda::traceplot(data_mcmc)
