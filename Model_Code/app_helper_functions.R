#===================================================
## Helper functions for data application
## Programmer: SM Wu   
## Data: NHANES Application with Covariates  
## Date updated: 2023/07/14
#===================================================

# Load libraries
library(plyr)
library(dplyr)
library(tidyr)
library(forcats)
library(LaplacesDemon)
library(truncnorm)
library(stringr)
library(R.matlab)
library(fastDummies)
library(matrixStats)
library(Matrix)
library(gtools)
library(e1071)
library(rstan)
library(survey)
library(Rcpp)
library(RcppArmadillo)
library(RcppTN)
library(ggpubr)

# `process_data` reads in data and converts it into a processed format for the 
# 'WSOLCA_app_covs_Rcpp' function
# Input:
#   data_path: String specifying path for input data
#   covs: String vector of covariates to include in regression. Default is NULL
#   formula: String specifying formula for regression model
# Output: list 'data_vars' containing the following objects:
#   x_mat: Matrix of categorical exposure variables for all individuals; nxp
#   y_all: Vector of binary outcome for all individuals; nx1 
#   s_all: Vector of stratum indicator for all individuals; nx1
#   clus_id_all: Vector of cluster indicator for all individuals; nx1
#   sample_wt: Vector of sampling weights for all individuals; nx1
#   V: Model matrix for regression; nxq
#   V_data: Matrix of cleaned variables to use in regression; nx(num_vars)
# Example: process_data(data_path = "C:/Documents", 
#                       covs = c("age_cat", "racethnic", "smoker", "physactive"),
#                       formula = "~ age_cat + racethnic + smoker + physactive")
process_data <- function(data_path, covs = NULL, formula = NULL) {
  # Read in data
  data_vars_f_low <- read.csv(data_path)
  # Drop those with refused or unknown education data
  # Complete n = 2003
  data_vars_f_low <- data_vars_f_low %>% filter(!DMDEDUC2 %in% c(7, 9))
      # # Drop those with race/ethnic equal to Other/Mixed 
      # # Complete n = 1921
      # data_vars_f_low <- data_vars_f_low %>% filter(!RIDRETH3 == 7)
  # Drop legumes (vegs) because duplicate of legumes (proteins), just computed
  # as cup eq vs oz eq
  data_vars_f_low <- data_vars_f_low %>% select(-leg_veg)

  # Obtain exposure and outcome data
  x_mat <- as.matrix(data_vars_f_low %>% select(citrus:drinks))
  y_all <- data_vars_f_low$BP_flag
  
  # Get stratum IDs
  s_all <- data_vars_f_low$SDMVSTRA  # 30 strata
  # Create unique nested PSU IDs using stratum IDs
  clus_id_all <- s_all * 10 + data_vars_f_low$SDMVPSU
  # Get sampling weights
  sample_wt <- data_vars_f_low$dietwt4yr
  
  # Other covariates to be included in the probit model 
  if (is.null(covs)) {
    n <- length(y_all)
    V <- matrix(1, nrow = n) 
    q <- 1
  } else {
    V_data <- data_vars_f_low %>% 
      select(RIDAGEYR, RIDRETH3, DMDEDUC2, smoker, Phys_Active) %>%
      mutate(
        age_cat = factor(case_when(  
          20 <= RIDAGEYR & RIDAGEYR <= 39 ~ 1,
          40 <= RIDAGEYR & RIDAGEYR <= 59 ~ 2,
          RIDAGEYR >= 60 ~ 3)),
        # RIDRETH3: Race/Hispanic origin w/ NH Asian: 1=Mex_Amer, 2=Other_Hisp,
        # 3=NH_White, 4=NH_Black, 6=NH_Asian, 7=Other/Mixed
        racethnic = factor(case_when(
          RIDRETH3 == 3 ~ 1,  # NH White
          RIDRETH3 == 4 ~ 2,  # NH Black
          RIDRETH3 == 6 ~ 3,  # NH Asian
          # RIDRETH3 == 1 ~ 4,  # Mexican-American
          # RIDRETH3 == 2 ~ 5,  # Other Hispanic
          RIDRETH3 %in% c(1, 2) ~ 4,  # Mexican-American/Other Hispanic
          RIDRETH3 == 7 ~ 5,  # Other/Mixed
          .default = NA)),
        # DMDEDUC2: Education level for adults 20+: 1= <9th, 2=9-11th, 3=HS/GED,
        # 4=Some college/AA, 5=college grad or above, 7=refused, 9=don't know
        educ = factor(case_when(
          DMDEDUC2 %in% c(4, 5) ~ 1,  # At least some college
          DMDEDUC2 == 3 ~ 2,  # HS/GED
          DMDEDUC2 %in% c(1, 2) ~ 3,  # Less than HS
          .default = NA)),
        smoker = factor(smoker),
        physactive = factor(Phys_Active),
        .keep = "unused")
    V_data <- V_data %>% select(covs)
    
    # Regression design matrix without class assignment, nxq
    # Exclude stratifying variable as well
    if (is.null(formula)) {
      stop("Error: Need to specify a string formula for the probit regression.")
    }
    V <- model.matrix(as.formula(formula), data = V_data)
  }
  
  # Return processed data
  data_vars <- list(x_mat = x_mat, y_all = y_all, s_all = s_all, 
                    clus_id_all = clus_id_all, sample_wt = sample_wt, 
                    V = V)
  if (!is.null(covs)) {
    data_vars$V_data <- V_data
  }
  return(data_vars)
}
