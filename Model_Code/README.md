# Model Code
This folder includes files for running the SWOLCA, SOLCA, and WOLCA models on simulated datasets and on the NHANES dataset.

### Main model code for simulations
* `WSOLCA_main_Rcpp.R`: Runs SWOLCA on a simulated dataset and outputs MCMC and posterior results
* `SOLCA_main_Rcpp.R`: Runs unweighted SOLCA on a simulated dataset and outputs MCMC and posterior results
* `WOLCA_main_Rcpp.R`: Runs two-step WOLCA on a simulated dataset and outputs MCMC and posterior results

### Main model code for data application
* `WSOLCA_application_covs.R`: Runs SWOLCA on NHANES data
* `SOLCA_application_covs.R`: Runs unweighted SOLCA on NHANES data
* `WOLCA_application_covs.R`: Runs two-step WOLCA on NHANES data

### Helper files 
* `helper_functions.R`: Helper functions for running the models
* `main_Rcpp_functions.cpp`: Rcpp helper functions for running the Gibbs sampler
* `WSOLCA_main.stan`: Specification of the SWOLCA model in Stan to compute the posterior Hessian for the post-processing variance correction
* `R_batch_submit.sh`: Bash script for submitting R script as a cluster job array to run the models over multiple simulation iterations

### `Old_Code` folder
This folder contains old code for implementing SWOLCA in Stan and in R (without C++). 
* `WSOLCA.stan`: Specification of the model in Stan
* `stan_wsOFMM.R`: Runs SWOLCA in Stan on a simulated dataset and outputs MCMC and posterior results
* `stan_batch_submit.sh`: Bash script for submitting the Stan implementation as a cluster job array to run over multiple simulation iterations
* `WSOLCA_main.R`: Runs SWOLCA completely in R without any interface with C++
* `troubleshooting.R`: Testing congruence of R and Rcpp implementations
