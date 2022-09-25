# Model Code
This folder includes files for running the sOFMM and wsOFMM models on simulated datasets.

### The following files are used to run the sOFMM and wsOFMM models
* `sOFMM_main.m`: Runs the supervised overfitted finite mixture model on a given dataset.
    * Input: Simulated dataset from "wsOFMM/Data" folder.
    * Outputs: MCMC and posterior results output in "wsOFMM/Results" folder.
* `wsOFMM_main.m`: Runs the weighted supervised overfitted finite mixture model on a given dataset.
    * Input: Simulated dataset from "wsOFMM/Data" folder.
    * Outputs: MCMC and posterior results output in "wsOFMM/Results" folder.

### Additional helper files 
* `analyze_results.m`: Helper function used to obtain true number of latent classes, posterior parameter estimates, and analysis metrics. Used by "sOFMM_main.m" and "wsOFMM_main.m."
* `drchrnd.m`: Helper function that generates a random sample from a Dirichlet distribution. Used by "sOFMM_main.m" and "wsOFMM_main.m."
* `get_data_vars.m`: Helper function that takes in sample data and outputs relevant variables for the sOFMM model. Used by "sOFMM_main.m."
* `init_OFMM_params.m`: Helper function that initializes priors and variables for the OFMM model. Used by "sOFMM_main.m" and "wsOFMM_main.m."
* `init_probit_params.m`: Helper function that initializes priors and variables for the probit model. Used by "sOFMM_main.m" and "wsOFMM_main.m."
* `post_process.m`: Helper function that conducts post-processing to remove extraneous empty classes and relabel class assignments. Used by "sOFMM_main.m" and "wsOFMM_main.m."
* `run_MCMC.m`: Helper function that runs the Gibbs Sampler MCMC algorithm to obtain posteriors for the sOFMM model. Used by "sOFMM_main.m."
* `truncnormrnd.m`: Helper function that generates a random sample from a truncated normal distribution. Used by "sOFMM_main.m" and "wsOFMM_main.m."
* `update_MCMC.m`: Helper function that updates the posterior distributions of the parameters and variables for the sOFMM model. Used by "sOFMM_main.m."
* `wtd_get_dara_vars.m`: Helper function that takes in sample data and outputs relevant variables for the wsOFMM model. Used by "wsOFMM_main.m."
* `wtd_run_MCMC.m`: Helper function that runs the Gibbs Sampler MCMC algorithm to obtain posteriors for the wsOFMM model. Used by "wsOFMM_main.m."
* `wtd_update_MCMC.m`: Helper function that updates the posterior distributions of the parameters and variables for the wsOFMM model. Used by "wsOFMM_main.m."
* `matlab_batchsubmit.sh`: Bash script for submitting a cluster job array to run the models on sample data over multiple samples.
* `matlab_batchsubmit_pop.sh`: Bash script for submitting a cluster job array to run the models on population data over multiple iterations.
