# Supervised Weighted Overfitted Latent Class Analysis (SWOLCA)

This repo contains code used to run the simulation and data application analyses demonstrating the utility of **supervised weighted overfitted latent class analysis (SWOLCA)**, proposed in the paper, "_Derivation of outcome-dependent dietary patterns for low-income women obtained from survey data using a Supervised Weighted Overfitted Latent Class Analysis,_" by Stephanie M. Wu and Briana Stephenson. SWOLCA is applied to characterize dietary patterns associated with hypertensive outcomes among low-income women in the United States, using dietary intake and hypertension data from the 2015-2018 National Health and Nutrition Examination Surveys (NHANES). 

SWOLCA is a Bayesian model-based clustering technique for survey data. It assigns individuals to latent exposure-outcome patterns based on a multivariate categorical exposure and a binary outcome.  Complex survey design such as stratification, clustering, and informative sampling are adjusted for by using a pseudo-likelihood approach that integrates sampling weights into latent pattern elicitation for consistent estimation, alongside a post-processing projection for valid inference. A mixture reference coding scheme is used to incorporate interaction effects when characterizing the association between exposure and outcome while accommodating label switching. Parameters are estimated using a two-stage sampling algorithm. The first stage involves an _adaptive_ sampler to estimate the appropriate number of latent patterns. This step may be skipped if the number of latent clusters is known a priori. The second stage involves a _fixed_ sampler, where the SWOLCA model is re-run with the number of patterns derived from the adaptive sampler.

Two alternative methods are used for comparison in the simulation study: 1) an unweighted supervised overfitted latent class analysis (SOLCA) that ignores survey design, and 2) a two-step unsupervised weighted overfitted latent class analysis (WOLCA) that performs pattern derivation and outcome regression in two separate steps.

### This repo contains the following folders:
 * `Simulation_Code`: Code used to generate simulated data.
 * `Data_Code`: Code and raw data used to create the NHANES data application dataset. 
 * `Data`: Examples of simulated datasets, as well as the NHANES data application dataset. 
 * `Model_Code`: Code for the SWOLCA, SOLCA, and WOLCA models.
 * `Results`: Example MCMC and posterior results from running the models on the simulated data and the NHANES dataset.
 * `Analysis_Code`: Code used to summarize the MCMC output from the simulation study and create tables and figures.


