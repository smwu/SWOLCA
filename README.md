# Weighted Supervised Overfitted Finite Mixture Model (wsOFMM)

This repo contains code and partial results from a simulation study of the wsOFMM. View the description of the folders below to help navigate to relevant code and materials.

The **supervised overfitted finite mixture model (sOFMM)** is a Bayesian clustering technique that assigns individuals to latent classes based on categorical exposure data and binary outcome data. It includes an overfitted finite mixture model for reducing the dimension of a categorical exposure, combined with a probit regression model for measuring association with a binary outcome. The model sorts subjects into underlying classes with similar exposure and outcome patterns. The **weighted supervised overfitted finite mixture model (wsOFMM)** extends the sOFMM by using a weighted pseudolikelihood method to incorporate survey weights. It can be used to adjust for selection bias when using survey data. 

The parameters of both models are estimated using a two-step sampling algorithm: 1) In the adaptive sampling algorithm, the model is run to determine the number of underlying latent classes. This step may be skipped if the number of latent clusters is known a priori. 2) In the fixed sampling algorithm, model parameters are estimated by re-running the model with the fixed number of latent clusters determined by the adaptive sampling algorithm.

### This repo contains the following folders:
 * `Simulation_Code`: Includes code used to generate simulated data.
 * `Data`: Includes simulated datasets generated using functions in the "Simulation_Code" folder.
 * `Model_Code`: Includes code used to run the sOFMM and wsOFMM models on simulated data in the "Data" folder.
 * `Results`: Includes MCMC and posterior results output from running functions in the "Model_Code" folder.
 * `Summary_Results`: Includes code and results that summarize output over multiple iterations and samples. Used for analysis. 


