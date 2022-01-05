# Weighted Supervised Overfitted Finite Mixture Model (wsOFMM)

The weighted supervsied overfitted finite mixture model (wsOFMM) is a Bayesian clustering technique that assigns individuals to latent clusters based on categorical exposure data and binary outcome data. It uses a weighted pseudolikelihood method to incorporate survey weights and adjust for selection bias. An overfitted finite mixture model, similar to a latent class model, is jointly modeled with a probit regression model in order to generate cluster profiles that are informed by the probability of a positive binary outcome.

The parameters of the wsOFMM model are estimated using a two-step sampling algorithm:

1. In the adaptive sampling algorithm, the wsOFMM model is run to determine the appropriate number of latent clusters. This step may be skipped if the number of latent clusters is known a priori.
2. In the fixed sampling algorithm, model parameters are estimated by re-runing the wsOFMM model with a fixed number of latent clusters determined by step 1.

# Code and Materials

This repo contains the following MATLAB files for creating simulated data and running the wsOFMM model:

- `wsOFMM_main.m`: This includes main and local functions needed to run the wsOFMM model. To run the example data, one of the simulated data files, titled "simdata...." is needed.
- `sim_wsRPC_scen1.m': Code to generate simulation data and weights where sampling follows Scenario 1: equal subpopulations; sample 100\% from each subpopulation (full population sampled).
- `sim_wsRPC_scen2.m': Code to generate simulation data and weights where sampling follows Scenario 2: unequal subpopulations; sample 100\% from each subpopulation (full population sampled).
- `sim_wsRPC_scen3.m': Code to generate simulation data and weights where sampling follows Scenario 3: equal subpopulations; sample 5\% from each subpopulation (proportional allocation).
- `sim_wsRPC_scen4.m': Code to generate simulation data and weights where sampling follows Scenario 4: unequal subpopulations; sample 5\% from each subpopulation (proportional allocation).
- `sim_wsRPC_scen5.m': Code to generate simulation data and weights where sampling follows Scenario 5: unequal subpopulations; sample 1000 from each subpopulation (equal allocation).
- `simdata_wsRPC_scen1_iter1.mat`: One simulated dataset for Scenario 1.
- `simdata_wsRPC_scen2_iter1.mat`: One simulated dataset for Scenario 2.
- `simdata_wsRPC_scen3_iter1.mat`: One simulated dataset for Scenario 3.
- `simdata_wsRPC_scen4_iter1.mat`: One simulated dataset for Scenario 4.
- `simdata_wsRPC_scen5_iter1.mat`: One simulated dataset for Scenario 5.

The "simdata..." datasets contain the following variables:
- `subpop_samp`: Subpopulation assignment for all sampled individuals
- `true_Ci`: True global class memberships for all sampled individuals
- 'true_Li`: True subpopualtion-specific local class memberships for all sampled individuals
- `sample_wt`: Survey weights for all sampled individuals
- `norm_const`: Normalization constants for the weights in each subpopulation
- `true_global_patterns`: True consumption patterns for global classes
- `true_local_patterns`: True consumption patterns for local classes
- `global_thetas`: True global item-response probabilities
- `local_thetas`: True local item-response probabilities
- `nu`: Subpopulation-specific probability of global assignment
- `true_G`: True global/local assignments for each item and subpopulation
- `sample_data`: Observed consumption data for all sampled individuals
- `true_xi`: True probit model coefficients
- `true_Phi`: True probit model mean
- `true_y`: True binary outcome values for all sampled individuals
