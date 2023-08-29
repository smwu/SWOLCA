# Simulated Data

This folder includes the NHANES application data, as well as some examples of simulated datasets.

### NHANES application data
`nhanes1518_adult_low_f_12jul2023.csv` contains the final cleaned and processed dataset, which includes includes data for the categorical dietary intake exposure, the binary hypertensive outcome, and the additional sociodemographic confounders.

### Simulated data
Due to storage constraints, only a few examples of the simulated datasets used for the simulation study are displayed here. For each scenario of the simulation study, models were run on 100 samples generated independently of each other, using different seeds for random number generation. 
* `simdata_scen111211_iter1_samp100.RData`: Default scenario with stratified sampling
* `simdata_scen111212_iter1_samp100.RData`: Scenario with stratified cluster sampling
* `simdata_scen111213_iter1_samp100.RData`: Scenario with simple random sampling
