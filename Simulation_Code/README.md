# Simulation Code
This folder includes files for creating simulated datasets.

### The following files are used to generate simulated data for an entire population of size N=8000:
* `sim_no_local.m`: Generates population data for Scenario 1 (see DAG 1 below), where subpopulation membership is independent of class membership and no local latent classes exist.
* `sim_no_local_pi_subpop.m`: Generates population data for Scenario 2 (see DAG 2 below), where subpopulation membership affects class membership and no local latent classes exist.
* `sim_uneq.m`: Generates population data for Scenario 3 (see DAG 3 below), where subpopulation membership is independent of class membership but affects local latent classes.
* `sim_uneq_pi_subpop.m`: Generates population data for Scenario 4 (see DAG 4 below), where subpopulation membership affects class membership and affects local latent classes.

### The following files are used to generate simulated data for a sample of size n=400 drawn from an existing population:
* `sample_SRS.m`: Draws a simple random sample (SRS) from a pre-existing population dataset. Scenarios 5-8 correspond to draws from population datasets for Scenarios 1-4, respectively, generated using the above functions.
* `sample_strat_prop.m`: Draws a stratified random sample with proportional allocation by subpopulation, with 5% of each subpopulation sampled from a pre-existing population dataset. Scenarios 9-12 correspond to draws from population datasets for Scenarios 1-4, respectively, generated using the above functions.
* `sample_strat_eq.m`: Draws a stratified random sample with equal allocation by subpopulation, with 100 subjects per subpopulation sampled from a pre-existing population dataset. Scenarios 9-12 correspond to draws from population datasets for Scenarios 1-4, respectively, generated using the above functions.

### Additional helper files:
* `sample_indivs.m`: Helper function used in all of the above files. Samples from a population.
* `create_consump_data.m`: Helper function that generates observed consumption data <img src="https://latex.codecogs.com/gif.latex?X" />. Called by "sim_uneq.m" and "sim_uneq_pi_subpop.m" files.
* `create_consump_data_no_local.m`: Helper function that generates observed consumption data <img src="https://latex.codecogs.com/gif.latex?X" /> when there are no local latent classes. Called by "sim_no_local.m" and "sim_no_local_pi_subpop.m" files.
* `create_item_response_probs.m`: Helper function that generates item response probabilities <img src="https://latex.codecogs.com/gif.latex?\theta" />. Called by "sim_uneq.m" and "sim_uneq_pi_subpop.m" files.
* `create_item_response_probs_no_local.m`: Helper function that generates item response probabilities <img src="https://latex.codecogs.com/gif.latex?\theta" /> when there are no local latent classes. Called by "sim_no_local.m" and "sim_no_local_pi_subpop.m" files.
* `sim_batcharray.sh`: Bash script for creating an array of simulations using a cluster.
* `Cluster_SOP_2022_1_26.pdf`: Instructions for running simulations and models on a cluster.

### Directed Acyclic Graphs (DAGs) of Data-Generating Mechanisms
![DAGs](https://user-images.githubusercontent.com/33609713/152619726-0ffca12d-e5c4-4a76-ad23-03cca3185e6a.PNG)

