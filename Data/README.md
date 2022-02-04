# Simulated Data

This folder contains simulated datasets that can be used when running the sOFMM and wsOFMM models. Datasets can either include **population data**, where observations are simulated for an entire population of size N=8000, or include **sample data**, where a sample of size n=400 is drawn from one of the simulated population datasets using a particular sampling scheme such as stratified sampling. 

Dataset names are created using the following ordered components: 
 1. `simdata_scen<scenario number>`: Indicates the sampling design scenario for the simulation. There are 16 possible scenarios, described in the "wsOFMM/Simulation_Code" README file.
 2. `\_iter<iteration number>`: Indicates the iteration number for the simulation. Iterations are generated independently of each other, using different seeds for random number generation.
 3. (OPTIONAL) `\_samp<sample number>`: This will only appear for sample data, not for population data. For a sample taken from simulated population, the sample number indexes the sample drawn. Samples are drawn independently of each other, with different random number generation seeds.
 4. `.mat`: Indicates dataset type as MATLAB data.
