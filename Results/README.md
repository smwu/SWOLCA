# Results from Model Runs

This folder contains output from running the wsOFMM or sOFMM model on the simulated datasets. Output is either **MCMC** output generated right after running the MCMC sampling algorithm, or **results** output generated after post-processing and including posterior estimates. 

File names are created using the following ordered components: 
 1. `wsOFMM` or `sOFMM`: Indicates if the wsOFMM or sOFMM model was used.
 2. `_MCMC` or `_results`: Indicates if output is MCMC output or posterior results output.
 3. `_scen<scenario number>`: Indicates the sampling design scenario for the simulation. There are 16 possible scenarios, described in the "wsOFMM/Simulation_Code" README file.
 4. `_iter<iteration number>`: Indicates the iteration number for the simulation. Iterations are generated independently of each other, using different seeds for random number generation.
 5. (OPTIONAL) `_samp<sample number>`: This will only appear for sample data, not for population data. For a sample taken from simulated population, the sample number indexes the sample drawn. Samples are drawn independently of each other, with different random number generation seeds.
 6. `.mat`: Indicates dataset type as MATLAB data.
