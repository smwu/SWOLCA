# Results Summarized over Iterations and Samples

This folder contains summary output calculated from combining results from multiple iterations or samples. It includes the following files:
* `diagnostics.R`: R code for obtaining graphics using the summarized results for multiple scenarios.
* `testing_code.m`: MATLAB code for troubleshooting and testing functions
* `sim_summary_wsOFMM.m`: MATLAB code to create summary output by combining results from multiple iterations or samples . 
  * Input: files from the "wsOFMM/Results" folder. If population data was used (Scenarios 1-4), all files must have the same scenario, but the iteration index may change. If sample data was used (Scenarios 5-16), all files must have the same scenario and iteration but the sample number index may change. 
  * Output:  (See next bullet point.) If population data was used (Scenarios 1-4), a "summary_...\_pop.mat" file. If sample data was used (Scenarios 5-16), a "summary_...\_sample.mat" file.
* `summary_....mat` files: Summary output from the "sim_summary_wsOFMM.m" file. 
  * File names are created using the following ordered components: 
     1. `summary_wsOFMM` or `summary_sOFMM`: Indicates if the wsOFMM or sOFMM model was used.
     3. `_scen<scenario number>`: Indicates the sampling design scenario for the simulation. There are 16 possible scenarios, described in the "wsOFMM/Simulation_Code" README file.
     5. `_pop` or `_sample`: Indicates whether sample or population data was used. For population data (Scenarios 1-4), output is summarized over multiple iterations. For sample data (Scenarios 5-16), output is summarized over multiple samples. 
     6. `.mat`: Indicates dataset type as MATLAB data.
