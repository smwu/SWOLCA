# Analysis Code
This folder includes the `summarize_results.R` file for summarizing model output into tables and figures, and the `summarize_results_functions.R` file that includes helper functions. 

Results are summarized across 100 independent iterations, with model performance assessed using mean absolute bias, width of 95\% credible interval, and coverage of 95\% credible interval. Results for each parameter are averaged across the dimensions of the parameter.

To display results, functions are available to create tables comparing the bias, variance, and coverage of parameter estimates for the SWOLCA, SOLCA, and WOLCA models for various scenarios. Additional code is provided for creating boxplots of RMSE, bias, and parameter estimation for each iteration, as well as for plotting coverage, modal exposure category probabilities, and outcome regression probabilities. 

### `Old_Code` folder
This folder includes the `summarize_results_new.R` file, which contains old code that summarizes results for various extraneous experiments.
