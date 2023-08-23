# Simulation Code
This folder includes the `simulate_data.R` file, which contains code used to generate simulated population and sample data.

### Simulated population of size N=80000
Simulated population scenarios are specified with a four-digit number. The first digit specifies the strength of the pattern (mode 85\% or 55\%,). The second digit specifies the dietary pattern (disjoint or overlapping). The third digit specifies the data-generating model (stratum as confounder, stratum as effect modifier, or additional confounders). The fourth digit specifies clustering in the outcome (no clustering or 1000 culsters of size 80). 

### Simulated samples
100 independent samples are drawn from the population for each sample scenario, specified with a six-digit number. The first four digits carry over from the population scenario. The fifth digit specifies sample size (1\%, 5\%, or 10%). The sixth digit specifies sampling design (simple random sample, stratified sampling, or stratified cluster sampling).

