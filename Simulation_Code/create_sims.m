%% Create simulations
scen = 2;            % Data-generating scenario
sim_n = 1;           % Only use 1 population iteration
num_samples = 100;   % Number of samples

% Create population data
sim_confounder(sim_n, scen);

% Create sample data
for samp_n = 1:num_samples
    % SRS
    sample_SRS(sim_n, scen, samp_n);
    % Stratified random sample with equal allocation per stratum
    sample_strat_eq(sim_n, scen, samp_n);
end 