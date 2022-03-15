
% sample_SRS reads in a simulated population dataset and takes a simple
% random sample, then saves the sample. Corresponds to Scenarios 5-8.
% Inputs:
%   scenario: Population simulation scenario (1, 2, 3, or 4)
%   sim_n: Simulation iteration number
%   samp_n: Sample iteration number. Default is 1, so one sample is taken per iteration
% Outputs a simulated sample dataset
function sample_SRS(scenario, sim_n, samp_n)
    %% Specifications
    if samp_n >1  % If more than one sample per iteration, set seed based on sample
        rng(samp_n, 'twister');                        
    else          % If only one sample per iteration, set seed based on iteration
        rng(sim_n, 'twister');
    end    
    in_dir = "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/Toy_Example/";   % Input directory
    out_dir = "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/Toy_Example/";  % Output directory
    % Load population simulated dataset
    load(strcat(in_dir, 'simdata_scen', num2str(scenario), '_iter', num2str(sim_n), '.mat'), 'sim_data')  
    
    %% Create and save sample
    S = length(unique(sim_data.true_Si));     % Number of subpopulations
    N_s = zeros(1, S);                        % Subpopulation sizes
    for s = 1:S
        N_s(s) = sum(sim_data.sample_wt(sim_data.true_Si == s));
    end    
    N = sum(N_s);                             % Population size
    if exist('sim_data.true_Li', 'var') == 0  % Check existence of local classes
        sim_data.true_Li = zeros(N, 1);
    end
    
    n_s = 400;                                % Sample size
    % Obtain sample
    sample_data = sample_indivs(N, n_s, S, false, sim_data.true_Si, sim_data.true_Ci, sim_data.true_Li, sim_data.X_data, sim_data.Y_data, sim_data.true_Phi, sim_data.true_K, sim_data);
    sim_data = sample_data; 
    % Save sample data
    save(strcat(out_dir, 'simdata_scen', num2str(scenario + 4),'_iter', num2str(sim_n), '_samp', num2str(samp_n)), 'sim_data');
end

  