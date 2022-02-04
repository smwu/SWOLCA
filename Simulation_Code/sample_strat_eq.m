
% sample_strat_eq reads in a simulated population dataset and takes a 
% simple random sample, then saves the sample. Corresponds to Scenarios 13-16.
function sample_strat_eq(scenario, sim_n, samp_n)
    %% Specifications
    rng(samp_n, 'twister');                        % Set seed
    in_dir = "/n/home01/stephwu18/wsOFMM/data/";   % Input directory
    out_dir = "/n/home01/stephwu18/wsOFMM/data/";  % Output directory
    in_dir = strcat(pwd, "\Data\");
    out_dir = strcat(pwd, "\Data\");   
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
    
    n_s = [100, 100, 100, 100];               % Sample sizes for each subpop
    % Obtain sample
    sample_data = sample_indivs(N_s, n_s, S, true, sim_data.true_Si, sim_data.true_Ci, sim_data.true_Li, sim_data.X_data, sim_data.Y_data, sim_data.true_Phi, sim_data.true_K, sim_data);
    sim_data = sample_data; 
    % Save sample data
    save(strcat(out_dir, 'simdata_scen', num2str(scenario + 12),'_iter', num2str(sim_n), '_samp', num2str(samp_n)), 'sim_data');
end  