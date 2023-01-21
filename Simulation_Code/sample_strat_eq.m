
% sample_strat_eq reads in a simulated population dataset and takes a 
% stratified random sample with equal numbers sampled from each 
% subpopulation, then saves the sample. Corresponds to Scenarios 13-16.
% Inputs:
%   sim_n: Simulation iteration number
%   scenario: Population simulation scenario (1, 2, 3, or 4)
%   samp_n: Sample iteration number. Default is 1, so one sample is taken per iteration
% Outputs a simulated sample dataset
function sample_strat_eq(sim_n, scenario, samp_n)
    %% Specifications
    rng(samp_n, 'twister');                        
    in_dir = "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/Data/";   % Input directory
    out_dir = "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/Data/";  % Output directory 
%          in_dir = strcat(pwd, "/");
%          out_dir = strcat(pwd, "/");    
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
    
    n_s = [2000, 2000];               % Sample sizes for each subpop
    % Obtain sample
    sample_data = sample_indivs(N_s, n_s, S, true, sim_data.true_Si, sim_data.true_Ci, sim_data.true_Li, sim_data.X_data, sim_data.Y_data, sim_data.true_Phi, sim_data.true_K, sim_data);
    sim_data = sample_data; 
    % Save sample data
    save(strcat(out_dir, 'simdata_scen', num2str(scenario + 200),'_iter', num2str(sim_n), '_samp', num2str(samp_n)), 'sim_data');
end  

% Testing Code
%     in_dir = "C:/Users/Lang/Documents/Harvard/Research/Briana/supRPC/wsOFMM/Toy_Example/";   % Input directory
%     out_dir = "C:/Users/Lang/Documents/Harvard/Research/Briana/supRPC/wsOFMM/Toy_Example/";  % Output directory 
%   
% for samp_n = 1:20
%     sample_strat_eq(20, 1, samp_n);
% end  
