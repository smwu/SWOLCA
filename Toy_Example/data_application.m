% NHANES Data Pre-processing
% Prepare diet and CVD data for low-income women
% Inputs: 'nhanesallage_frscores1118.csv' file
% Outputs: Structural array named 'samp_data' with the following columns:
%   X_data: food intake data as a matrix
%   Y_data: outcome data as a vector
%   sample_wt: survey weights as a vector

%% Load data
clear; clc; 
nhanes = readtable('nhanesallage_frscores1118.csv');

%% Subsetting 
% Select for women
nhanes_f = nhanes(strcmp(nhanes.sex, 'female'), :);
% Select for low-income
nhanes_f_low = nhanes_f(strcmp(nhanes_f.poverty, 'At or Below'), :);

%% Select key variables
samp_data.X_data = table2array(nhanes_f_low(:, 57:85));
samp_data.Y_data = nhanes_f_low(:, 'bp_flag');
samp_data.sample_wt = nhanes_f_low(:, 'dietwt8yr');
samp_data.true_Si = 
