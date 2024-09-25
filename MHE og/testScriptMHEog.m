%% Moving Horizon Estimation
clc
clear all
close all
addpath('/Users/marcodelloro/Downloads/casadi-3.6.3-osx64-matlab2018b')
import casadi.*
load('SIDTTHE_data_DEF.mat');
set(0,'DefaultFigureWindowStyle','docked')

%%  Data Loading and Initialization
tic

Npop = 59240329; % Total Population of Italy
N = 399;         % MPC horizon length
N_mhe = 21;      % Estimation horizon (3 weeks)
date = SIDTTHE_data{1,1}.date;

I_data = SIDTTHE_data{1,1}.data / Npop;     
D_data = SIDTTHE_data{2,1}.data / Npop;     
T_data = (SIDTTHE_data{3,1}.data + SIDTTHE_data{4,1}.data) / Npop;
H_data = SIDTTHE_data{7,1}.data / Npop;
H_dataAug = SIDTTHE_data{6,1}.data / Npop;
E_data = SIDTTHE_data{5,1}.data  / Npop;
S_data = ones(length(I_data),1)' - (I_data + D_data + T_data + H_dataAug + E_data );

ymeas = [S_data; I_data; D_data; T_data; H_dataAug; E_data]; % creation of the measurement vector
n_meas = size(ymeas,1);

T_sim = 70;
%% Bayesian optimization of Hyperparams workflow
% see the BO steps @ https://nl.mathworks.com/help/stats/bayesian-optimization-workflow.html

%                         ----- Step 1 - OPTIMIZATION VARIABLES PREPARATION ----- %
% For each variable in my objective function, I create a variable for the BO
%
%   Z1  --> 6 elements vector  ( weights on the difference between consecutive sts )
%   Z2  --> 5 elements vector  ( weights on the difference between consecutive param. estiamte )
%   Z3  --> 6 elements vector  ( weights on the process noise )
%

% Z1 vector elements 
Z1_1 = optimizableVariable('z1_1', [1, 1e4], 'Type', 'integer');
Z1_2 = optimizableVariable('z1_2', [1, 1e4], 'Type', 'integer');
Z1_3 = optimizableVariable('z1_3', [1, 1e4], 'Type', 'integer');

% Z2 vector elements 
Z2_1 = optimizableVariable('z2_1', [1, 1e4], 'Type', 'integer');
Z2_2 = optimizableVariable('z2_2', [1, 1e4], 'Type', 'integer');
Z2_3 = optimizableVariable('z2_3', [1, 1e4], 'Type', 'integer');
Z2_4 = optimizableVariable('z2_4', [1, 1e4], 'Type', 'integer');
Z2_5 = optimizableVariable('z2_5', [1, 1e4], 'Type', 'integer');

% Z3 vector elements 
Z3_1  = optimizableVariable('z3_1', [1, 1e4], 'Type', 'integer');
Z3_2  = optimizableVariable('z3_2', [1, 1e4], 'Type', 'integer'); 
Z3_3  = optimizableVariable('z3_3', [1, 1e4], 'Type', 'integer');

weights = [ Z1_1, Z1_2, Z1_3,...
            Z2_1, Z2_2, Z2_3, Z2_4, Z2_5,...
            Z3_1, Z3_2, Z3_3];


%                      ----- Step 2 - Built of the objective function   ----- %

objectiveFcn = @(weights) bayesObjMHEog(weights, N_mhe, T_sim, ymeas);

%                      ----- Step 3 - Bayesian optimization and results ----- %

results = bayesopt(objectiveFcn, weights, ...
                   'Verbose', 0, ...
                   'MaxObjectiveEvaluations', 30);

optResults = struct();
optResults.fullResults = results;

% Save the struct to a .mat file if needed
save('BayesResult_ThursNight.mat', 'optResults');  

toc
