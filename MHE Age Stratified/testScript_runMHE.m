%% Moving Horizon Estimation - Test on age groups
clc
clear all
close all
addpath(genpath('casadi-3.6.5-linux64-matlab2018b'))
addpath('/Users/marcodelloro/Desktop/Pandemic-System-Identification/Reconstructed Datasets')
import casadi.*
set(0,'DefaultFigureWindowStyle','docked')

%%  Data Loading and Initialization

N_u40 = 23205875; % Total under40 pop
N_mid = 18142711; % Total 40-60 pop
N_old = 13408810; % Total 60-80 pop
N_ger = 3951057;  % Total 80+ pop
Npop = N_u40 + N_mid + N_old + N_ger; % Total Italy population 
N_sim = 140;       % MPC horizon length
N_mhe = 21;        % Estimation horizon (3 weeks)
Ts = 1;            % Integration time step  

S_data = readtable('Reconstructed_S.csv'); 
I_data = readtable('Infected_I.csv');     
D_data = readtable('Reconstructed_D.csv');     
T1_data = readtable('Reconstructed_T1.csv');
T2_data = readtable('Reconstructed_T2.csv');
H_data = readtable('Reconstructed_H.csv');
E_data = readtable('Reconstructed_E.csv');

c_struct.home =  [  19.7896   7.35764  1.43232   1.3;
             4.73514   5.62234  0.445607  0.44;
             3.26361   2.32613  4.3079    4.00; 
             3.00      2.00     4.00      3.00  ]; 

c_struct.schl =  [  24.9874    2.216       0.0889567   0.0389567;
             8.94216    1.1419      0.0793523   0.0293523;
             0.666502   0.0661153   0.054174    0.014174;
             0.0566502  0.0361153   0.014174    0.004174 ];


c_struct.work =  [  13.2266     8.39663   0.0329148    0.0229148;
             9.40362     8.72244   0.0357202    0.0157202;
             0.0930825   0.154919  0.000840544  0.00040544;
             0.030825    0.015491  0.00040544   0.00000001 ]; 


c_struct.othr =  [  44.1768     12.4868   4.61979   4.61979
             11.3875     8.96874   4.04791   4.04791
             6.62855     9.04265   8.95355   8.00000 
             3.62855     6.04265   9.00000   9.95355 ];

% c = c_home + c_schl + c_work + c_othr;       % How the contact matrices work with each other

ymeas.u40 = [ S_data.u40';  [I_data.u40(1) I_data.u40']./Npop; D_data.u40'./Npop; (T1_data.u40 + T2_data.u40)'./Npop; H_data.u40'./Npop; E_data.u40'./Npop ];  
ymeas.mid = [ S_data.mid';  [I_data.mid(1) I_data.mid']./Npop; D_data.mid'./Npop; (T1_data.mid + T2_data.mid)'./Npop; H_data.mid'./Npop; E_data.mid'./Npop ];  
ymeas.old = [ S_data.old';  [I_data.old(1) I_data.old']./Npop; D_data.old'./Npop; (T1_data.old + T2_data.mid)'./Npop; H_data.old'./Npop; E_data.old'./Npop ];  
ymeas.ger = [ S_data.ger';  [I_data.ger(1) I_data.ger']./Npop; D_data.ger'./Npop; (T1_data.ger + T2_data.mid)'./Npop; H_data.ger'./Npop; E_data.ger'./Npop ]; 


[e,dltp] = runMHE(N_mhe, N_sim, Ts, ymeas, c_struct);