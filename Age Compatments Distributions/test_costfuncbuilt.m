clc
clear all 
close all
load 'MHE Age Stratified'/AgeOpti-Results2.mat
addpath '/Users/marcodelloro/Desktop/Pandemic-System-Identification'/'Reconstructed Datasets'/
load 'MHE Age Stratified'/AgeOpti-Results2.mat

N_u40 = 23205875; % Total under40 pop
N_mid = 18142711; % Total 40-60 pop
N_old = 13408810; % Total 60-80 pop
N_ger = 3951057;  % Total 80+ pop
Npop = N_u40 + N_mid + N_old + N_ger; % Total Italy population 
T_sim = 147;       % MPC horizon length (NOTE: This must be a multiplle integer of N_mhe)
N_mhe = 21;        % Estimation horizon (3 weeks)
Ts = 1;            % Integration time step  

S_data = readtable('Reconstructed_S.csv'); 
I_data = readtable('Infected_I.csv');     
D_data = readtable('Reconstructed_D.csv');     
T1_data = readtable('Reconstructed_T1.csv');
T2_data = readtable('Reconstructed_T2.csv');
H_data = readtable('Reconstructed_H.csv');
E_data = readtable('Reconstructed_E.csv');

ymeas.u40 = [ S_data.u40';  [I_data.u40(1) I_data.u40']./Npop; D_data.u40'./Npop; (T1_data.u40 + T2_data.u40)'./Npop; H_data.u40'./Npop; E_data.u40'./Npop ];  
ymeas.mid = [ S_data.mid';  [I_data.mid(1) I_data.mid']./Npop; D_data.mid'./Npop; (T1_data.mid + T2_data.mid)'./Npop; H_data.mid'./Npop; E_data.mid'./Npop ];  
ymeas.old = [ S_data.old';  [I_data.old(1) I_data.old']./Npop; D_data.old'./Npop; (T1_data.old + T2_data.mid)'./Npop; H_data.old'./Npop; E_data.old'./Npop ];  
ymeas.ger = [ S_data.ger';  [I_data.ger(1) I_data.ger']./Npop; D_data.ger'./Npop; (T1_data.ger + T2_data.mid)'./Npop; H_data.ger'./Npop; E_data.ger'./Npop ]; 


% modification of the 'y_meas' struct to be used into the obj function
y_meas = [ymeas.u40' ymeas.mid' ymeas.old' ymeas.ger'];
y_meas = y_meas(N_mhe:N_mhe+T_sim-1,:);
sts_mat = table2array(results.sts(1:T_sim,:));

%% CONSTRUCTION OF THE COST FUNCTION

e = abs(y_meas - sts_mat); % residuals calculation (estimation error)
r = []; 
N = length(y_meas);

alpha = 0.05;                       % significance level (95%)
df = N;                             % Degrees of freedom
z_alpha = chi2inv(1 - alpha, df);   % Chi-squared threshold value


%                      ----- Construction of estimated covariance function  ----- %

for tau = 1:N
    r_hat = 0;  % Initialize sum for the current tau "r_hat" value
    for t = 1:N - tau
        r_hat = r_hat + e(t + tau,:).* e(t,:);
    end
    r(tau,:) = r_hat ./ N;  % Store autocorrelation value for current tau
end

r_hat0 = sum(e.^2, 1)./N;      % value of the autocorrelation when tau = 0


%                     ----- Construction of test quantity NrTr/r_hat(0)^2 (Distributed as Chi-squared)----- %

testq = N .* sum(r.^2, 1)./(r_hat0).^2;

%                     ----- Construction of Cost Function For Bayes ----- %


