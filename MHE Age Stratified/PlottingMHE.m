clc
clear all 
close all
load 'MHE Age Stratified'/AgeOpti-Results_newBayes.mat
addpath('/Users/marcodelloro/Desktop/Pandemic-System-Identification/Reconstructed Datasets/')
S_data = readtable('Reconstructed_S.csv'); 
I_raw = readtable('Infected_I.csv');     
D_data = readtable('Reconstructed_D.csv');     
T1_data = readtable('Reconstructed_T1.csv');
T2_data = readtable('Reconstructed_T2.csv');
H_data = readtable('Reconstructed_H.csv');
E_data = readtable('Reconstructed_E.csv');
N_u40 = 23205875; N_mid = 18142711; 
N_old = 13408810; N_ger = 3951057;  
N_mhe = 21;  N = 399; Npop = N_u40 + N_mid + N_old + N_ger; 
date = datetime('31-Aug-2020', 'InputFormat', 'dd-MMM-yyyy'):datetime('03-Oct-2021', 'InputFormat', 'dd-MMM-yyyy');

% fixing the values of data for the plots
I_data.u40 = [I_raw.u40(1) I_raw.u40']./Npop;         I_data.mid = [I_raw.mid(1) I_raw.mid']./Npop;
I_data.old = [I_raw.old(1) I_raw.old']./Npop;         I_data.ger = [I_raw.ger(1) I_raw.ger']./Npop;

T_data.u40 = (T1_data.u40  + T2_data.u40)'./Npop;       T_data.mid = (T1_data.mid  + T2_data.mid)'./Npop;     
T_data.old = (T1_data.old  + T2_data.u40)'./Npop;       T_data.ger = (T1_data.ger  + T2_data.ger)'./Npop;  

%% Plotting

policy_idx = [1 40 69 116 141 190 242 294 399];   % Values taken from italian policies applied for the pandemic

% for ii=1:length(policy_idx)
%     policy_dates(ii) = [date(policy_idx(ii))];
% end
% 
% customColors2 = {   [1, 0.6, 0.2],
%                     [1, 0.3, 0.05],
%                     [1, 0.2, 0.05],
%                     [0.8, 0.2, 0.1],
%                     [1, 0.2, 0.05],            
%                     [0.8, 0.2, 0.1],
%                     [1, 0.3, 0.05],
%                     [1, 0.6, 0.2]
%                 };
% 
% for ii = 1:length(policy_dates)-1
%     area.x(ii, :) = [policy_dates(ii) policy_dates(ii) policy_dates(ii+1) policy_dates(ii+1)];
%     area.y_alpha(ii, :) = [0 max(results.par.alpha)*1.5 max(results.par.alpha)*1.5 0];
% end

figure(1)
% Subplot 1: Under 40
subplot(2, 2, 1);
scatter(date(N_mhe:N-N_mhe), S_data.u40(N_mhe:N-N_mhe), 20, 'filled');
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("S-u40"), 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
sgtitle('\textbf{S - Susceptible}', 'Interpreter', 'latex');


% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
scatter(date(N_mhe:N-N_mhe), S_data.mid(N_mhe:N-N_mhe), 20, 'filled');
grid on, box on 
hold on
plot(date(N_mhe:N-N_mhe), results.sts.("S-mid"), 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
scatter(date(N_mhe:N-N_mhe), S_data.old(N_mhe:N-N_mhe), 20, 'filled');
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("S-old"), 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
scatter(date(N_mhe:N-N_mhe), S_data.ger(N_mhe:N-N_mhe), 20, 'filled');
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("S-ger"), 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('Real Data', 'Model', 'Interpreter', 'latex', 'location', 'northeast');
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
lgd = legend('Real Data', 'Model', 'Interpreter', 'latex', 'location', 'northeast');


figure(2)
% Subplot 1: Under 40
subplot(2, 2, 1);
scatter(date(N_mhe:N-N_mhe), I_data.u40(N_mhe:N-N_mhe), 20, 'filled');
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("I-u40"), 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,6e-3])
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
scatter(date(N_mhe:N-N_mhe), I_data.mid(N_mhe:N-N_mhe), 20, 'filled');
hold on 
plot(date(N_mhe:N-N_mhe), results.sts.("I-mid"), 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,5e-3])
hold on
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
scatter(date(N_mhe:N-N_mhe), I_data.old(N_mhe:N-N_mhe), 20, 'filled');
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("I-old"), 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
hold on
ylim([0,3e-3])
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
scatter(date(N_mhe:N-N_mhe), I_data.ger(N_mhe:N-N_mhe), 20, 'filled');
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("I-ger"), 'LineWidth', 1.5, 'MarkerSize', 5);
hold on
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,18e-4])
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('Real Data', 'Model', 'Interpreter', 'latex', 'location', 'northeast');
sgtitle('\textbf{I - Infected}', 'Interpreter', 'latex');



figure(3)
% Subplot 1: Under 40
subplot(2, 2, 1);
scatter(date(N_mhe:N-N_mhe), D_data.u40(N_mhe:N-N_mhe)./Npop, 20, 'filled');
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("D-u40"), 'LineWidth', 1.5, 'MarkerSize', 5);
hold on
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
scatter(date(N_mhe:N-N_mhe), D_data.mid(N_mhe:N-N_mhe)./Npop, 20, 'filled');
hold on 
plot(date(N_mhe:N-N_mhe), results.sts.("D-mid"), 'LineWidth', 1.5, 'MarkerSize', 5);
hold on
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
scatter(date(N_mhe:N-N_mhe), D_data.old(N_mhe:N-N_mhe)./Npop, 20, 'filled');
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("D-old"), 'LineWidth', 1.5, 'MarkerSize', 5);
hold on
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
scatter(date(N_mhe:N-N_mhe), D_data.ger(N_mhe:N-N_mhe)./Npop, 20, 'filled');
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("D-ger"), 'LineWidth', 1.5, 'MarkerSize', 5);
hold on
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('Real Data', 'Model', 'Interpreter', 'latex', 'location', 'northeast');
lgd.FontSize = 18; 
sgtitle('\textbf{D - Detected}', 'Interpreter', 'latex');

figure(4)
% Subplot 1: Under 40
subplot(2, 2, 1);
scatter(date(N_mhe:N-N_mhe), T_data.u40(N_mhe:N-N_mhe), 20, 'filled');
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("T-u40"), 'LineWidth', 1.5, 'MarkerSize', 5);
hold on
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
scatter(date(N_mhe:N-N_mhe), T_data.mid(N_mhe:N-N_mhe), 20, 'filled');
hold on 
plot(date(N_mhe:N-N_mhe), results.sts.("T-mid"), 'LineWidth', 1.5, 'MarkerSize', 5);
hold on
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
scatter(date(N_mhe:N-N_mhe), T_data.old(N_mhe:N-N_mhe), 20, 'filled');
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("T-old"), 'LineWidth', 1.5, 'MarkerSize', 5);
hold on
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
scatter(date(N_mhe:N-N_mhe), T_data.ger(N_mhe:N-N_mhe), 20, 'filled');
hold on, box on 
plot(date(N_mhe:N-N_mhe), results.sts.("T-ger"), 'LineWidth', 1.5, 'MarkerSize', 5);
hold on
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('Real Data', 'Model', 'Interpreter', 'latex', 'location', 'northeast');
sgtitle('\textbf{T - Hospitalised}', 'Interpreter', 'latex');


figure(5)
% Subplot 1: Under 40
subplot(2, 2, 1);
scatter(date(N_mhe:N-N_mhe), H_data.u40(N_mhe:N-N_mhe)./Npop, 20, 'filled');
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("H-u40"), 'LineWidth', 1.5, 'MarkerSize', 5);
hold on
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
scatter(date(N_mhe:N-N_mhe), H_data.mid(N_mhe:N-N_mhe)./Npop, 20, 'filled');
hold on 
plot(date(N_mhe:N-N_mhe), results.sts.("H-mid"), 'LineWidth', 1.5, 'MarkerSize', 5);
hold on
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
scatter(date(N_mhe:N-N_mhe), H_data.old(N_mhe:N-N_mhe)./Npop, 20, 'filled');
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("H-old"), 'LineWidth', 1.5, 'MarkerSize', 5);
hold on
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
scatter(date(N_mhe:N-N_mhe), H_data.ger(N_mhe:N-N_mhe)./Npop, 20, 'filled');
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("H-ger"), 'LineWidth', 1.5, 'MarkerSize', 5);
hold on
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('Real Data', 'Model', 'Interpreter', 'latex', 'location', 'northeast');
lgd.FontSize = 18; 

figure(6)
% Subplot 1: Under 40
subplot(2, 2, 1);
scatter(date(N_mhe:N-N_mhe), E_data.u40(N_mhe:N-N_mhe)./Npop, 20, 'filled');
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("E-u40"), 'LineWidth', 1.5, 'MarkerSize', 5);
hold on
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
scatter(date(N_mhe:N-N_mhe), E_data.mid(N_mhe:N-N_mhe)./Npop, 20, 'filled');
hold on 
plot(date(N_mhe:N-N_mhe), results.sts.("E-mid"), 'LineWidth', 1.5, 'MarkerSize', 5);
hold on
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
scatter(date(N_mhe:N-N_mhe), E_data.old(N_mhe:N-N_mhe)./Npop, 20, 'filled');
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("E-old"), 'LineWidth', 1.5, 'MarkerSize', 5);
hold on
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
scatter(date(N_mhe:N-N_mhe), E_data.ger(N_mhe:N-N_mhe)./Npop, 20, 'filled');
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("E-ger"), 'LineWidth', 1.5, 'MarkerSize', 5);
hold on
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('Real Data', 'Model', 'Interpreter', 'latex', 'location', 'northeast');
lgd.FontSize = 18; 
sgtitle('\textbf{E - Deceased}', 'Interpreter', 'latex');


%% 
figure(11)
% Subplot 1: Under 40
subplot(2, 2, 1);
plot(date(N_mhe:N-N_mhe), results.par.("alp-u40"), 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
plot(date(N_mhe:N-N_mhe), results.par.("alp-mid"), 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
plot(date(N_mhe:N-N_mhe), results.par.("alp-old"), 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
plot(date(N_mhe:N-N_mhe), results.par.("alp-ger"), 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
sgtitle('$\alpha$ - \textbf{Transmission rate}', 'Interpreter', 'latex');


figure(12)
% Subplot 1: Under 40
subplot(2, 2, 1);
plot(date(N_mhe:N-N_mhe), results.par.("gam-u40"), 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
plot(date(N_mhe:N-N_mhe), results.par.("gam-mid"), 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
plot(date(N_mhe:N-N_mhe), results.par.("gam-old"), 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
plot(date(N_mhe:N-N_mhe), results.par.("gam-ger"), 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
sgtitle('$\gamma$ - \textbf{Detection rate}', 'Interpreter', 'latex');



figure(13)
% Subplot 1: Under 40
subplot(2, 2, 1);
plot(date(N_mhe:N-N_mhe), results.par.("dlt-u40"), 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
plot(date(N_mhe:N-N_mhe), results.par.("dlt-mid"), 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
plot(date(N_mhe:N-N_mhe), results.par.("dlt-old"), 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
plot(date(N_mhe:N-N_mhe), results.par.("dlt-ger"), 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
sgtitle('$\delta$ - \textbf{Aggravation rate}', 'Interpreter', 'latex');

 
figure(14)
% Subplot 1: Under 40
subplot(2, 2, 1);
plot(date(N_mhe:N-N_mhe), results.par.("tau-u40"), 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
plot(date(N_mhe:N-N_mhe), results.par.("tau-mid"), 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
plot(date(N_mhe:N-N_mhe), results.par.("tau-old"), 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
plot(date(N_mhe:N-N_mhe), results.par.("tau-ger"), 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
sgtitle('$\tau$ - \textbf{Mortality rate}', 'Interpreter', 'latex');



figure(15)
% Subplot 1: Under 40
subplot(2, 2, 1);
plot(date(N_mhe:N-N_mhe), results.par.("sgm-u40"), 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
plot(date(N_mhe:N-N_mhe), results.par.("sgm-mid"), 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
plot(date(N_mhe:N-N_mhe), results.par.("sgm-old"), 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
plot(date(N_mhe:N-N_mhe), results.par.("sgm-ger"), 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
sgtitle('$\sigma$ - \textbf{(Hospital) Healing rate}', 'Interpreter', 'latex');
