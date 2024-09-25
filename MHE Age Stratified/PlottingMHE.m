clc
clear all 
close all
load '/Users/marcodelloro/Desktop'/SatBayesAGE_ResWeightsUppati.mat
load SmoothVariantasIta.mat
load '/Users/marcodelloro/Desktop/Pandemic-System-Identification/MHE Age Stratified/FinalResults'/AgeOpti-ResBleau.mat
addpath('/Users/marcodelloro/Desktop/Pandemic-System-Identification/Reconstructed Datasets/')
set(0,'DefaultFigureWindowStyle','docked')
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

D_data.u40 = D_data.u40 ./Npop;       D_data.mid = D_data.mid ./Npop;      
D_data.old = D_data.old ./Npop;       D_data.ger = D_data.ger ./Npop;  

H_data.u40 = H_data.u40 ./Npop;       H_data.mid = H_data.mid ./Npop;      
H_data.old = H_data.old ./Npop;       H_data.ger = H_data.ger ./Npop;  

E_data.u40 = E_data.u40 ./Npop;       E_data.mid = E_data.mid ./Npop;      
E_data.old = E_data.old ./Npop;       E_data.ger = E_data.ger ./Npop;  

%% Plotting Settings

policy_idx = [1 40 69 116 141 190 242 294 399];   % Values taken from italian policies applied for the pandemic

for ii=1:length(policy_idx)
    policy_dates(ii) = [date(policy_idx(ii))];
end

customColors2 = {   [1, 0.6, 0.2],
                    [1, 0.3, 0.05],
                    [1, 0.2, 0.05],
                    [0.8, 0.2, 0.1],
                    [1, 0.2, 0.05],            
                    [0.8, 0.2, 0.1],
                    [1, 0.3, 0.05],
                    [1, 0.6, 0.2]
                };

for ii = 1:length(policy_dates)-1
    areaNPI.x(ii, :) = [policy_dates(ii) policy_dates(ii) policy_dates(ii+1) policy_dates(ii+1)];
    areaNPI.y_alpha(ii, :) = [0 max(results.par.("alp-u40"))*1.5 max(results.par.("alp-u40"))*1.5 0];
end


% Choice of the color to plot the variants
colors = {[0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250], [ 0.4940, 0.1840, 0.5560], [0.4660, 0.6740, 0.1880]};

%% I Compartment ALL plots

figure(2)
% Subplot 1: Under 40
subplot(2, 2, 1);
scatter(date(N_mhe:N-N_mhe), I_data.u40(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("I-u40"), 'LineWidth', 2, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,max(results.sts.("I-u40"))*1.05])
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
scatter(date(N_mhe:N-N_mhe), I_data.mid(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on 
plot(date(N_mhe:N-N_mhe), results.sts.("I-mid"), 'LineWidth', 2, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,max(results.sts.("I-mid"))*1.05])
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
scatter(date(N_mhe:N-N_mhe), I_data.old(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("I-old"), 'LineWidth', 2, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,max(results.sts.("I-old"))*1.05])
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
scatter(date(N_mhe:N-N_mhe), I_data.ger(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("I-ger"), 'LineWidth', 2, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,max(results.sts.("I-ger"))*1.05])
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('Real Data', 'Model', 'Interpreter', 'latex', 'location', 'northeast');
sgtitle('\textbf{I - Infected}', 'Interpreter', 'latex');


%                                        ---- Plot with NPIs ----                  %

figure(3)
% Subplot 1: Under 40
subplot(2, 2, 1);
scatter(date(N_mhe:N-N_mhe), I_data.u40(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("I-u40"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,max(results.sts.("I-u40"))*1.05])
hold on 
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
scatter(date(N_mhe:N-N_mhe), I_data.mid(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on 
plot(date(N_mhe:N-N_mhe), results.sts.("I-mid"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,max(results.sts.("I-mid"))*1.05])
hold on
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
scatter(date(N_mhe:N-N_mhe), I_data.old(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("I-old"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,3e-3])
hold on
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on, box on 
ylim([0,max(results.sts.("I-old"))*1.05])
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
scatter(date(N_mhe:N-N_mhe), I_data.ger(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("I-ger"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,max(results.sts.("I-ger"))*1.05])
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('Real Data', 'Model', 'Interpreter', 'latex', 'location', 'northeast');
sgtitle('\textbf{I - Infected}', 'Interpreter', 'latex');


%                                        ---- Plot with variants ----                  %
figure(4)
% Subplot 1: Under 40
subplot(2, 2, 1);
scatter(date(N_mhe:N-N_mhe), I_data.u40(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5,'HandleVisibility', 'off');
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("I-u40"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0], 'HandleVisibility', 'off');
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,max(results.sts.("I-u40"))*1.05])
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
ylabel('Variant percentage $\%$', 'Interpreter', 'latex');
ylim([0, 100]);
yyaxis left
ylabel('I-data u40', 'Interpreter', 'latex');
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2, 'HandleVisibility', 'off');
end
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('SARS CoV-2', 'Alpha var.', 'Gamma var.', 'Delta var.', 'Interpreter', 'latex', 'Location', 'northeast');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
scatter(date(N_mhe:N-N_mhe), I_data.mid(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("I-mid"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,max(results.sts.("I-mid"))*1.05])
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('I-data mid', 'Interpreter', 'latex');
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
scatter(date(N_mhe:N-N_mhe), I_data.old(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("I-old"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,max(results.sts.("I-old"))*1.05])
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('I-data old','Interpreter', 'latex');  % Label for the primary y-axis
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on;
box on;
ylim([0,max(results.sts.("I-old"))*1.05])
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
scatter(date(N_mhe:N-N_mhe), I_data.ger(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("I-ger"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,max(results.sts.("I-ger"))*1.05])
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('I-data ger','Interpreter', 'latex');  % Label for the primary y-axis
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
sgtitle('\textbf{I - Infected}', 'Interpreter', 'latex');

%% D Compartment ALL plots

figure(5)
% Subplot 1: Under 40
subplot(2, 2, 1);
scatter(date(N_mhe:N-N_mhe), D_data.u40(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("D-u40"), 'LineWidth', 2, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,6e-3])
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
scatter(date(N_mhe:N-N_mhe), D_data.mid(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on 
plot(date(N_mhe:N-N_mhe), results.sts.("D-mid"), 'LineWidth', 2, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,5e-3])
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
scatter(date(N_mhe:N-N_mhe), D_data.old(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("D-old"), 'LineWidth', 2, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,3e-3])
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
scatter(date(N_mhe:N-N_mhe), D_data.ger(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("D-ger"), 'LineWidth', 2, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,18e-4])
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('Real Data', 'Model', 'Interpreter', 'latex', 'location', 'northeast');
sgtitle('\textbf{D - Detected}', 'Interpreter', 'latex');


%                                        ---- Plot with NPIs ----                  %

figure(6)
% Subplot 1: Under 40
subplot(2, 2, 1);
scatter(date(N_mhe:N-N_mhe), D_data.u40(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("D-u40"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,6e-3])
hold on 
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
scatter(date(N_mhe:N-N_mhe), D_data.mid(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on 
plot(date(N_mhe:N-N_mhe), results.sts.("D-mid"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,5e-3])
hold on
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
scatter(date(N_mhe:N-N_mhe), D_data.old(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("D-old"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,3e-3])
hold on
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on, box on 
ylim([0,30e-4])
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
scatter(date(N_mhe:N-N_mhe), D_data.ger(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("D-ger"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,18e-4])
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('Real Data', 'Model', 'Interpreter', 'latex', 'location', 'northeast');
sgtitle('\textbf{D - Detected}', 'Interpreter', 'latex');


%                                        ---- Plot with variants ----                  %

figure(7)
% Subplot 1: Under 40
subplot(2, 2, 1);
scatter(date(N_mhe:N-N_mhe), D_data.u40(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5,'HandleVisibility', 'off');
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("D-u40"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0], 'HandleVisibility', 'off');
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0, 6e-3]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
ylabel('Variant percentage $\%$', 'Interpreter', 'latex');
ylim([0, 100]);
yyaxis left
ylabel('D-data u40', 'Interpreter', 'latex');
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2, 'HandleVisibility', 'off');
end
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('SARS CoV-2', 'Alpha var.', 'Gamma var.', 'Delta var.', 'Interpreter', 'latex', 'Location', 'northeast');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
scatter(date(N_mhe:N-N_mhe), D_data.mid(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("D-mid"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0, 5e-3]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('D-data mid', 'Interpreter', 'latex');
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
scatter(date(N_mhe:N-N_mhe), D_data.old(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("D-old"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0, 3e-3]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('D-data old','Interpreter', 'latex');  % Label for the primary y-axis
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on;
box on;
ylim([0, 30e-4]);
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
scatter(date(N_mhe:N-N_mhe), D_data.ger(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on
plot(date(N_mhe:N-N_mhe), results.sts.("D-ger"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0, 18e-4]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('D-data ger','Interpreter', 'latex');  % Label for the primary y-axis
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
sgtitle('\textbf{D - Detected}', 'Interpreter', 'latex');

%% T Compartment ALL plots

figure(8)
% Subplot 1: Under 40
subplot(2, 2, 1);
scatter(date(N_mhe:N-N_mhe), T_data.u40(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("T-u40"), 'LineWidth', 2, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
scatter(date(N_mhe:N-N_mhe), T_data.mid(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on 
plot(date(N_mhe:N-N_mhe), results.sts.("T-mid"), 'LineWidth', 2, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
scatter(date(N_mhe:N-N_mhe), T_data.old(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("T-old"), 'LineWidth', 2, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
scatter(date(N_mhe:N-N_mhe), T_data.ger(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("T-ger"), 'LineWidth', 2, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('Real Data', 'Model', 'Interpreter', 'latex', 'location', 'northeast');
sgtitle('\textbf{T - Hospitalised}', 'Interpreter', 'latex');


%                                        ---- Plot with NPIs ----                  %

figure(9)
% Subplot 1: Under 40
subplot(2, 2, 1);
scatter(date(N_mhe:N-N_mhe), T_data.u40(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("T-u40"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
hold on 
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on, box on 
ylim([0, max(T_data.u40)*1.2]);
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
scatter(date(N_mhe:N-N_mhe), T_data.mid(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on 
plot(date(N_mhe:N-N_mhe), results.sts.("T-mid"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
hold on
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
ylim([0, max(T_data.mid)*1.2]);
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
scatter(date(N_mhe:N-N_mhe), T_data.old(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("T-old"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
hold on
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on, box on 
ylim([0, max(T_data.old)*1.2]);
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
scatter(date(N_mhe:N-N_mhe), T_data.ger(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("T-ger"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on, box on 
ylim([0, max(T_data.ger)*1.2]);
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('Real Data', 'Model', 'Interpreter', 'latex', 'location', 'northeast');
sgtitle('\textbf{T - Hospitalised}', 'Interpreter', 'latex');


%                                        ---- Plot with variants ----                  %

figure(10)
% Subplot 1: Under 40
subplot(2, 2, 1);
scatter(date(N_mhe:N-N_mhe), T_data.u40(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5,'HandleVisibility', 'off');
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("T-u40"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0], 'HandleVisibility', 'off');
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
ylabel('Variant percentage $\%$', 'Interpreter', 'latex');
ylim([0, 100]);
yyaxis left
ylabel('T-data u40', 'Interpreter', 'latex');
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2, 'HandleVisibility', 'off');
end
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('SARS CoV-2', 'Alpha var.', 'Gamma var.', 'Delta var.', 'Interpreter', 'latex', 'Location', 'northeast');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
scatter(date(N_mhe:N-N_mhe), T_data.mid(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("T-mid"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('T-data mid', 'Interpreter', 'latex');
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
scatter(date(N_mhe:N-N_mhe), T_data.old(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("T-old"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('T-data old','Interpreter', 'latex');  % Label for the primary y-axis
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
scatter(date(N_mhe:N-N_mhe), T_data.ger(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on
plot(date(N_mhe:N-N_mhe), results.sts.("T-ger"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('T-data ger','Interpreter', 'latex');  % Label for the primary y-axis
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
sgtitle('\textbf{T - Hospitalised}', 'Interpreter', 'latex');

%% E Compartment ALL plots

figure(11)
% Subplot 1: Under 40
subplot(2, 2, 1);
scatter(date(N_mhe:N-N_mhe), E_data.u40(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("E-u40"), 'LineWidth', 2, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
scatter(date(N_mhe:N-N_mhe), E_data.mid(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on 
plot(date(N_mhe:N-N_mhe), results.sts.("E-mid"), 'LineWidth', 2, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
scatter(date(N_mhe:N-N_mhe), E_data.old(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("E-old"), 'LineWidth', 2, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
scatter(date(N_mhe:N-N_mhe), E_data.ger(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("E-ger"), 'LineWidth', 2, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('Real Data', 'Model', 'Interpreter', 'latex', 'location', 'northeast');
sgtitle('\textbf{E - Deceased}', 'Interpreter', 'latex');


%                                        ---- Plot with NPIs ----                  %

figure(12)
% Subplot 1: Under 40
subplot(2, 2, 1);
scatter(date(N_mhe:N-N_mhe), E_data.u40(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("E-u40"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
hold on 
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on, box on 
ylim([0, max(E_data.u40)*1.2]);
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
scatter(date(N_mhe:N-N_mhe), E_data.mid(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on 
plot(date(N_mhe:N-N_mhe), results.sts.("E-mid"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
hold on
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
ylim([0, max(E_data.mid)*1.2]);
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
scatter(date(N_mhe:N-N_mhe), E_data.old(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("E-old"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
hold on
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on, box on 
ylim([0, max(E_data.old)*1.2]);
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
scatter(date(N_mhe:N-N_mhe), E_data.ger(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("E-ger"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on, box on 
ylim([0, max(E_data.ger)*1.2]);
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('Real Data', 'Model', 'Interpreter', 'latex', 'location', 'northeast');
sgtitle('\textbf{E - Deceased}', 'Interpreter', 'latex');


%                                        ---- Plot with variants ----                  %

figure(13)
% Subplot 1: Under 40
subplot(2, 2, 1);
scatter(date(N_mhe:N-N_mhe), E_data.u40(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5,'HandleVisibility', 'off');
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("E-u40"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0], 'HandleVisibility', 'off');
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
ylabel('Variant percentage $\%$', 'Interpreter', 'latex');
ylim([0, 100]);
yyaxis left
ylabel('E-data u40', 'Interpreter', 'latex');
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2, 'HandleVisibility', 'off');
end
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('SARS CoV-2', 'Alpha var.', 'Gamma var.', 'Delta var.', 'Interpreter', 'latex', 'Location', 'northeast');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
scatter(date(N_mhe:N-N_mhe), E_data.mid(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("E-mid"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('E-data mid', 'Interpreter', 'latex');
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
scatter(date(N_mhe:N-N_mhe), E_data.old(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("E-old"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('E-data old','Interpreter', 'latex');  % Label for the primary y-axis
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
scatter(date(N_mhe:N-N_mhe), E_data.ger(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on
plot(date(N_mhe:N-N_mhe), results.sts.("E-ger"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('E-data ger','Interpreter', 'latex');  % Label for the primary y-axis
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
sgtitle('\textbf{E - Deceased}', 'Interpreter', 'latex');

%% H Compartment ALL plots

figure(50)
% Subplot 1: Under 40
subplot(2, 2, 1);
scatter(date(N_mhe:N-N_mhe), H_data.u40(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("H-u40"), 'LineWidth', 2, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
scatter(date(N_mhe:N-N_mhe), H_data.mid(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on 
plot(date(N_mhe:N-N_mhe), results.sts.("H-mid"), 'LineWidth', 2, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
scatter(date(N_mhe:N-N_mhe), H_data.old(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("H-old"), 'LineWidth', 2, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
scatter(date(N_mhe:N-N_mhe), H_data.ger(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("H-ger"), 'LineWidth', 2, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('Real Data', 'Model', 'Interpreter', 'latex', 'location', 'northeast');
sgtitle('\textbf{H - Healed}', 'Interpreter', 'latex');


%                                        ---- Plot with NPIs ----                  %

figure(51)
% Subplot 1: Under 40
subplot(2, 2, 1);
scatter(date(N_mhe:N-N_mhe), H_data.u40(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("H-u40"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
hold on 
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on, box on 
ylim([0, max(H_data.u40)*1.2]);
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
scatter(date(N_mhe:N-N_mhe), H_data.mid(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on 
plot(date(N_mhe:N-N_mhe), results.sts.("H-mid"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
hold on
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
ylim([0, max(H_data.mid)*1.2]);
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
scatter(date(N_mhe:N-N_mhe), H_data.old(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("H-old"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
hold on
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on, box on 
ylim([0, max(H_data.old)*1.2]);
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
scatter(date(N_mhe:N-N_mhe), H_data.ger(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("H-ger"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on, box on 
ylim([0, max(H_data.ger)*1.2]);
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('Real Data', 'Model', 'Interpreter', 'latex', 'location', 'northeast');
sgtitle('\textbf{H - Healed}', 'Interpreter', 'latex');


%                                        ---- Plot with variants ----                  %

figure(52)
% Subplot 1: Under 40
subplot(2, 2, 1);
scatter(date(N_mhe:N-N_mhe), H_data.u40(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5,'HandleVisibility', 'off');
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("H-u40"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0], 'HandleVisibility', 'off');
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
ylabel('Variant percentage $\%$', 'Interpreter', 'latex');
ylim([0, 100]);
yyaxis left
ylabel('H-data u40', 'Interpreter', 'latex');
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2, 'HandleVisibility', 'off');
end
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('SARS CoV-2', 'Alpha var.', 'Gamma var.', 'Delta var.', 'Interpreter', 'latex', 'Location', 'northeast');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
scatter(date(N_mhe:N-N_mhe), H_data.mid(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("H-mid"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('H-data mid', 'Interpreter', 'latex');
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
scatter(date(N_mhe:N-N_mhe), H_data.old(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("H-old"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('H-data old','Interpreter', 'latex');  % Label for the primary y-axis
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
scatter(date(N_mhe:N-N_mhe), H_data.ger(N_mhe:N-N_mhe), 50, 'filled', 'MarkerFaceAlpha', 0.5);
hold on
plot(date(N_mhe:N-N_mhe), results.sts.("H-ger"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('H-data ger','Interpreter', 'latex');  % Label for the primary y-axis
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
sgtitle('\textbf{H - Healed}', 'Interpreter', 'latex');

%%  ALPHA Parameter - All plots

%                                        ---- Plot with NPIs ----                  %

figure(14)
% Subplot 1: Under 40
subplot(2, 2, 1);
plot(date(N_mhe:N-N_mhe), results.par.("alp-u40"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
hold on 
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on, box on 
ylim([min(results.par.("alp-u40"))*0.8, max(results.par.("alp-u40"))*1.2]);
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
plot(date(N_mhe:N-N_mhe), results.par.("alp-mid"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
hold on
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
ylim([min(results.par.("alp-mid"))*0.8, max(results.par.("alp-mid"))*1.2]);
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
plot(date(N_mhe:N-N_mhe), results.par.("alp-old"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
hold on
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on, box on 
ylim([min(results.par.("alp-old"))*0.8, max(results.par.("alp-old"))*1.2]);
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
plot(date(N_mhe:N-N_mhe), results.par.("alp-ger"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
hold on 
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on, box on 
ylim([min(results.par.("alp-ger"))*0.8, max(results.par.("alp-ger"))*1.05]);
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('Real Data', 'Model', 'Interpreter', 'latex', 'location', 'southeast');
sgtitle('\textbf{$\alpha$ - Transmission rate}', 'Interpreter', 'latex');


%                                        ---- Plot with variants ----                  %

figure(15)
% Subplot 1: Under 40
subplot(2, 2, 1);
plot(date(N_mhe:N-N_mhe), results.par.("alp-u40"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0], 'HandleVisibility', 'off');
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
ylabel('Variant percentage $\%$', 'Interpreter', 'latex');
ylim([0, 100]);
yyaxis left
ylabel('E-data u40', 'Interpreter', 'latex');
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2, 'HandleVisibility', 'off');
end
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('SARS CoV-2', 'Alpha var.', 'Gamma var.', 'Delta var.', 'Interpreter', 'latex', 'Location', 'northeast');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
plot(date(N_mhe:N-N_mhe), results.par.("alp-mid"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('E-data mid', 'Interpreter', 'latex');
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
plot(date(N_mhe:N-N_mhe), results.par.("alp-old"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('E-data old','Interpreter', 'latex');  % Label for the primary y-axis
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
plot(date(N_mhe:N-N_mhe), results.par.("alp-ger"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('E-data ger','Interpreter', 'latex');  % Label for the primary y-axis
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
sgtitle('\textbf{$\alpha$ - Transmission Rate}', 'Interpreter', 'latex');


%%  GAMMA Parameter - All plots

%                                        ---- Plot with NPIs ----                  %

figure(16)
% Subplot 1: Under 40
subplot(2, 2, 1);
plot(date(N_mhe:N-N_mhe), results.par.("gam-u40"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
hold on 
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on, box on 
ylim([min(results.par.("gam-u40"))*0.8, max(results.par.("gam-u40"))*1.2]);
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
plot(date(N_mhe:N-N_mhe), results.par.("gam-mid"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
hold on
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
ylim([min(results.par.("gam-mid"))*0.8, max(results.par.("gam-mid"))*1.2]);
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
plot(date(N_mhe:N-N_mhe), results.par.("gam-old"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
hold on
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on, box on 
ylim([min(results.par.("gam-old"))*0.8, max(results.par.("gam-old"))*1.2]);
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
plot(date(N_mhe:N-N_mhe), results.par.("gam-ger"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
hold on 
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on, box on 
ylim([min(results.par.("gam-ger"))*0.8, max(results.par.("gam-ger"))*1.05]);
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('Real Data', 'Model', 'Interpreter', 'latex', 'location', 'southeast');
sgtitle('\textbf{$\gamma$ - Detection rate}', 'Interpreter', 'latex');


%                                        ---- Plot with variants ----                  %

figure(17)
% Subplot 1: Under 40
subplot(2, 2, 1);
plot(date(N_mhe:N-N_mhe), results.par.("gam-u40"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0], 'HandleVisibility', 'off');
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
ylabel('Variant percentage $\%$', 'Interpreter', 'latex');
ylim([0, 100]);
yyaxis left
ylabel('E-data u40', 'Interpreter', 'latex');
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2, 'HandleVisibility', 'off');
end
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('SARS CoV-2', 'Alpha var.', 'Gamma var.', 'Delta var.', 'Interpreter', 'latex', 'Location', 'northeast');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
plot(date(N_mhe:N-N_mhe), results.par.("gam-mid"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('E-data mid', 'Interpreter', 'latex');
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
plot(date(N_mhe:N-N_mhe), results.par.("gam-old"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('E-data old','Interpreter', 'latex');  % Label for the primary y-axis
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
plot(date(N_mhe:N-N_mhe), results.par.("gam-ger"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('E-data ger','Interpreter', 'latex');  % Label for the primary y-axis
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
sgtitle('\textbf{$\gamma$ - Detection Rate}', 'Interpreter', 'latex');

%%  DELTA Parameter - All plots

%                                        ---- Plot with NPIs ----                  %

figure(18)
% Subplot 1: Under 40
subplot(2, 2, 1);
plot(date(N_mhe:N-N_mhe), results.par.("dlt-u40"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
hold on 
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on, box on 
ylim([min(results.par.("dlt-u40"))*0.8, max(results.par.("dlt-u40"))*1.2]);
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
plot(date(N_mhe:N-N_mhe), results.par.("dlt-mid"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
hold on
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
ylim([min(results.par.("dlt-mid"))*0.8, max(results.par.("dlt-mid"))*1.2]);
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
plot(date(N_mhe:N-N_mhe), results.par.("dlt-old"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
hold on
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on, box on 
ylim([min(results.par.("dlt-old"))*0.8, max(results.par.("dlt-old"))*1.2]);
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
plot(date(N_mhe:N-N_mhe), results.par.("dlt-ger"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
hold on 
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on, box on 
ylim([min(results.par.("dlt-ger"))*0.8, max(results.par.("dlt-ger"))*1.05]);
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('Real Data', 'Model', 'Interpreter', 'latex', 'location', 'southeast');
sgtitle('\textbf{$\delta$ - Aggravation rate}', 'Interpreter', 'latex');


%                                        ---- Plot with variants ----                  %

figure(19)
% Subplot 1: Under 40
subplot(2, 2, 1);
plot(date(N_mhe:N-N_mhe), results.par.("dlt-u40"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0], 'HandleVisibility', 'off');
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
ylabel('Variant percentage $\%$', 'Interpreter', 'latex');
ylim([0, 100]);
yyaxis left
ylabel('E-data u40', 'Interpreter', 'latex');
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2, 'HandleVisibility', 'off');
end
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('SARS CoV-2', 'Alpha var.', 'Gamma var.', 'Delta var.', 'Interpreter', 'latex', 'Location', 'northeast');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
plot(date(N_mhe:N-N_mhe), results.par.("dlt-mid"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('E-data mid', 'Interpreter', 'latex');
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
plot(date(N_mhe:N-N_mhe), results.par.("dlt-old"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('E-data old','Interpreter', 'latex');  % Label for the primary y-axis
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
plot(date(N_mhe:N-N_mhe), results.par.("dlt-ger"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('E-data ger','Interpreter', 'latex');  % Label for the primary y-axis
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
sgtitle('\textbf{$\delta$ - Aggravation Rate}', 'Interpreter', 'latex');

%%  SIGMA Parameter - All plots

%                                        ---- Plot with NPIs ----                  %

figure(20)
% Subplot 1: Under 40
subplot(2, 2, 1);
plot(date(N_mhe:N-N_mhe), results.par.("sgm-u40"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
hold on 
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on, box on 
ylim([min(results.par.("sgm-u40"))*0.8, max(results.par.("sgm-u40"))*1.2]);
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
plot(date(N_mhe:N-N_mhe), results.par.("sgm-mid"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
hold on
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
ylim([min(results.par.("sgm-mid"))*0.8, max(results.par.("sgm-mid"))*1.2]);
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
plot(date(N_mhe:N-N_mhe), results.par.("sgm-old"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
hold on
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on, box on 
ylim([min(results.par.("sgm-old"))*0.8, max(results.par.("sgm-old"))*1.2]);
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
plot(date(N_mhe:N-N_mhe), results.par.("sgm-ger"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
hold on 
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on, box on 
ylim([min(results.par.("sgm-ger"))*0.8, max(results.par.("sgm-ger"))*1.05]);
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('Real Data', 'Model', 'Interpreter', 'latex', 'location', 'southeast');
sgtitle('\textbf{$\sigma$ - Healing rate}', 'Interpreter', 'latex');


%                                        ---- Plot with variants ----                  %

figure(21)
% Subplot 1: Under 40
subplot(2, 2, 1);
plot(date(N_mhe:N-N_mhe), results.par.("sgm-u40"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0], 'HandleVisibility', 'off');
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
ylabel('Variant percentage $\%$', 'Interpreter', 'latex');
ylim([0, 100]);
yyaxis left
ylabel('E-data u40', 'Interpreter', 'latex');
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2, 'HandleVisibility', 'off');
end
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('SARS CoV-2', 'Alpha var.', 'Gamma var.', 'Delta var.', 'Interpreter', 'latex', 'Location', 'northeast');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
plot(date(N_mhe:N-N_mhe), results.par.("sgm-mid"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('E-data mid', 'Interpreter', 'latex');
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
plot(date(N_mhe:N-N_mhe), results.par.("sgm-old"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('E-data old','Interpreter', 'latex');  % Label for the primary y-axis
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
plot(date(N_mhe:N-N_mhe), results.par.("sgm-ger"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('E-data ger','Interpreter', 'latex');  % Label for the primary y-axis
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
sgtitle('\textbf{$\sigma$ - Healing Rate}', 'Interpreter', 'latex');

%%  TAU Parameter - All plots

%                                        ---- Plot with NPIs ----                  %

figure(22)
% Subplot 1: Under 40
subplot(2, 2, 1);
plot(date(N_mhe:N-N_mhe), results.par.("tau-u40"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
hold on 
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on, box on 
ylim([min(results.par.("tau-u40"))*0.8, max(results.par.("tau-u40"))*1.2]);
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
plot(date(N_mhe:N-N_mhe), results.par.("tau-mid"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
hold on
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
ylim([min(results.par.("tau-mid"))*0.8, max(results.par.("tau-mid"))*1.2]);
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
plot(date(N_mhe:N-N_mhe), results.par.("tau-old"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
hold on
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on, box on 
ylim([min(results.par.("tau-old"))*0.8, max(results.par.("tau-old"))*1.2]);
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
plot(date(N_mhe:N-N_mhe), results.par.("tau-ger"), 'LineWidth', 2, 'MarkerSize', 5, 'Color',[0 0 0]);
hold on 
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    fill(areaNPI.x(ii, :) ,areaNPI.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on, box on 
ylim([min(results.par.("tau-ger"))*0.8, max(results.par.("tau-ger"))*1.05]);
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('Real Data', 'Model', 'Interpreter', 'latex', 'location', 'southeast');
sgtitle('\textbf{$\tau$ - Mortality rate}', 'Interpreter', 'latex');


%                                        ---- Plot with variants ----                  %

figure(23)
% Subplot 1: Under 40
subplot(2, 2, 1);
plot(date(N_mhe:N-N_mhe), results.par.("tau-u40"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0], 'HandleVisibility', 'off');
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
ylabel('Variant percentage $\%$', 'Interpreter', 'latex');
ylim([0, 100]);
yyaxis left
ylabel('E-data u40', 'Interpreter', 'latex');
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2, 'HandleVisibility', 'off');
end
title('\textbf{Under 40}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
lgd = legend('SARS CoV-2', 'Alpha var.', 'Gamma var.', 'Delta var.', 'Interpreter', 'latex', 'Location', 'northeast');
% Subplot 2: Middle aged (40-59)
subplot(2, 2, 2);
plot(date(N_mhe:N-N_mhe), results.par.("tau-mid"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('E-data mid', 'Interpreter', 'latex');
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Middle aged (40-59)}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 3: Senior (60-80)
subplot(2, 2, 3);
plot(date(N_mhe:N-N_mhe), results.par.("tau-old"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('E-data old','Interpreter', 'latex');  % Label for the primary y-axis
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Senior (60-80)}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
subplot(2, 2, 4);
plot(date(N_mhe:N-N_mhe), results.par.("tau-ger"), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
yyaxis right
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
ylabel('E-data ger','Interpreter', 'latex');  % Label for the primary y-axis
grid on;
hold off;
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
end
title('\textbf{Geriatric (80+)}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
sgtitle('\textbf{$\tau$ - Mortality Rate}', 'Interpreter', 'latex');