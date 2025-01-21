clc
clear 
close all
load '/Users/marcodelloro/Desktop/Pandemic-System-Identification/MHE Age Stratified/SmoothVariantasIta.mat'
load '/Users/marcodelloro/Desktop/Pandemic-System-Identification/MHE Age Stratified - Single Alpha'/Domenica.mat
addpath('/Users/marcodelloro/Desktop/Pandemic-System-Identification/Reconstructed Datasets/')

%% Data Loading

S_data = readtable('Reconstructed_S.csv'); I_raw = readtable('Infected_I.csv'); D_data = readtable('Reconstructed_D.csv');     
T1_data = readtable('Reconstructed_T1.csv'); T2_data = readtable('Reconstructed_T2.csv'); H_data = readtable('Reconstructed_H.csv');
E_data = readtable('Reconstructed_E.csv');

N_u40 = 23205875; N_mid = 18142711;  N_old = 13408810; N_ger = 3951057;  
N_mhe = 21;  N = 399; Npop = N_u40 + N_mid + N_old + N_ger; 
date = datetime('31-Aug-2020', 'InputFormat', 'dd-MMM-yyyy'):datetime('03-Oct-2021', 'InputFormat', 'dd-MMM-yyyy');

I_data.u40 = [I_raw.u40(1) I_raw.u40'];         I_data.mid = [I_raw.mid(1) I_raw.mid'];
I_data.old = [I_raw.old(1) I_raw.old'];         I_data.ger = [I_raw.ger(1) I_raw.ger'];

T_data.u40 = (T1_data.u40  + T2_data.u40);       T_data.mid = (T1_data.mid  + T2_data.mid);     
T_data.old = (T1_data.old  + T2_data.u40);       T_data.ger = (T1_data.ger  + T2_data.ger);  

%% Plotting Settings

policy_idx = [1 40 69 116 141 190 242 294 399];

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
    areaNPI.y_alpha(ii, :) = [0 max(results.par.("alp-u40"))*1.5*Npop max(results.par.("alp-u40"))*1.5*Npop 0];
end

colors = {[0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250], [ 0.4940, 0.1840, 0.5560], [0.4660, 0.6740, 0.1880]};

specDates = [date(21) date(141) date(261)];
scatsize = 20;
linesize = 2;


%% MHE example plots

tinit = 200;
scatsize = 70;
ahead = 10;
resample_interval = 7;
resampled_indices = tinit:resample_interval:tinit+2*N_mhe;
resampled_indices2 = tinit:resample_interval:tinit+2*N_mhe+7;
% first plot box
last_index = resampled_indices(end);
second_last_index = resampled_indices(end-3); 
x_rect = [date(second_last_index), date(last_index), date(last_index), date(second_last_index)];
y_rect = [1850, 1850, 2500, 2500]; 
% second plot box 
last_index2 = resampled_indices2(end);
second_last_index2 = resampled_indices2(end-3); 
x_rect2 = [date(second_last_index2), date(last_index2), date(last_index2), date(second_last_index2)];
y_rect2 = [1850, 1850, 2500, 2500]; 

figure(100)
t = tiledlayout(2,1, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile
fill(x_rect, y_rect, [1, 0.5, 0.0], 'EdgeColor',[1, 0.5, 0.0] , 'FaceAlpha',0.3 );
hold on
scatter(date(resampled_indices), T_data.u40(resampled_indices), scatsize, 'filled', 'MarkerFaceColor', [0, 0.4470, 0.7410],'MarkerEdgeColor', [0, 0.4470, 0.7410]);
hold on;
plot(date(resampled_indices), results.sts.("T-u40")(resampled_indices-N_mhe)*Npop, '-gs',Color='red',MarkerFaceColor='red',MarkerSize=8,LineWidth=1.5)
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
xlim([date(resampled_indices2(1)) date(last_index2)])
ylim([1600 2600])
title('Time instant $t=k$','Interpreter', 'latex');
set(gca, 'XTickLabel', [], 'YTickLabel', []);
nexttile 
fill(x_rect, y_rect, [1, 1, 1] ,'EdgeColor', [1, 0.5, 0.0], 'LineStyle', '--','FaceAlpha',0 );  
hold on
fill(x_rect2, y_rect2, [1, 0.5, 0.0], 'EdgeColor',[1, 0.5, 0.0] , 'FaceAlpha',0.3 );
hold on
scatter(date(resampled_indices2), T_data.u40(resampled_indices2), scatsize, 'filled','MarkerFaceColor', [0, 0.4470, 0.7410],'MarkerEdgeColor', [0, 0.4470, 0.7410]);
hold on;
plot(date(resampled_indices2), results.sts.("T-u40")(resampled_indices2-N_mhe)*Npop, '-gs',Color='red',MarkerFaceColor='red',MarkerSize=8,LineWidth=1.5)
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
xlim([date(resampled_indices2(1)) date(last_index2)])
ylim([1600 2600])
title('Time instant $t=k+1$','Interpreter', 'latex')
set(gca, 'XTickLabel', [], 'YTickLabel', []);
annotation('rectangle', [0 0 1 1], 'LineWidth', 2, 'Color', 'black');   


%% D Compartment ALL plots

figure(1)
t = tiledlayout(2,2, 'Padding', 'compact', 'TileSpacing', 'compact');
ylabel(t, 'Detected Individuals', 'Interpreter', 'latex','FontSize', 20);
nexttile
% Subplot 1: Under 40
scatter(date(N_mhe:N-N_mhe), D_data.u40(N_mhe:N-N_mhe), scatsize, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("D-u40")*Npop, 'LineWidth', linesize, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,max(results.sts.("D-u40"))*1.05*Npop])
title('\textbf{$<$ 40}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
xticks(specDates); 
xticklabels(datestr(specDates, 'mm-yy')); 
% Subplot 2: Middle aged (40-59)
nexttile
scatter(date(N_mhe:N-N_mhe), D_data.mid(N_mhe:N-N_mhe), scatsize, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on 
plot(date(N_mhe:N-N_mhe), results.sts.("D-mid")*Npop, 'LineWidth', linesize, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,max(results.sts.("D-mid"))*1.05*Npop])
title('\textbf{40-59}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
xticks(specDates); 
xticklabels(datestr(specDates, 'mm-yy')); 
% Subplot 3: Senior (60-80)
nexttile
scatter(date(N_mhe:N-N_mhe), D_data.old(N_mhe:N-N_mhe), scatsize, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("D-old")*Npop, 'LineWidth', linesize, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,max(results.sts.("D-old"))*1.05*Npop])
title('\textbf{60-79}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
xticks(specDates); 
xticklabels(datestr(specDates, 'mm-yy')); 
% Subplot 4: Geriatric (80+)
nexttile
scatter(date(N_mhe:N-N_mhe), D_data.ger(N_mhe:N-N_mhe), scatsize, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("D-ger")*Npop, 'LineWidth', linesize, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,max(results.sts.("D-ger"))*1.05*Npop])
title('\textbf{$\geq$ 80}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
xticks(specDates); 
xticklabels(datestr(specDates, 'mm-yy')); 
% lgd = legend('Real Data', 'Model', 'Interpreter', 'latex', 'location', 'northeast');
% sgtitle('\textbf{D - Detected}', 'Interpreter', 'latex');


%% T Compartment ALL plots

figure(2)
t = tiledlayout(2,2, 'Padding', 'compact', 'TileSpacing', 'compact');
ylabel(t, 'Threatened Individuals', 'Interpreter', 'latex','FontSize', 20);
nexttile
% Subplot 1: Under 40
scatter(date(N_mhe:N-N_mhe), T_data.u40(N_mhe:N-N_mhe), scatsize, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("T-u40")*Npop, 'LineWidth', linesize, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,max(results.sts.("T-u40"))*1.05*Npop])
title('\textbf{$<$ 40}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
xticks(specDates); 
xticklabels(datestr(specDates, 'mm-yy'));
ax = gca;
ax.YRuler.Exponent = 3;
% Subplot 2: Middle aged (40-59)
nexttile
scatter(date(N_mhe:N-N_mhe), T_data.mid(N_mhe:N-N_mhe), scatsize, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on 
plot(date(N_mhe:N-N_mhe), results.sts.("T-mid")*Npop, 'LineWidth', linesize, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,max(results.sts.("T-mid"))*1.05*Npop])
title('\textbf{40-59}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
xticks(specDates); 
xticklabels(datestr(specDates, 'mm-yy')); 
ax = gca;
ax.YRuler.Exponent = 3; 
% Subplot 3: Senior (60-80)
nexttile
scatter(date(N_mhe:N-N_mhe), T_data.old(N_mhe:N-N_mhe), scatsize, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("T-old")*Npop, 'LineWidth', linesize, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,max(results.sts.("T-old"))*1.05*Npop])
title('\textbf{60-79}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
xticks(specDates); 
xticklabels(datestr(specDates, 'mm-yy')); 
ax = gca;
ax.YRuler.Exponent = 3;
% Subplot 4: Geriatric (80+)
nexttile
scatter(date(N_mhe:N-N_mhe), T_data.ger(N_mhe:N-N_mhe), scatsize, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("T-ger")*Npop, 'LineWidth', linesize, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,max(results.sts.("T-ger"))*1.05*Npop])
title('\textbf{$\geq$ 80}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
xticks(specDates); 
xticklabels(datestr(specDates, 'mm-yy')); 
ax = gca;
ax.YRuler.Exponent = 3;
lgd = legend('Data', 'Model', 'Interpreter', 'latex', 'location', 'northeast');
% sgtitle('\textbf{D - Detected}', 'Interpreter', 'latex');

%% I Compartment ALL plots 

figure(3)
t = tiledlayout(2,2, 'Padding', 'compact', 'TileSpacing', 'compact');
ylabel(t, 'Infected Individuals', 'Interpreter', 'latex','FontSize', 20);
nexttile
scatter(date(N_mhe:N-N_mhe), I_data.u40(N_mhe:N-N_mhe), scatsize, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("I-u40")*Npop, 'LineWidth', linesize);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,max(results.sts.("I-u40"))*1.05*Npop])
title('\textbf{$<$ 40}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
xticks(specDates); 
xticklabels(datestr(specDates, 'mm-yy')); 
ax = gca;
ax.YRuler.Exponent = 3;
% Subplot 2: Middle aged (40-59)
nexttile
scatter(date(N_mhe:N-N_mhe), I_data.mid(N_mhe:N-N_mhe), scatsize, 'filled', 'MarkerFaceAlpha', 0.5);
hold on 
plot(date(N_mhe:N-N_mhe), results.sts.("I-mid")*Npop, 'LineWidth', linesize);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,max(results.sts.("I-mid"))*1.05*Npop])
title('\textbf{40-59}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
xticks(specDates); 
xticklabels(datestr(specDates, 'mm-yy')); 
ax = gca;
ax.YRuler.Exponent = 3;
% Subplot 3: Senior (60-80)
nexttile
scatter(date(N_mhe:N-N_mhe), I_data.old(N_mhe:N-N_mhe), scatsize, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("I-old")*Npop, 'LineWidth', linesize);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{60-79}', 'Interpreter', 'latex');
grid on, box on 
ylim([0,max(results.sts.("I-old"))*1.05*Npop])
set(gca, 'TickLabelInterpreter', 'Latex');
xticks(specDates); 
xticklabels(datestr(specDates, 'mm-yy')); 
ax = gca;
ax.YRuler.Exponent = 3;
% Subplot 4: Geriatric (80+)
nexttile
scatter(date(N_mhe:N-N_mhe), I_data.ger(N_mhe:N-N_mhe), scatsize, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("I-ger")*Npop, 'LineWidth', linesize);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0,max(results.sts.("I-ger"))*1.05*Npop])
title('\textbf{$\geq$ 80}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
xticks(specDates); 
xticklabels(datestr(specDates, 'mm-yy')); 
ax = gca;
ax.YRuler.Exponent = 3;
% lgd = legend('Estimated Data', 'Model Estimation', 'Interpreter', 'latex', 'location', 'southoutside', 'orientation', 'horizontal');
% sgtitle('\textbf{I - Infected}', 'Interpreter', 'latex');

%% E Compartment ALL plots

figure(4)
t = tiledlayout(2,2, 'Padding', 'compact', 'TileSpacing', 'compact');
ylabel(t, 'Deceased Individuals', 'Interpreter', 'latex','FontSize', 20);
nexttile
scatter(date(N_mhe:N-N_mhe), E_data.u40(N_mhe:N-N_mhe), scatsize, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("E-u40")*Npop, 'LineWidth', linesize);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{$<$ 40}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
xticks(specDates); 
xticklabels(datestr(specDates, 'mm-yy')); 
ax = gca;
ax.YRuler.Exponent = 3;
% Subplot 2: Middle aged (40-59)
nexttile
scatter(date(N_mhe:N-N_mhe), E_data.mid(N_mhe:N-N_mhe), scatsize, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on 
plot(date(N_mhe:N-N_mhe), results.sts.("E-mid")*Npop, 'LineWidth', linesize);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{40-59}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
xticks(specDates); 
xticklabels(datestr(specDates, 'mm-yy')); 
ax = gca;
ax.YRuler.Exponent = 3;
% Subplot 3: Senior (60-80)
nexttile
scatter(date(N_mhe:N-N_mhe), E_data.old(N_mhe:N-N_mhe), scatsize, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("E-old")*Npop, 'LineWidth', linesize);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{60-79}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
xticks(specDates); 
xticklabels(datestr(specDates, 'mm-yy')); 
ax = gca;
ax.YRuler.Exponent = 3;
% Subplot 4: Geriatric (80+)
nexttile
scatter(date(N_mhe:N-N_mhe), E_data.ger(N_mhe:N-N_mhe), scatsize, 'filled', 'MarkerFaceAlpha', 0.5); % Increased size and added transparency
hold on;
plot(date(N_mhe:N-N_mhe), results.sts.("E-ger")*Npop, 'LineWidth', linesize);
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{$\geq$ 80}', 'Interpreter', 'latex');
grid on, box on 
set(gca, 'TickLabelInterpreter', 'Latex');
xticks(specDates); 
xticklabels(datestr(specDates, 'mm-yy')); 
ax = gca;
ax.YRuler.Exponent = 3;
% lgd = legend('Real Data', 'Model', 'Interpreter', 'latex', 'location', 'northeast');
% sgtitle('\textbf{E - Deceased}', 'Interpreter', 'latex');

%% alpha plots

figure(5);
yyaxis right;
set(gca, 'YColor', [0 0 0]);
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
%ylabel('Variant percentage $\%$', 'Interpreter', 'latex');
ylim([0, 100]);
yyaxis left;
set(gca, 'YColor', [0 0 0]);
plot3(date(N_mhe:N-N_mhe), results.par.("alp-u40"), 100000*ones(size(results.par.("alp-u40"))), 'LineWidth', 2, 'MarkerSize', 5, 'Color', [0 0 0], 'HandleVisibility', 'off');
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylabel('$\alpha$ parameter', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
xticks(specDates); 
xticklabels(datestr(specDates, 'mm-yy')); 
lgd = legend('Wild Type', 'Alpha var.', 'Gamma var.', 'Delta var.', 'Interpreter', 'latex', 'Location', 'southoutside', 'Orientation', 'horizontal');
% lgd.NumColumns = 2;
lgd.NumColumns = 4;
ylim([0, max(results.par.("alp-u40"))*1.3]);
hold off;


%% Delta plots 

figure(6)
% Subplot 1: Under 40
t = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
ylabel(t, '$\delta$ parameter', 'Interpreter', 'latex','FontSize', 20)

nexttile 
plot(date(N_mhe:N-N_mhe), results.par.("dlt-u40"), 'LineWidth', linesize, 'MarkerSize', 5, 'Color', [0 0 0], 'HandleVisibility', 'off');
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0, max(results.par.("dlt-u40"))*1.3])
ax = gca;
ax.YRuler.Exponent = -3;
xticks(specDates); 
xticklabels(datestr(specDates, 'mm-yy'));
set(gca, 'TickLabelInterpreter', 'latex')
yyaxis right
set(gca, 'YColor', [0 0 0]);
set(gca, 'YTickLabel', []); 
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
ylim([0, 100]);
yyaxis left
grid on;
hold off;
% for idx = policy_idx
%     xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2, 'HandleVisibility', 'off');
% end
title('\textbf{ $<$ 40}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 2: Middle aged (40-59)
nexttile 
plot(date(N_mhe:N-N_mhe), results.par.("dlt-mid"), 'LineWidth', linesize, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0, max(results.par.("dlt-mid"))*1.3])
xticks(specDates); 
xticklabels(datestr(specDates, 'mm-yy'));
set(gca, 'TickLabelInterpreter', 'latex')
ax = gca;
ax.YRuler.Exponent = -3;
yyaxis right
set(gca, 'YColor', [0 0 0]);
set(gca, 'YTickLabel', []); 
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
grid on;
hold off;
% for idx = policy_idx
%     xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
% end
title('\textbf{40-59}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 3: Senior (60-80)
nexttile 
plot(date(N_mhe:N-N_mhe), results.par.("dlt-old"), 'LineWidth', linesize, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0, max(results.par.("dlt-old"))*1.3])
ax = gca;
ax.YRuler.Exponent = -3;
yyaxis right
xticks(specDates); 
xticklabels(datestr(specDates, 'mm-yy'));
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
grid on;
hold off;
% for idx = policy_idx
%     xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
% end
title('\textbf{60-79}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
nexttile 
plot(date(N_mhe:N-N_mhe), results.par.("dlt-ger"), 'LineWidth', linesize, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0, max(results.par.("dlt-ger"))*1.3])
ax = gca;
ax.YRuler.Exponent = -3;
yyaxis right
xticks(specDates); 
xticklabels(datestr(specDates, 'mm-yy'));
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
grid on;
hold off;
% for idx = policy_idx
%     xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
% end
title('\textbf{$\geq$ 80}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');


%% Tau plots 

figure(7)
% Subplot 1: Under 40
t = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
ylabel(t,'$\tau$ parameter', 'Interpreter', 'latex','FontSize', 20);
nexttile 
plot(date(N_mhe:N-N_mhe), results.par.("tau-u40"), 'LineWidth', linesize, 'MarkerSize', 5, 'Color', [0 0 0], 'HandleVisibility', 'off');
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0, max(results.par.("tau-u40"))*1.3])
ax = gca;
ax.YRuler.Exponent = -3;
xticks(specDates); 
xticklabels(datestr(specDates, 'mm-yy'));
set(gca, 'TickLabelInterpreter', 'latex')
yyaxis right
set(gca, 'YColor', [0 0 0]);
set(gca, 'YTickLabel', []); 
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
ylim([0, 100]);
yyaxis left
grid on;
hold off;
% for idx = policy_idx
%     xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2, 'HandleVisibility', 'off');
% end
title('\textbf{ $<$ 40}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 2: Middle aged (40-59)
nexttile 
plot(date(N_mhe:N-N_mhe), results.par.("tau-mid"), 'LineWidth', linesize, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0, max(results.par.("tau-mid"))*1.3])
xticks(specDates); 
xticklabels(datestr(specDates, 'mm-yy'));
set(gca, 'TickLabelInterpreter', 'latex')
ax = gca;
ax.YRuler.Exponent = -3;
yyaxis right
set(gca, 'YColor', [0 0 0]);
set(gca, 'YTickLabel', []); 
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
grid on;
hold off;
% for idx = policy_idx
%     xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
% end
title('\textbf{40-59}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 3: Senior (60-80)
nexttile 
plot(date(N_mhe:N-N_mhe), results.par.("tau-old"), 'LineWidth', linesize, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0, max(results.par.("tau-old"))*1.3])
ax = gca;
ax.YRuler.Exponent = -3;
yyaxis right
xticks(specDates); 
xticklabels(datestr(specDates, 'mm-yy'));
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
grid on;
hold off;
% for idx = policy_idx
%     xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
% end
title('\textbf{60-79}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
nexttile 
plot(date(N_mhe:N-N_mhe), results.par.("tau-ger"), 'LineWidth', linesize, 'MarkerSize', 5, 'Color', [0 0 0]);
xlim([date(1+N_mhe), date(end-N_mhe)]);
ylim([0, max(results.par.("tau-ger"))*1.3])
ax = gca;
ax.YRuler.Exponent = -3;
yyaxis right
xticks(specDates); 
xticklabels(datestr(specDates, 'mm-yy'));
area(SmoothVar.date, SmoothVar.SC2, 'FaceAlpha', 0.3, 'FaceColor', colors{1}, 'EdgeColor', colors{1});
hold on;
area(SmoothVar.date, SmoothVar.alpha, 'FaceAlpha', 0.3, 'FaceColor', colors{2}, 'EdgeColor', colors{2});
hold on;
area(SmoothVar.date, SmoothVar.gamma, 'FaceAlpha', 0.3, 'FaceColor', colors{3}, 'EdgeColor', colors{3});
hold on;
area(SmoothVar.date, SmoothVar.delta, 'FaceAlpha', 0.3, 'FaceColor', colors{4}, 'EdgeColor', colors{4});
yyaxis left
grid on;
hold off;
% for idx = policy_idx
%     xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 2);
% end
title('\textbf{$\geq$ 80}', 'Interpreter', 'latex');
grid on;
box on;
set(gca, 'TickLabelInterpreter', 'Latex');

%% Gamma parameter

figure(8)
t = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
ylabel(t,'$\gamma$ parameter', 'Interpreter', 'latex','FontSize', 20);
nexttile 
plot(date(N_mhe:N-N_mhe), results.par.("gam-u40"), 'LineWidth', linesize, 'MarkerSize', 5, 'Color',[0 0 0]);
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
title('\textbf{ $<$ 40}', 'Interpreter', 'latex');
grid on, box on 
ylim([0, max(results.par.("gam-u40"))*1.5]);
ax = gca;
ax.YRuler.Exponent = -2;
xticks(specDates); 
xticklabels(datestr(specDates, 'mm-yy'));
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 2: Middle aged (40-59)
nexttile
plot(date(N_mhe:N-N_mhe), results.par.("gam-mid"), 'LineWidth', linesize, 'MarkerSize', 5, 'Color',[0 0 0]);
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
title('\textbf{40-59}', 'Interpreter', 'latex');
grid on, box on 
ax = gca;
ax.YRuler.Exponent = -2;
xticks(specDates); 
xticklabels(datestr(specDates, 'mm-yy'));
ylim([0, max(results.par.("gam-mid"))*1.2]);
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 3: Senior (60-80)
nexttile
plot(date(N_mhe:N-N_mhe), results.par.("gam-old"), 'LineWidth', linesize, 'MarkerSize', 5, 'Color',[0 0 0]);
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
title('\textbf{60-79}', 'Interpreter', 'latex');
grid on, box on 
ylim([0, max(results.par.("gam-old"))*1.2]);
ax = gca;
ax.YRuler.Exponent = -2;
xticks(specDates); 
xticklabels(datestr(specDates, 'mm-yy'));
set(gca, 'TickLabelInterpreter', 'Latex');
% Subplot 4: Geriatric (80+)
nexttile
plot(date(N_mhe:N-N_mhe), results.par.("gam-ger"), 'LineWidth', linesize, 'MarkerSize', 5, 'Color',[0 0 0]);
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
hold on
for idx = policy_idx
    xline(date(idx), '--', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1.5);
end
xlim([date(1+N_mhe), date(end-N_mhe)]);
title('\textbf{$\geq$ 80}', 'Interpreter', 'latex');
grid on, box on 
ylim([0, max(results.par.("gam-ger"))*1.5]);
ax = gca;
ax.YRuler.Exponent = -2;
xticks(specDates); 
xticklabels(datestr(specDates, 'mm-yy'));
set(gca, 'TickLabelInterpreter', 'Latex');
% lgd = legend('Real Data', 'Model', 'Interpreter', 'latex', 'location', 'southeast');
% sgtitle('\textbf{$\gamma$ - Detection rate}', 'Interpreter', 'latex');


