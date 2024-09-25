%% Plotting BayesOpti MHE og Results
clc
clear all
close all
addpath('/Users/marcodelloro/Downloads/casadi-3.6.3-osx64-matlab2018b')
import casadi.*
load('SIDTTHE_data_DEF.mat');
set(0,'DefaultFigureWindowStyle','docked')
load('BayesResult_ThursNight.mat')
load("WeightsOpti-ResLambdaUppato.mat")


Npop = 59240329; % Total Population of Italy
N = 399;         % MPC horizon length
N_mhe = 21;      % Estimation horizon (3 weeks)
date = SIDTTHE_data{1,1}.date;
Ts = 1;

I_data = SIDTTHE_data{1,1}.data / Npop;     
D_data = SIDTTHE_data{2,1}.data / Npop;     
T_data = (SIDTTHE_data{3,1}.data + SIDTTHE_data{4,1}.data) / Npop;
H_data = SIDTTHE_data{7,1}.data / Npop;
H_dataAug = SIDTTHE_data{6,1}.data / Npop;
E_data = SIDTTHE_data{5,1}.data  / Npop;
S_data = ones(length(I_data),1)' - (I_data + D_data + T_data + H_dataAug + E_data );

ymeas = [S_data; I_data; D_data; T_data; H_dataAug; E_data]; % creation of the measurement vector
n_meas = size(ymeas,1);

%% Plotting section

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
    area.x(ii, :) = [policy_dates(ii) policy_dates(ii) policy_dates(ii+1) policy_dates(ii+1)];
    area.y_alpha(ii, :) = [0 max(results.par.alp)*1.5 max(results.par.alp)*1.5 0];
end

% States
figure(1)
scatter(date(N_mhe:N-N_mhe), S_data(N_mhe:N-N_mhe),20,'filled')
hold on 
plot(date(N_mhe:N-N_mhe),results.sts.S, 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)])
title('\textbf{S}','Interpreter','latex')
grid on
set(gca, 'TickLabelInterpreter', 'Latex')
lgd = legend('Real Data', 'MHE Fitted Data','Interpreter','latex','location','northeast');
lgd.FontSize = 18;

figure(2)
sct = scatter(date(N_mhe:N-N_mhe), I_data(N_mhe:N-N_mhe),100,'filled','MarkerEdgeAlpha', 0.3,'MarkerFaceAlpha', 0.3); % Set face transparency
colorsct=[0, 0.4470, 0.7410];
hold on 
plot(date(N_mhe:N-N_mhe),results.sts.I, 'LineWidth', 2, 'MarkerSize', 5,'Color',[0.8500 0.3250 0.0980]);
xlim([date(1+N_mhe), date(end-N_mhe)])
% title('\textbf{\textit{I} - Infected individuals}','Interpreter','latex')
yax = ylabel('Normalized Population','Interpreter','latex');
grid on
set(gca, 'TickLabelInterpreter', 'Latex')
lgd = legend('Estimated Data', '95\% confidence interval','MHE Fitted Data','Interpreter','latex','location','northeast');
box on 

figure(3)
scatter(date(N_mhe:N-N_mhe), D_data(N_mhe:N-N_mhe),100,'filled','MarkerEdgeAlpha', 0.8,'MarkerFaceAlpha', 0.8)
hold on 
plot(date(N_mhe:N-N_mhe),results.sts.D, 'LineWidth', 2, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)])
% title('\textbf{\textit{D} - Detected individuals}','Interpreter','latex')
grid on
set(gca, 'TickLabelInterpreter', 'Latex')
yax = ylabel('Normalized Population','Interpreter','latex');
lgd = legend('Real Data', 'MHE Fitted Data','Interpreter','latex','location','northeast');
box on 

figure(4)
scatter(date(N_mhe:N-N_mhe), T_data(N_mhe:N-N_mhe),100,'filled','MarkerEdgeAlpha', 0.8,'MarkerFaceAlpha', 0.8)
hold on 
plot(date(N_mhe:N-N_mhe),results.sts.T, 'LineWidth', 2, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)])
% title('\textbf{\textit{T} - Hospitalised and ICUs individuals}','Interpreter','latex')
grid on
set(gca, 'TickLabelInterpreter', 'Latex')
yax = ylabel('Normalized Population','Interpreter','latex');
box on 

figure(5)
scatter(date(N_mhe:N-N_mhe), H_data(N_mhe:N-N_mhe),100,'filled','MarkerEdgeAlpha', 0.8,'MarkerFaceAlpha', 0.8)
hold on 
% scatter(date(N_mhe:N-N_mhe), H_dataAug(N_mhe:N-N_mhe),20,'filled')
% hold on 
plot(date(N_mhe:N-N_mhe),results.sts.H, 'LineWidth', 2, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)])
grid on
set(gca, 'TickLabelInterpreter', 'Latex')
yax = ylabel('Normalized Population','Interpreter','latex');
box on 

figure(6)
scatter(date(N_mhe:N-N_mhe), E_data(N_mhe:N-N_mhe),100,'filled','MarkerEdgeAlpha', 0.8,'MarkerFaceAlpha', 0.8)
hold on 
plot(date(N_mhe:N-N_mhe),results.sts.E, 'LineWidth', 2, 'MarkerSize', 5);
hold on 
xlim([date(1+N_mhe), date(end-N_mhe)])
% title('\textbf{\textit{D} - Expired individuals}','Interpreter','latex')
grid on
set(gca, 'TickLabelInterpreter', 'Latex')
yax = ylabel('Normalized Population','Interpreter','latex');
box on 

figure(10)
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    
    fill(area.x(ii, :) ,area.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
hold on
plot(date(N_mhe:N-N_mhe), results.par.alp,'k','LineWidth',2, 'DisplayName', '$\alpha$','HandleVisibility', 'off')
yax = ylabel('Coefficient Value','Interpreter','latex');
% title('\textbf{$\alpha$ coefficient}','Interpreter','latex')
grid on
lgd = legend('Mild Restrictions', 'NPIs \& Social Limitations','Curfew','Total Lockdown','Interpreter','latex','location','northeast');
title(lgd, '\textbf{Policy Level}');
xlim([date(1+N_mhe), date(end-N_mhe)])
ylim([0, max(results.par.alp)*1.5])
set(gca, 'TickLabelInterpreter', 'Latex')

figure(11)
plot(date(N_mhe:N-N_mhe),results.par.gam, 'LineWidth', 2, 'Color','k')
xlim([date(1+N_mhe), date(end-N_mhe)])
ylim([0, max(results.par.gam)*1.5])
grid on
set(gca, 'TickLabelInterpreter', 'Latex')
yax = ylabel('Coefficient Value','Interpreter','latex');
title('\textbf{$\gamma$}','Interpreter','latex')

figure(12)
plot(date(N_mhe:N-N_mhe),results.par.dlt, 'LineWidth', 2, 'Color','k')
xlim([date(1+N_mhe), date(end-N_mhe)])
ylim([0, max(results.par.dlt)*1.2])
grid on
set(gca, 'TickLabelInterpreter', 'Latex')
yax = ylabel('Coefficient Value','Interpreter','latex');
title('\textbf{$\delta$}','Interpreter','latex')

figure(13)
plot(date(N_mhe:N-N_mhe),results.par.sgm, 'LineWidth', 2, 'Color','k')
xlim([date(1+N_mhe), date(end-N_mhe)])
ylim([0, max(results.par.sgm)*1.5])
grid on
set(gca, 'TickLabelInterpreter', 'Latex')
yax = ylabel('Coefficient Value','Interpreter','latex');
title('\textbf{$\sigma$}','Interpreter','latex')

figure(14)
plot(date(N_mhe:N-N_mhe),results.par.tau, 'LineWidth', 2, 'Color','k')
xlim([date(1+N_mhe), date(end-N_mhe)])
ylim([0, max(results.par.tau)*1.5])
grid on
set(gca, 'TickLabelInterpreter', 'Latex')
yax = ylabel('Coefficient Value','Interpreter','latex');
title('\textbf{$\tau$}','Interpreter','latex')


figure(90)
scatter(date(N_mhe:N-N_mhe), 50*T_data(N_mhe:N-N_mhe),100,'filled','MarkerEdgeAlpha', 0.8,'MarkerFaceAlpha', 0.8)
hold on 
plot(date(N_mhe:N-N_mhe),results.par.dlt, 'LineWidth', 2, 'MarkerSize', 5);
hold on 
xlim([date(1+N_mhe), date(end-N_mhe)])
hold on
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    
    fill(area.x(ii, :) ,0.1*area.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
% title('\textbf{\textit{D} - Expired individuals}','Interpreter','latex')
grid on
set(gca, 'TickLabelInterpreter', 'Latex')
yax = ylabel('Normalized Population','Interpreter','latex');
box on 