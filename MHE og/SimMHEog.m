%% Moving Horizon Estimation
clc
clear all
close all
addpath('/Users/marcodelloro/Downloads/casadi-3.6.3-osx64-matlab2018b')
import casadi.*
load('SIDTTHE_data_DEF.mat');
set(0,'DefaultFigureWindowStyle','docked')

load('BayesResult_ThursNight.mat')

%%  Data Loading and Initialization
tic

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

%% MHE simulation
    
    %                        ----- CasADi states & dynamics initialization -----
    
    s = casadi.SX.sym('s',1,1);   alp = casadi.SX.sym('alp',1,1);   v1 = casadi.SX.sym('v1',1,1);
    i = casadi.SX.sym('i',1,1);   gma = casadi.SX.sym('gma',1,1);   v2 = casadi.SX.sym('v2',1,1);
    d = casadi.SX.sym('d',1,1);   dlt = casadi.SX.sym('dlt',1,1);   v3 = casadi.SX.sym('v3',1,1);
    t = casadi.SX.sym('t',1,1);   sgm = casadi.SX.sym('sgm',1,1);   v4 = casadi.SX.sym('v4',1,1);
    h = casadi.SX.sym('h',1,1);   tau = casadi.SX.sym('tau',1,1);   v5 = casadi.SX.sym('v5',1,1);
    e = casadi.SX.sym('e',1,1);   lambda = 1/10;                    v6 = casadi.SX.sym('v6',1,1);
    
    x = [s; i; d; t; h; e; alp; gma; dlt; sgm; tau];
    n_states = length(x);

   % ------------------------------------------------------------------------------------------------------- %

   %                      ----- States equations posing and computatio of the dynamics constraint ----

   eqn = [    -x(1) * (x(7) * x(2));...
             x(1) * (x(7) * x(2)) - (x(8)+lambda) * x(2);...
             x(2) * x(8) - x(3) * (lambda + x(9));...
             x(9) * x(3) - ( x(11) + x(10) )*x(4);...
             lambda * x(3) + x(4) * x(10) + lambda * x(2);...
             x(11) * x(4);...
             0;...
             0;...
             0;...
             0;...
             0          ];
    
    f = casadi.Function('f',{x},{eqn});

   % ------------------------------------------------------------------------------------------------------- %

   %                                   ----- Optmization Problem -----
    
    opti = casadi.Opti();
    X = opti.variable(n_states, N_mhe);           % Augmented state
    v = opti.variable(6,N_mhe-1);                 % Process noise matrix
    X_tilde = opti.parameter(6,1);                % Arrival Cost on the states (6 states, 4 age groups)
    P_tilde = opti.parameter(5,1);                % Arrival Cost on the params (5 params, 4 age groups)
    Ymeas = opti.parameter(6, N_mhe);             % Available measurement matrix 


    %                  ---- Declarations of the weights of the cost function ----

  
    Z1 = diag([optResults.fullResults.XAtMinObjective.z1_1, optResults.fullResults.XAtMinObjective.z1_2, optResults.fullResults.XAtMinObjective.z1_1,...  
               optResults.fullResults.XAtMinObjective.z1_3, optResults.fullResults.XAtMinObjective.z1_1, optResults.fullResults.XAtMinObjective.z1_1]);
    
    Z2 = diag([optResults.fullResults.XAtMinObjective.z2_1, optResults.fullResults.XAtMinObjective.z2_2, optResults.fullResults.XAtMinObjective.z2_3,...
               optResults.fullResults.XAtMinObjective.z2_4, optResults.fullResults.XAtMinObjective.z2_5]);
    
    Z3 = diag([optResults.fullResults.XAtMinObjective.z3_1, optResults.fullResults.XAtMinObjective.z3_2, optResults.fullResults.XAtMinObjective.z3_3,... 
               optResults.fullResults.XAtMinObjective.z3_3, optResults.fullResults.XAtMinObjective.z3_1, optResults.fullResults.XAtMinObjective.z3_3]);
        

    %                        ---- Computation of the dynamics constraints (RK45) ----

    for jj=1:N_mhe-1 

        k1 = f(X(:,jj));
        k2 = f(X(:,jj)+Ts/2*k1);
        k3 = f(X(:,jj)+Ts/2*k2);
        k4 = f(X(:,jj)+Ts*k3);
        x_plus = Ts/6*(k1+2*k2+2*k3+k4);
        opti.subject_to(X(:,jj+1)== X(:,jj) + x_plus + [v(:,jj); zeros(5,1)]);

    end

    opti.subject_to(X(:) > 0); 
    
    for ii=1:N_mhe
         opti.subject_to( sum(X(1:6,ii)) == 1 )
    end
    
    %                                       ---- Disturbances bounds ----
    
    stdevs = std(ymeas,1,2);
    c = 0.005;
    opti.subject_to(-stdevs(1)*c < v(1,:) < stdevs(1)*c);
    opti.subject_to(-stdevs(2)*c < v(2,:) < stdevs(2)*c);
    opti.subject_to(-stdevs(3)*c < v(3,:) < stdevs(3)*c);
    opti.subject_to(-stdevs(4)*c < v(4,:) < stdevs(4)*c);
    opti.subject_to(-stdevs(5)*c < v(5,:) < stdevs(5)*c);
    opti.subject_to(-stdevs(6)*c < v(6,:) < stdevs(6)*c);
    
  
    %                                       ---- Parameters Bounds ----

    opti.subject_to( 0.05  <= X(7,:)  <= 0.8 ); % alpha constr.
    opti.subject_to( 0.005 <= X(8,:)  <= 0.6 ); % gamma constr.
    opti.subject_to( 1e-4  <= X(9,:)  <= 0.6 ); % delta constr.
    opti.subject_to( 1e-4  <= X(10,:) <= 0.5 ); % sigma constr.
    opti.subject_to( 1e-4  <= X(11,:) <= 0.5 ); % tau constr.
    
    
    %                                        ---- Objective function ----

    obj = sumsqr( Z1*(X(1:6,:) - X_tilde) ) + sumsqr( Z2*(X(7:11,:) - P_tilde) ) + sumsqr( Z3*(Ymeas - X(1:6,:))./max(ymeas, [], 2)  ) + sumsqr(v); 
    opti.minimize(obj);
    
    % Solver Options
    opts.print_time = 0;
    opts.ipopt.print_level = 0;
    opti.solver('ipopt', opts);


   % ------------------------------------------------------------------------------------------------------- %

    %                                     ---- Moving Horizon FOR cycle ----
   
    matrix_sts = []; % matrix to store simulation results (states)
    matrix_par = []; % matrix to store simulation results (states)
 
    for ii = N_mhe:1:N-N_mhe
    
        if ii == N_mhe   
            % First iteration "arrival cost parameters" setting
            opti.set_value(X_tilde, ymeas(1:6,1) );    
            opti.set_value(P_tilde, [0.25; 0.12; 0.01; 0.02; 0.02]);    
        else 
            % Every other iteration "arrival cost parameters"
            opti.set_value(X_tilde, currentState)   
            opti.set_value(P_tilde, currentParam)
        end
    
        % measurement update 
        opti.set_value(Ymeas, ymeas(:,ii-N_mhe+1:ii));
    
        sol = opti.solve();
        
        % Save of the optimization values for STATES and PARAMETERS
        states_dyn = opti.value(X(1:6,end));
        params_dyn= opti.value(X(7:11,end));

        matrix_sts = [matrix_sts, states_dyn];
        matrix_par = [matrix_par,params_dyn];

        % Arrival cost update
        currentState = opti.value(X(1:6,2));
        currentParam = opti.value(X(7:11,1));
    
    end

%                                     ---- Optimization Results extraction ----

column_names_sts = {'S', 'I', 'D', 'T', 'H', 'E'};
column_names_par = {'alp', 'gam', 'dlt', 'sgm', 'tau'};

table_par = array2table(matrix_par', 'VariableNames', column_names_par);
table_sts = array2table(matrix_sts', 'VariableNames', column_names_sts);
results.par = table_par;
results.sts = table_sts;
save('AgeOpti-Results.mat', 'results');
