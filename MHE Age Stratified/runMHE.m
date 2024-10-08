function [res, par] = runMHEog(weightsvar, N_mhe, T_sim, ymeas_struct, c_struct, normcoef)

    % The following FUNCTION "runMHE_2" performs a MHE optimization on a reduced interval of 70 days (10 weeks)
    % in order to evaluate the residuals between simulation and measured data.
    % The active parameters which influence the results of the algorithm are the weights of the cost fucntion,
    % which are imposed as random variables changing every time the function is called. 

    % Compared to runMHE, runMHEog is a function running the ORIGINAL (homogeneous and deterministic) SIDTHE model

   % ------------------------------------------------------------------------------------------------------- %
    
    % INPUTS:
    %   N_mhe              - Estimation horizon (e.g., 21 for 3 weeks)
    %   T_sim              - Days over whcih the MHE algorithm is optimsed
    %   ymeas_struct       - Structure containing measured data
    %   c_struct           - Structure containing the contact matrices

    % FIXED/KNOWN VALUES INSIDE OF THE FUNCTION
    %   lambda = 1/10;     - Average recovery rate of the population
    %   Ts = 1;            - Integration timestep (1 day)

    % OUTPUT: 
    %   res                - Matrix of the MHE simulated STATES [6 (states) x 70 (time instants)]
    %   par                - Matrix of the MHE simulated PARAMETERS [20 (n parameters) x 140 (time instants)]

 % --------------------------------------------------------------------------------------------------------- %
   if T_sim < 70 || T_sim > 399
    error('Invalid value for N: %d. N must be between 70 and 399.', N);
   end
   
   Ts = 1;
   % ------------------------------------------------------------------------------------------------------- %
    
    %                        ----- CasADi states & dynamics initialization -----
    
    % "1" index referred to the population under 40
    s1 = casadi.SX.sym('s1',1,1);   alp1 = casadi.SX.sym('alp1',1,1);   v11 = casadi.SX.sym('v11',1,1);
    i1 = casadi.SX.sym('i1',1,1);   gma1 = casadi.SX.sym('gma1',1,1);   v12 = casadi.SX.sym('v12',1,1);
    d1 = casadi.SX.sym('d1',1,1);   dlt1 = casadi.SX.sym('dlt1',1,1);   v13 = casadi.SX.sym('v13',1,1);
    t1 = casadi.SX.sym('t1',1,1);   sgm1 = casadi.SX.sym('sgm1',1,1);   v14 = casadi.SX.sym('v14',1,1);
    h1 = casadi.SX.sym('h1',1,1);   tau1 = casadi.SX.sym('tau1',1,1);   v15 = casadi.SX.sym('v15',1,1);
    e1 = casadi.SX.sym('e1',1,1);                                       v16 = casadi.SX.sym('v16',1,1);
    
    % "2" index referred to the population between 40-60
    s2 = casadi.SX.sym('s2',1,1);   alp2 = casadi.SX.sym('alp2',1,1);   v21 = casadi.SX.sym('v21',1,1);
    i2 = casadi.SX.sym('i2',1,1);   gma2 = casadi.SX.sym('gma2',1,1);   v22 = casadi.SX.sym('v22',1,1);
    d2 = casadi.SX.sym('d2',1,1);   dlt2 = casadi.SX.sym('dlt2',1,1);   v23 = casadi.SX.sym('v23',1,1);
    t2 = casadi.SX.sym('t2',1,1);   sgm2 = casadi.SX.sym('sgm2',1,1);   v24 = casadi.SX.sym('v24',1,1);
    h2 = casadi.SX.sym('h2',1,1);   tau2 = casadi.SX.sym('tau2',1,1);   v25 = casadi.SX.sym('v25',1,1);
    e2 = casadi.SX.sym('e2',1,1);                                       v26 = casadi.SX.sym('v26',1,1);
    
    % "3" index referred to the population between 60-80
    s3 = casadi.SX.sym('s3',1,1);   alp3 = casadi.SX.sym('alp3',1,1);   v31 = casadi.SX.sym('v31',1,1);
    i3 = casadi.SX.sym('i3',1,1);   gma3 = casadi.SX.sym('gma3',1,1);   v32 = casadi.SX.sym('v32',1,1);
    d3 = casadi.SX.sym('d3',1,1);   dlt3 = casadi.SX.sym('dlt3',1,1);   v33 = casadi.SX.sym('v33',1,1);
    t3 = casadi.SX.sym('t3',1,1);   sgm3 = casadi.SX.sym('sgm3',1,1);   v34 = casadi.SX.sym('v34',1,1);
    h3 = casadi.SX.sym('h3',1,1);   tau3 = casadi.SX.sym('tau3',1,1);   v35 = casadi.SX.sym('v35',1,1);
    e3 = casadi.SX.sym('e3',1,1);                                       v36 = casadi.SX.sym('v36',1,1);
    
    % "4" index referred to the population over 80
    s4 = casadi.SX.sym('s4',1,1);   alp4 = casadi.SX.sym('alp4',1,1);   v41 = casadi.SX.sym('v41',1,1);
    i4 = casadi.SX.sym('i4',1,1);   gma4 = casadi.SX.sym('gma4',1,1);   v42 = casadi.SX.sym('v42',1,1);
    d4 = casadi.SX.sym('d4',1,1);   dlt4 = casadi.SX.sym('dlt4',1,1);   v43 = casadi.SX.sym('v43',1,1);
    t4 = casadi.SX.sym('t4',1,1);   sgm4 = casadi.SX.sym('sgm4',1,1);   v44 = casadi.SX.sym('v44',1,1);
    h4 = casadi.SX.sym('h4',1,1);   tau4 = casadi.SX.sym('tau4',1,1);   v45 = casadi.SX.sym('v45',1,1);
    e4 = casadi.SX.sym('e4',1,1);                                       v46 = casadi.SX.sym('v46',1,1); 

    lambda1 = 1/7; lambda2 = 1/12; lambda3 = 1/15; lambda4 = 1/15; 

    x = [   s1; i1; d1; t1; h1; e1; alp1; gma1; dlt1; sgm1; tau1; ...
            s2; i2; d2; t2; h2; e2; alp2; gma2; dlt2; sgm2; tau2; ...
            s3; i3; d3; t3; h3; e3; alp3; gma3; dlt3; sgm3; tau3; ...
            s4; i4; d4; t4; h4; e4; alp4; gma4; dlt4; sgm4; tau4      ];

    c = casadi.SX.sym('c',4,4); % contact matrix "c" time switching
    n_states = length(x);

   % ------------------------------------------------------------------------------------------------------- %
    %                      ----- States equations posing and computatio of the dynamics constraint ----

    eqns = [    
            %                      ----- set of equation for u40 -----
    
                -s1 * alp1 * ( c(1,1)*i1 + c(2,1)*i2 + c(3,1)*i3 + c(4,1)*i4 );...                          % dS1/dt
                 s1 * alp1 * ( c(1,1)*i1 + c(2,1)*i2 + c(3,1)*i3 + c(4,1)*i4 ) - (gma1+lambda1) * i1;...    % dI1/dt
                 i1 * gma1 - d1 * (lambda1 + dlt1);...                                                      % dD1/dt
                 dlt1 * d1 - (sgm1 + tau1) * t1;...                                                         % dT1/dt
                 (i1 + d1) * lambda1 + t1 * sgm1;...                                                        % dH1/dt
                 t1 * tau1;...                                                                              % dE1/dt
                 0;...
                 0;...
                 0;...
                 0;...
                 0;...         
                 
                 %                  ----- set of equation for 40-60 (middle aged) -----
                             
                -s2 * alp2 * ( c(1,2)*i1 + c(2,2)*i2 + c(3,2)*i3 + c(4,2)*i4 );...                          % dS2/dt
                 s2 * alp2 * ( c(1,2)*i1 + c(2,2)*i2 + c(3,2)*i3 + c(4,2)*i4 ) - (gma2+lambda2) * i2;...    % dI2/dt
                 i2 * gma2 - d2 * (lambda2 + dlt2);...                                                      % dD2/dt
                 dlt2 * d2 - (sgm2 + tau2) * t2;...                                                         % dT2/dt
                 (i2 + d2) * lambda2 + t2 * sgm2;...                                                        % dH2/dt
                 t2 * tau2;...                                                                              % dE2/dt
                 0;...
                 0;...
                 0;...
                 0;...
                 0;...     
    
                 %                  ----- set of equation for 60-80 (old aged) -----
                             
                -s3 * alp3 * ( c(1,3)*i1 + c(2,3)*i2 + c(3,3)*i3 + c(4,3)*i4 );...                          % dS3/dt
                 s3 * alp3 * ( c(1,3)*i1 + c(2,3)*i2 + c(3,3)*i3 + c(4,3)*i4 ) - (gma3+lambda3) * i3;...    % dI3/dt
                 i3 * gma3 - d3 * (lambda3 + dlt3);...                                                      % dD3/dt
                 dlt3 * d3 - (sgm3 + tau3) * t3;...                                                         % dT3/dt
                 (i3 + d3) * lambda3 + t3 * sgm3;...                                                        % dH3/dt
                 t3 * tau3;...                                                                              % dE3/dt
                 0;...
                 0;...
                 0;...
                 0;...
                 0;...  
    
    
                 %                  ----- set of equation for 80+ -----
                             
                -s4 * alp4 * ( c(1,4)*i1 + c(2,4)*i2 + c(3,4)*i3 + c(4,4)*i4 );...                          % dS4/dt
                 s4 * alp4 * ( c(1,4)*i1 + c(2,4)*i2 + c(3,4)*i3 + c(4,4)*i4 ) - (gma4+lambda4) * i4;...    % dI4/dt
                 i4 * gma4 - d4 * (lambda4 + dlt4);...                                                      % dD4/dt
                 dlt3 * d3 - (sgm4 + tau4) * t4;...                                                         % dT4/dt
                 (i4 + d4) * lambda4 + t4 * sgm4;...                                                        % dH4/dt
                 t4 * tau4;...                                                                              % dE4/dt
                 0;...
                 0;...
                 0;...
                 0;...
                 0
                 
                 ];
    
    f = casadi.Function('f',{x,c},{eqns});

   % ------------------------------------------------------------------------------------------------------- %

   %                              ----- Optmization Problem -----
    
    opti = casadi.Opti();
    X = opti.variable(n_states, N_mhe);           % Augmented state
    V = opti.variable(24,N_mhe-1);                % Process noise matrix
    C = opti.parameter(4,4*N_mhe);                % Contact matrix
    X_tilde = opti.parameter(6,4);                % Arrival Cost on the states (6 states, 4 age groups)
    P_tilde = opti.parameter(5,4);                % Arrival Cost on the params (5 params, 4 age groups)
    Ymeas_u40 = opti.parameter(6, N_mhe);
    Ymeas_mid = opti.parameter(6, N_mhe);
    Ymeas_old = opti.parameter(6, N_mhe);
    Ymeas_ger = opti.parameter(6, N_mhe);

    % LUT for the switching matrix for every time step
    c_LUT = [];
    for kk = 1:T_sim+N_mhe-1
        switch true
            case ismember(kk, 1:39)  || ismember(kk, 294:399)
                daily_c = (c_struct.home + c_struct.schl + c_struct.work + c_struct.othr)./normcoef';  
                
            case ismember(kk, 40:68) || ismember(kk, 242:293)
                daily_c = (c_struct.home + c_struct.work)./normcoef';   
                
            case ismember(kk, 69:115) || ismember(kk, 141:189)
                daily_c = (c_struct.home)./normcoef';  
                
            case ismember(kk, 116:140) || ismember(kk, 190:241)
                daily_c = (c_struct.home/2)./normcoef';  
        end
        c_LUT(:, :, kk) = daily_c;
    end

    %                  ---- Declarations of the weights of the cost function ----

    % Note that "weightsvar" is the actual optimization variable of the Bayesian Optimization

    % weightsAr= [  weightsvar{:, 1},  weightsvar{:, 2},  weightsvar{:, 3}, ...
    %               weightsvar{:, 4},  weightsvar{:, 5},  weightsvar{:, 6}, weightsvar{:, 7},  weightsvar{:, 8},...
    %               weightsvar{:, 9}, weightsvar{:, 10}, weightsvar{:, 11}   ];

    weightsAr= [  weightsvar{:, 1},...
                  weightsvar{:, 2},  weightsvar{:, 3},  weightsvar{:, 4}, weightsvar{:, 5},  weightsvar{:, 6},...
                  weightsvar{:, 7}, weightsvar{:, 8}, weightsvar{:, 9}   ];

    % Z1 = diag([weightsAr(:,1), weightsAr(:,2), weightsAr(:,1),...     % weight on the states consecutive estimate
    %            weightsAr(:,3), weightsAr(:,1), weightsAr(:,1)]);
    
    Z1 = diag([weightsAr(:,1), weightsAr(:,1), weightsAr(:,1),...     % weight on the states consecutive estimate
               weightsAr(:,1), weightsAr(:,1), weightsAr(:,1)]);
    
    % Z2 = diag([weightsAr(:,4), weightsAr(:,5)...                      % weight on the params consecutive estimate
    %            weightsAr(:,6), weightsAr(:,7), weightsAr(:,8)]);

    Z2 = diag([weightsAr(:,2), weightsAr(:,3)...                      % weight on the params consecutive estimate
               weightsAr(:,4), weightsAr(:,5), weightsAr(:,6)]);

    % Z3 = diag([weightsAr(:,9), weightsAr(:,10), weightsAr(:,11),...
    %            weightsAr(:,11), weightsAr(:,9), weightsAr(:,11)]);   % weight on the process noise on the ODE states dynamics
 
    Z3 = diag([weightsAr(:,7), weightsAr(:,8), weightsAr(:,9),...     % weight on the params consecutive estimate
               weightsAr(:,9), weightsAr(:,7), weightsAr(:,9)]);

    % Constraints Setting and Solver Options
    rr = 1;
    for jj=1:N_mhe-1 
    
        k1 = f(X(:,jj),C(:,rr:rr+3));
        k2 = f(X(:,jj)+Ts/2*k1,C(:,rr:rr+3));
        k3 = f(X(:,jj)+Ts/2*k2,C(:,rr:rr+3));
        k4 = f(X(:,jj)+Ts*k3,C(:,rr:rr+3));
        x_plus = Ts/6*(k1+2*k2+2*k3+k4);
        opti.subject_to(X(:,jj+1)== X(:,jj) + x_plus + [ V(1:6,jj); zeros(5,1);   V(7:12,jj);  zeros(5,1);...
                                                         V(13:18,jj); zeros(5,1); V(19:24,jj); zeros(5,1) ]);
    
        rr = rr+4; % 
    end

    opti.subject_to(X(:) > 0); 
    
    for ii=1:N_mhe
        opti.subject_to( sum([X(1:6,ii); X(12:17,ii); X(23:28,ii); X(34:39,ii)]) == 1 )
    end
    
    %                                       ---- Disturbances bounds ----
    std_u40 = std(ymeas_struct.u40,1,2);
    std_mid = std(ymeas_struct.mid,1,2);
    std_old = std(ymeas_struct.old,1,2);
    std_ger = std(ymeas_struct.ger,1,2);
    
    w = 0.005;
    
    opti.subject_to(-std_u40(1)*w < V(1,:) < std_u40(1)*w);   opti.subject_to(-std_mid(1)*w < V(7,:)  < std_mid(1)*w);
    opti.subject_to(-std_u40(2)*w < V(2,:) < std_u40(2)*w);   opti.subject_to(-std_mid(2)*w < V(8,:)  < std_mid(2)*w);
    opti.subject_to(-std_u40(3)*w < V(3,:) < std_u40(3)*w);   opti.subject_to(-std_mid(3)*w < V(9,:)  < std_mid(3)*w);
    opti.subject_to(-std_u40(4)*w < V(4,:) < std_u40(4)*w);   opti.subject_to(-std_mid(4)*w < V(10,:) < std_mid(4)*w);
    opti.subject_to(-std_u40(5)*w < V(5,:) < std_u40(5)*w);   opti.subject_to(-std_mid(5)*w < V(11,:) < std_mid(5)*w);
    opti.subject_to(-std_u40(6)*w < V(6,:) < std_u40(6)*w);   opti.subject_to(-std_mid(6)*w < V(12,:) < std_mid(6)*w);
    
    opti.subject_to(-std_old(1)*w < V(13,:) < std_old(1)*w);   opti.subject_to(-std_ger(1)*w < V(19,:) < std_ger(1)*w);
    opti.subject_to(-std_old(2)*w < V(14,:) < std_old(2)*w);   opti.subject_to(-std_ger(2)*w < V(20,:) < std_ger(2)*w);
    opti.subject_to(-std_old(3)*w < V(15,:) < std_old(3)*w);   opti.subject_to(-std_ger(3)*w < V(21,:) < std_ger(3)*w);
    opti.subject_to(-std_old(4)*w < V(16,:) < std_old(4)*w);   opti.subject_to(-std_ger(4)*w < V(22,:) < std_ger(4)*w);
    opti.subject_to(-std_old(5)*w < V(17,:) < std_old(5)*w);   opti.subject_to(-std_ger(5)*w < V(23,:) < std_ger(5)*w);
    opti.subject_to(-std_old(6)*w < V(18,:) < std_old(6)*w);   opti.subject_to(-std_ger(6)*w < V(24,:) < std_ger(6)*w);
    
    
    %                                       ---- Parameters Bounds ----
    % alpha bounds 
    opti.subject_to( 0.05  <=  X(7,:)  <= 100 );     opti.subject_to( 0.05  <= X(29,:)  <= 100 );
    opti.subject_to( 0.005 <= X(18,:)  <= 100 );     opti.subject_to( 0.05  <= X(40,:)  <= 100 );
    
    % gamma bounds 
    opti.subject_to( 0.005 <= X(8,:)   <= 0.6 );     opti.subject_to( 0.005 <= X(30,:)  <= 0.6 );
    opti.subject_to( 0.005 <= X(19,:)  <= 0.6 );     opti.subject_to( 0.005 <= X(41,:)  <= 0.6 );
    
    % delta bounds 
    opti.subject_to( 1e-4  <= X(9,:)   <= 0.6 );     opti.subject_to( 1e-4  <= X(31,:)  <= 0.6 );
    opti.subject_to( 1e-4  <= X(20,:)  <= 0.6 );     opti.subject_to( 1e-4  <= X(42,:)  <= 0.6 );
    
    % sigma bounds
    opti.subject_to( 1e-4  <= X(10,:)  <= 0.5 );     opti.subject_to( 1e-4  <= X(32,:)  <= 0.5 ); 
    opti.subject_to( 1e-4  <= X(21,:)  <= 0.5 );     opti.subject_to( 1e-4  <= X(43,:)  <= 0.5 ); 
    
    % tau bounds
    opti.subject_to( 1e-4  <= X(11,:)  <= 0.5 );     opti.subject_to( 1e-4  <= X(33,:)  <= 0.5 ); 
    opti.subject_to( 1e-4  <= X(22,:)  <= 0.5 );     opti.subject_to( 1e-4  <= X(44,:)  <= 0.5 ); 
    
    
    %                                        ---- Objective function ----
    
    obj = sumsqr( Z1*(X(1:6,:)   - X_tilde(:,1)) ) + sumsqr( Z2*(X(7:11,:)  - P_tilde(:,1)) ) + sumsqr( Z3*(Ymeas_u40 - X(1:6,:))./max(ymeas_struct.u40, [], 2) ) + ... % u40 cost 
          sumsqr( Z1*(X(12:17,:) - X_tilde(:,2)) ) + sumsqr( Z2*(X(18:22,:) - P_tilde(:,2)) ) + sumsqr( Z3*(Ymeas_mid - X(12:17,:))./max(ymeas_struct.mid, [], 2) ) + ... % middle aged cost 
          sumsqr( Z1*(X(23:28,:) - X_tilde(:,3)) ) + sumsqr( Z2*(X(29:33,:) - P_tilde(:,3)) ) + sumsqr( Z3*(Ymeas_old - X(23:28,:))./max(ymeas_struct.old, [], 2) ) + ... % senior cost
          sumsqr( Z1*(X(34:39,:) - X_tilde(:,4)) ) + sumsqr( Z2*(X(40:44,:) - P_tilde(:,4)) ) + sumsqr( Z3*(Ymeas_ger - X(34:39,:))./max(ymeas_struct.ger, [], 2) ) + ... % geriatric cost
          sumsqr(V);

    opti.minimize(obj);
    
    % Solver Options
    opts.print_time = 0;
    opts.ipopt.print_level = 0;
    opti.solver('ipopt', opts);


   % ------------------------------------------------------------------------------------------------------- %

    %                                     ---- Moving Horizon FOR cycle ----
   
    matrix_sts = []; % matrix to store simulation results (states)
    matrix_par = []; % matrix to store simulation results (states)
 
    for ii = N_mhe:1:N_mhe+T_sim-1
        % Matrix switching settings in the dynamics
        rr = 1;
        for jj = 1:N_mhe
        opti.set_value(C(:,rr:rr+3), c_LUT(:,:,ii-N_mhe+jj))
        rr = rr+4;
        end
    
        if ii == N_mhe   
            % First iteration "arrival cost parameters" setting
            opti.set_value(X_tilde, [ ymeas_struct.u40(:,1) ymeas_struct.mid(:,1) ymeas_struct.old(:,1) ymeas_struct.ger(:,1) ]);    
            opti.set_value(P_tilde, ( [ 2      1.5    1.5   1.5; ...
                                        0.3    0.13   0.12  0.12; ...
                                        0.002  0.004  0.01  0.01; ...
                                        0.04   0.04   0.02  0.01; ...
                                        0.002  0.01   0.01  0.02] ));    
        else 
            % Every other iteration "arrival cost parameters"
            opti.set_value(X_tilde, currentState)   
            opti.set_value(P_tilde, currentParam)
        end
    
        % measurement update 
        opti.set_value(Ymeas_u40, ymeas_struct.u40(:,ii-N_mhe+1:ii));
        opti.set_value(Ymeas_mid, ymeas_struct.mid(:,ii-N_mhe+1:ii));
        opti.set_value(Ymeas_old, ymeas_struct.old(:,ii-N_mhe+1:ii));
        opti.set_value(Ymeas_ger, ymeas_struct.ger(:,ii-N_mhe+1:ii));
    
        sol = opti.solve();
        
        % Save of the optimization values for every different group (PARAMETERS)
        params_dyn1 = opti.value(X(7:11, end));     params_dyn3 = opti.value(X(29:33,end));
        params_dyn2 = opti.value(X(18:22,end));     params_dyn4 = opti.value(X(40:44,end));
    
        % Save of the optimization values for every different group (STATES)
        states_dyn1 = opti.value(X(1:6,  end));      states_dyn3 = opti.value(X(23:28,end));
        states_dyn2 = opti.value(X(12:17,end));      states_dyn4 = opti.value(X(34:39,end));
    
        currentState = [states_dyn1, states_dyn2, states_dyn3, states_dyn4];
        currentParam = [params_dyn1, params_dyn2, params_dyn3, params_dyn4];
    
        matrix_sts = [matrix_sts; [states_dyn1', states_dyn2', states_dyn3', states_dyn4'] ];
        matrix_par = [matrix_par; [params_dyn1', params_dyn2', params_dyn3', params_dyn4'] ];
    
    end

    % First output computation: difference between data and estimated states
    res = matrix_sts;

    % Second output computation: difference between consecutive estimated parameters
    par = matrix_par;

end
