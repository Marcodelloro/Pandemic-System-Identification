function [res, par] = runMHEog(weightsvar, N_mhe, T_sim, ymeas)

    % The following FUNCTION "runMHE_2" performs a MHE optimization on a reduced interval of 70 days (10 weeks)
    % in order to evaluate the residuals between simulation and measured data.
    % The active parameters which influence the results of the algorithm are the weights of the cost fucntion,
    % which are imposed as random variables changing every time the function is called. 

    % Compared to runMHE, runMHEog is a function running the ORIGINAL (homogeneous and deterministic) SIDTHE model

   % ------------------------------------------------------------------------------------------------------- %
    
    % INPUTS:
    %   N_mhe              - Estimation horizon (e.g., 21 for 3 weeks)
    %   T_sim              - Days over whcih the MHE algorithm is optimsed
    %   ymeas              - Matrix containing the available compartments data
    %   c_struct           - Structure containing the contact matrices

    % FIXED/KNOWN VALUES INSIDE OF THE FUNCTION
    %   lambda = 1/10;     - Average recovery rate of the population
    %   Ts = 1;            - Integration timestep (1 day)

    % OUTPUT: 
    %   res                - Matrix of the MHE simulated STATES [6 (states) x 70 (time instants)]
    %   par                - Matrix of the MHE simulated PARAMETERS [5 (parameters) x 70 (time instants)]

 % --------------------------------------------------------------------------------------------------------- %
   if T_sim < 70 || T_sim > 399
    error('Invalid value for N: %d. N must be between 70 and 399.', N);
   end
   
   Ts = 1;
   % ------------------------------------------------------------------------------------------------------- %
    
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

    % Note that "weightsvar" is the actual optimization variable of the Bayesian Optimization

    weightsAr= [  weightsvar{:, 1},  weightsvar{:, 2},  weightsvar{:, 3}, ...
                  weightsvar{:, 4},  weightsvar{:, 5},  weightsvar{:, 6}, weightsvar{:, 7},  weightsvar{:, 8},...
                  weightsvar{:, 9}, weightsvar{:, 10}, weightsvar{:, 11}   ];

  
    Z1 = diag([weightsAr(:,1), weightsAr(:,2), weightsAr(:,1),...     % weight on the states consecutive estimate
               weightsAr(:,3), weightsAr(:,1), weightsAr(:,1)]);
    
    Z2 = diag([weightsAr(:,4), weightsAr(:,5)...                      % weight on the params consecutive estimate
               weightsAr(:,6), weightsAr(:,7), weightsAr(:,8)]);

    Z3 = diag([weightsAr(:,9), weightsAr(:,10), weightsAr(:,11),...
               weightsAr(:,11), weightsAr(:,9), weightsAr(:,11)]);   % weight on the process noise on the ODE states dynamics
 

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
         opti.subject_to( sum(X(1:6,ii)) == 1 ) % constraint on the I - D compartment
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
 
    for ii = N_mhe:1:N_mhe+T_sim-1
    
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

    % First output computation: difference between data and estimated states
    res = matrix_sts;

    % Second output computation: difference between consecutive estimated parameters
    par = matrix_par;

end
