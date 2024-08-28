function objFcn = bayesMHEObj(N_mhe, T_sim, ymeas_struct, c_struct, N_mc)
    
    % The following FUNCTION "bayesMHEObj" outputs the Bayesian Optimization funtion considering 
    % the MHE optimization results over a certain number of N Montecarlo simulations

   % ------------------------------------------------------------------------------------------------------- %
    
    % INPUTS:
    %   N_mhe              - Estimation horizon (e.g., 21 for 3 weeks)
    %   T_sim              - Days over whcih the MHE algorithm is optimsed
    %   ymeas_struct       - Structure containing measured data per age groups
    %   c_struct           - Structure containing the contact matrices
    %   N_mc               - Number of Montecarlo simulations one wants to perform

    % FIXED/KNOWN VALUES INSIDE OF THE FUNCTION
    %   x = 24             - total number of states (6 states x 4 age groups = 24 total states)
    %   p = 20             - total number of parameters (5 parameters x 4 age groups = 20 total parameters)

    % OUTPUT: 
    %   objFcn             - Objective function (numerical value)

    % --------------------------------------------------------------------------------------------------------- %
    x = 24;           
    p = 20;
    % --------------------------------------------------------------------------------------------------------- %

    e_mat = zeros(x, T_sim, N);    % initialization of the 3D matrix to store Montecarlo simulation errors
    dltp_mat = zeros(p, T_sim-1, N); % initialization of the 3D matrix to store Montecarlo simulation errors

    % Montecarlo simulations of the MHE algorithm over N_sim days (see "runMHE" fcn)
    for j = 1:N_mc
        [e_mat(:,:,j), dltp_mat(:,:,j)] = runMHE(N_mhe, T_sim, ymeas_struct, c_struct);
    end

    % --------------------------------------------------------------------------------------------------------- %
    
    % Evaluation of the variables inside of the cost function
    e_bar = mean(e_mat, 3);        % mean error for each time instant across all experiments (states)
    dltp_bar = mean(dltp_mat,3);   % mean error for each time instant across all experiments (parameters)
    
    e_tilde = mean(e_bar, 1);            % mean error over the time instats T_sim
    dltp_tilde = mean(dltp_bar, 1);      % mean error over the time instats T_sim

    objFcn = abs(log(e_tilde/m)) + abs(log(S_tilde/(2*m)));
end