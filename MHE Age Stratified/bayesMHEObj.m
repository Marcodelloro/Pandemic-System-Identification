function objFcn = bayesMHEObj(weightsvar, N_mhe, T_sim, ymeas_struct, c_struct)
    
    % The following FUNCTION "bayesMHEObj" outputs the Bayesian Optimization funtion considering 
    % the MHE optimization results over a certain number of N Montecarlo simulations

   % ------------------------------------------------------------------------------------------------------- %
    
    % INPUTS:
    %   weightsvar         - Struct containing the weights Z1, Z2, Z3 for the MHE simulation
    %   N_mhe              - Estimation horizon (e.g., 21 for 3 weeks)
    %   T_sim              - Days over whcih the MHE algorithm is optimsed
    %   ymeas_struct       - Structure containing measured data per age groups
    %   c_struct           - Structure containing the contact matrices

    % OUTPUT: 
    %   objFcn             - Objective function (numerical value)

   % ------------------------------------------------------------------------------------------------------- %
   
   % Run of the MHE algorithm
    [sts_mat, ~ ] = runMHE(weightsvar, N_mhe, T_sim, ymeas_struct, c_struct);

    % modification of the 'y_meas' struct to be used into the obj function
    y_meas = [ymeas_struct.u40' ymeas_struct.mid' ymeas_struct.old' ymeas_struct.ger'];
    y_meas = y_meas(N_mhe:N_mhe+T_sim-1,:);

    e = abs(y_meas - sts_mat); % residuals calculation (estimation error)
    r = []; 
    N = length(y_meas);

    alpha = 0.05;                       % significance level (95%)
    df = N;                             % Degrees of freedom
    z_alpha = chi2inv(1 - alpha, df);   % Chi-squared threshold value


    %                      ----- Construction of estimated covariance function  ----- %

    for tau = 1:N
        r_hat = 0;              % Initialize sum for the current tau "r_hat" value
        for t = 1:N - tau
            r_hat = r_hat + e(t + tau,:).* e(t,:);
        end
        r(tau,:) = r_hat ./ N;  % Store autocorrelation value for current tau
    end

    r_hat0 = sum(e.^2, 1)./N;   % value of the autocorrelation when tau = 0


    %                     ----- Construction of test quantity NrTr/r_hat(0)^2 (Distributed as Chi-squared)----- %

    testq = N .* sum(r.^2, 1)./(r_hat0).^2;


    %                     ----- Construction of Cost Function For Bayes ----- %

    objFcn = sumsqr( testq - ones(1,length(testq))*z_alpha );

end