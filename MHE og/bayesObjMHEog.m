function objFcn = bayesObjMHEog(weightsvar, N_mhe, T_sim, ymeas)
    
    % The following FUNCTION "bayesObjMHEog" is a modification of "bayesMHEObj" for a homogeneous model

   % ------------------------------------------------------------------------------------------------------- %
    
    % INPUTS:
    %   weightsvar         - Struct containing the weights Z1, Z2, Z3 for the MHE simulation
    %   N_mhe              - Estimation horizon (e.g., 21 for 3 weeks)
    %   T_sim              - Days over whcih the MHE algorithm is optimsed
    %   ymeas              - Structure containing measured data per age groups

    % OUTPUT: 
    %   objFcn             - Objective function (numerical value)

   % ------------------------------------------------------------------------------------------------------- %
   
   % Run of the MHE algorithm
    [sts_mat, ~ ] = runMHEog(weightsvar, N_mhe, T_sim, ymeas);

    % modification of the 'y_meas' struct to be used into the obj function
    y_meas = ymeas(:,N_mhe:N_mhe+T_sim-1);

    e = abs(y_meas - sts_mat); % residuals calculation (estimation error)
    r = []; 
    N = T_sim;

    alpha = 0.05;                       % significance level (95%)
    df = N;                             % Degrees of freedom
    z_alpha = chi2inv(1 - alpha, df);   % Chi-squared threshold value


    %                      ----- Construction of estimated covariance function  ----- %

    for tau = 1:N
        r_hat = 0;              % Initialize sum for the current tau "r_hat" value
        for t = 1:N - tau
            r_hat = r_hat + e(:,t + tau).* e(:,t);
        end
        r(tau,:) = r_hat ./ N;  % Store autocorrelation value for current tau
    end

    r_hat0 = sum(e.^2, 2)./N;   % value of the autocorrelation when tau = 0

    %                     ----- Construction of test quantity NrTr/r_hat(0)^2 (Distributed as Chi-squared)----- %
    
    testq = N .* sum(r.^2, 1)./(r_hat0').^2;


    %                     ----- Construction of Cost Function For Bayes ----- %

    % penaltyTerm = ( testq - z_alpha) .* max(0, (testq - z_alpha).^3 );
    penaltyTerm = ( testq - z_alpha) .* max(0, (testq - z_alpha) );
    objFcn = sum( penaltyTerm );

end