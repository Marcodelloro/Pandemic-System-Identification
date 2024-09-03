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
    [sts_mat, par_mat] = runMHE(weightsvar, N_mhe, T_sim, ymeas_struct, c_struct);

    % modification of the 'y_meas' struct to be used into the obj function
    y_meas = [ymeas_struct.u40' ymeas_struct.mid' ymeas_struct.old' ymeas_struct.ger'];
    y_meas = y_meas(N_mhe:N_mhe+T_sim-1,:);

    % cost function computation 
    objFcn = sumsqr(sts_mat - y_meas) + sumsqr( diff(par_mat)./max(abs(diff(par_mat)), [], 1) );
end