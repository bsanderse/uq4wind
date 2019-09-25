function sigmaSq = uq_Kriging_calc_sigmaSq( parameters, estim_method )
% UQ_KRIGING_CALC_SIGMASQ: returns the variance of the Guassian process
% while trying to calculate it in the most efficient way depending on the
% hyperparameters estimation method and in case of ML, the condition of the
% correlation matrix R
%
% See also UQ_KRIGING_CALCULATE, UQ_KRIGING_EVAL_J_OF_THETA_ML

switch lower(estim_method)
    case 'cv'
        sigmaSq =  sum(cell2mat(parameters.errors)./cell2mat(parameters.sigma_i_Sq)) ...
            / parameters.CV_K ;
    case 'ml_chol'
          N = parameters.N ;
          beta = parameters.beta ;
          Ytilde = parameters.Ytilde;
          Ftilde = parameters.Ftilde;
          
          z = Ytilde - Ftilde*beta ;
          sigmaSq = transpose(z) * z / N; 
    case 'ml_bypass_chol'
        % bypass the calculation of the regression coefficients
        N = parameters.N;
        Q1 = parameters.Q1;
        Ytilde = parameters.Ytilde;
        
        z = Ytilde - Q1*transpose(Q1)* Ytilde ;
        sigmaSq = transpose(z) * z / N;
     case {'ml_nochol', 'ml_bypass_nochol'}       
        % calculation of beta cannot be avoided in that case
        F = parameters.F;
        Y = parameters.Y;
        N = parameters.N ;
        
        beta = parameters.beta ;
        z = Y - F* beta ;
        sigmaSq = transpose(z) * parameters.Rinv * z  /N;
end
