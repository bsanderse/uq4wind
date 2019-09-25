function J = uq_Kriging_eval_J_of_theta_ML(theta, parameters)
% UQ_KRIGING_EVAL_J_OF_THETA_ML : Maximum Likelihood objective
% function
%
% See also: UQ_KRIGING_OPTIMIZER, UQ_KRIGING_EVAL_J_OF_THETA_CV,
% UQ_KRIGING_CALCULATE

X = parameters.X;
Y = parameters.Y;
N = parameters.N;
CorrOptions = parameters.CorrOptions;
evalR_handle = CorrOptions.Handle ;
F = parameters.F;
trend_type = parameters.trend_type;

try
    % calculate R and auxiliary matrices
    R = evalR_handle( X, X, theta, CorrOptions );
    auxMatrices = uq_Kriging_calc_auxMatrices( R, F, Y, 'ml_optimization' );

    if ~isnan(auxMatrices.cholR)
        logDetR = 2* sum(log(diag(auxMatrices.cholR))) ;
        parameters.Q1 = auxMatrices.Q1 ;
        parameters.Ytilde = auxMatrices.Ytilde ;
        sigmaSq = uq_Kriging_calc_sigmaSq( parameters, 'ml_bypass_chol' ) ;
    else
        eps = 1e-320; 
        logDetR = log(max(det(R), eps));
        
        parameters.beta = uq_Kriging_calc_beta( ...
                        F, trend_type, Y, 'standard', auxMatrices );
        parameters.Rinv = auxMatrices.Rinv ;
        sigmaSq = uq_Kriging_calc_sigmaSq( parameters, 'ml_bypass_NOchol' ) ;
    end
    
    J = 0.5*(N * log(2*pi*sigmaSq) + logDetR + N) ;
catch
    J = realmax ;
end