function J = uq_Kriging_eval_J_of_theta_CV(theta, parameters)
% UQ_KRIGING_EVAL_J_OF_THETA_CV : Leave-K-Out Cross Validation objective
% function
%
% See also: UQ_KRIGING_OPTIMIZER, UQ_KRIGING_EVAL_J_OF_THETA_ML,
% UQ_KRIGING_CALCULATE

X = parameters.X;
CorrOptions = parameters.CorrOptions;
evalR_handle = CorrOptions.Handle ;
CV_K = parameters.CV_K;

try
    R = evalR_handle( X, X, theta, CorrOptions );
    [ errors] = uq_Kriging_calc_leaveKout( CV_K, parameters.N, ...
        parameters.Y, R, parameters.F, 'default' ) ;
    
    J = sum(cell2mat(errors)) / CV_K ;
    
catch 
    % If something goes wrong (typically ill-conditioned R) return a
    % 'large' value
    J = realmax ; 
end






