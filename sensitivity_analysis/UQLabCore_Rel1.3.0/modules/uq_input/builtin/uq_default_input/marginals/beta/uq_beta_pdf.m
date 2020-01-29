function F = uq_beta_pdf( X, parameters )
% F = UQ_BETA_PDF(X, parameters):
%     calculates the Probability Density Function values of samples X  that 
%     follow a Beta distribution with parameters specified in the vector
%     'parameters'
%
% Input parameters:
%     'X'             The samples
%     'parameters'    [r s a b] where [a b] are the bounds of the
%                     distribution and r,s the other required parameters.
%
%     NOTE: 
%     Please keep in mind the differences between the MATLAB and UQLab
%     implementation
%     MATLAB : beta distribution always has support [0,1] 
%     UQLab  : generalised case with arbitrary support [a,b] 
% 
% See also: BETAPDF, UQ_EVAL_JACOBI

r = parameters(1);
s = parameters(2);

switch length(parameters)
    case 2 
        % the support is [0,1]
        X_rescaled = X;
        coeff = 1;
    case 4
        % the support is [a,b]
        a = parameters(3);
        b = parameters(4);
        X_rescaled = (X - a)/(b-a);
        coeff = 1/(b-a);
    otherwise
        error('Error: Beta distribution can only work with 2 or 4 parameters defined!')
end

F = coeff*betapdf(X_rescaled,  r, s);

