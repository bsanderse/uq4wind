function auxMatrices = uq_Kriging_calc_auxMatrices( R, F, Y, run_case )
% UQ_KRIGING_CALC_AUXMATRICES: calculates some auxiliary matrices that can be
% useful to have in order to speed up some calculations during training of a
% Kriging metamodel as well as during predictions
%
% auxMatrices = UQ_KRIGING_CALC_AUXMATRICES( R, F, Y, run_case )
% Given a correlation matrix R, and trend matrix F the following auxiliary
% matrices are returned inside the auxMatrices structure:
%   .cholR      : the Cholesky decomposition of R if exists otherwise NaN
%   .Rinv       : the pseudo inverse of R in case the Cholesky decomposion
%                 of it is not available due to R being ill conditioned. If
%                 the Cholesky decomposition is available then this is
%                 empty
%   .FTRinv     : the quantity F' * R^(-1)
%   .FTRinvF    : the quantity F' * R^(-1) * F
%   .FTRinvF_inv: the pseudo-inverse (F' * R^(-1) * F)^(-1) in case 
%                 F' * R^(-1) * F is ill-conditioned. Otherwise this field
%                 is empty
% ** aditional matrices that may be present in case of ML estimation
% method:
%   .Q          : (see documentation for description of this quantity)
%   .G          : ( -//-)
%   .Ytilde     : ( -//-)
% the above fields will be empty in case of ill-conditioned R matrix.
%
% See also UQ_KRIGING_CALC_LEAVEKOUT, UQ_KRIGING_EVAL_J_OF_THETA_ML, 
%  UQ_KRIGING_CALCULATE

%1) try to calculate the Cholesky decomposition of R. This will speed up
%   the calculation of the regression coefficients, beta, BUT it is only
%   possible given that R is not ill-conditioned
try
    L = chol(R);
    auxMatrices.cholR = L;
    auxMatrices.Rinv = [];
    FTRinv = (transpose(F) / L ) / transpose(L) ;
    if any(strcmpi(run_case, {'ml_optimization', 'ml_estimation'}))
        T = transpose(L);
        Ytilde = T \ Y ;
        [Q1,G] = qr(T \ F,0) ;
        auxMatrices.Q1 = Q1;
        auxMatrices.G = G;
        auxMatrices.Ytilde = Ytilde;
        auxMatrices.Ftilde = T \ F ;
        % this will occur when this function is called during optimization
        % of the ML objective function.
        % in this case no additional matrices are required so just return
        if strcmpi(run_case, 'ml_optimization')
            return
        end
    end
catch
    % The Cholesky decomposition failed
    auxMatrices.cholR = nan ;
    if strcmpi(run_case, 'ml_estimation')
        % in case of ML estimation method set QR-relevant matrices to empty
        % 'signalling' the ill-conditioned state of Chol(R)
        % this makes sure that subsequent calculations will avoid using the
        % QR decomposition results. 
        auxMatrices.Q1 = [];
        auxMatrices.G = [];
        auxMatrices.Ytilde = [];
        auxMatrices.Ftilde = [];
    end  
    
    % if R is ill conditioned store its pseudo-inverse instead
    Rinv = pinv(R);
    auxMatrices.Rinv = Rinv ;
    FTRinv = transpose(F)* Rinv ;

end
auxMatrices.FTRinv = FTRinv ;
FTRinvF = FTRinv * F;
auxMatrices.FTRinvF = FTRinvF ;
%2) try to calculate the Cholesky decomposition of F'RinvF. This will speed up
%   the calculation of the regression coefficients, beta, BUT it is only
%   possible given that F'RinvF is not ill-conditioned
if rcond(FTRinvF) > 1e-10
    auxMatrices.FTRinvF_inv = [];
else
    auxMatrices.FTRinvF_inv = pinv(FTRinvF) ;
end


