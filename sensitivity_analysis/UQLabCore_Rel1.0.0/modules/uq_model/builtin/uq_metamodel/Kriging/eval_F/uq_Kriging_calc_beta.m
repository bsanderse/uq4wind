function beta = uq_Kriging_calc_beta( ...
                        F, trend_type, Y, beta_estim_method, aux_matrices )
% UQ_KRIGING_CALC_BETA: computes the trend coefficients beta of a Kriging
% metamodel
%
% See also UQ_KRIGING_CALCULATE, UQ_KRIGING_EVAL_J_OF_THETA_ML

if isrow(Y)
    Y = transpose(Y) ;
end

if strcmpi(trend_type,'simple')
    % In simple kriging no regression takes place
    beta = ones(size(F,2),1);
    return
end

% In case of QR method make sure that all the required matrices are
% available, otherwise switch to standard method
if strcmpi(beta_estim_method, 'qr') && ...
        (~isfield(aux_matrices,'Q1') || isempty(aux_matrices.Q1))
        beta_estim_method = 'standard';
end


switch lower(beta_estim_method)
    case 'qr'
        % If the QR decomposition of the decorelated F(Ftilde) has already been
        % calculated use it to find beta
        Q1 = aux_matrices.Q1 ;
        G = aux_matrices.G ;
        Ytilde = aux_matrices.Ytilde ;
        if size(G,1)~= size(G,2) || rcond(G) < 1e-10
            beta = pinv(G) * transpose(Q1) * Ytilde ;
        else
            beta = G \ transpose(Q1) * Ytilde ;
        end
    case 'standard'
    % Calculate beta using the standard approach
        % retrieve auxiliary stored matrices
    FTRinv = aux_matrices.FTRinv ;

    if isempty(aux_matrices.FTRinvF_inv)
        FTRinvF = aux_matrices.FTRinvF ;
        beta = (FTRinvF \ FTRinv) * Y ;
    else
        FTRinvF_inv = aux_matrices.FTRinvF_inv ;
        beta = FTRinvF_inv * FTRinv * Y ;
    end    
end


    








    


