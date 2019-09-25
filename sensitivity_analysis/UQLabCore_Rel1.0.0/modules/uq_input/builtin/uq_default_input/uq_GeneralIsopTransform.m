function Y = uq_GeneralIsopTransform(X, X_Marginals, X_Copula, Y_Marginals, Y_Copula)
% Y = UQ_GENERALISOPTRANSFORM(X, X_Marginals, X_Copula, Y_Marginals, Y_Copula):
%     Allows to map samples of random vectors of arbitrary distributions 
%     to a different space. Both spaces need to be specified via the 
%     Marginals and Copula formalism as described in the UQLab Input module
%     user manual.
%
% See also:  UQ_INVNATAFTRANSFORM, UQ_ISOPTRANSFORM

M = size(X,2);
Y = zeros(size(X));

switch lower(X_Copula.Type) 
    %% X is described by independent copula  
    case 'independent'
        switch lower(Y_Copula.Type) % and independent output as well
            case 'independent'
                Y = uq_IsopTransform(X, X_Marginals, Y_Marginals);
            case 'gaussian' 
                % Specify the U marginals
                [UMarginals(1:M).Type] = deal('Gaussian');
                [UMarginals(1:M).Parameters] = deal([0 1]);
                % Map X to the standard normal space, i.e. get U
                U = uq_IsopTransform(X, X_Marginals, UMarginals);
                % Finally map U to Y (the output space)
                Y = uq_invNatafTransform(U, Y_Marginals, Y_Copula);
        end
    case 'gaussian' 
        %% X is described by Gaussian copula  
        % First map the samples to the standard normal space via the Nataf transform
        U = uq_NatafTransform(X, X_Marginals, X_Copula);
        % Specify the U marginals
        [UMarginals(1:M).Type] = deal('Gaussian');
        [UMarginals(1:M).Parameters] = deal([0 1]);
        
        % Finally map the samples to the output space
        switch lower(Y_Copula.Type)
            case 'independent'  
                Y = uq_IsopTransform(U, UMarginals, Y_Marginals);
            case 'gaussian'
                Y = uq_invNatafTransform(U, Y_Marginals, Y_Copula);
        end
    otherwise
        error('Error: the specified copula type "%s" cannot be handled. Please check the definition of your marginals.', X_Copula.Type)
end
