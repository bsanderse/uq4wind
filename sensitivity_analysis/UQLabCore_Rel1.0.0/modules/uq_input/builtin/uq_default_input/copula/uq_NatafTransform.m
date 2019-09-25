function U = uq_NatafTransform( Xin, marginals, copula )
% U = UQ_NATAFTRANSFORM(Xin, marginals, copula):
%     performs (Generalised) Nataf Transformation of samples (X) 
%     from the physical to the standard normal space (U)
%
% References:
%
%   Lebrun, R. and Dutfoy (2009), A., A generalization of the Nataf 
%   transformation to distributions with elliptical copula, 
%   Prob. Eng. Mech. 24(2), 172-178
%
% See also: UQ_INVNATAFTRANSFORM

%% Get the non-constant components of the random vector
% find the indices of the non-constant components
Types = {marginals(:).Type};
indConst = strcmpi(Types, 'constant'); 
indNonConst = ~indConst ;
% store non-constant components in X
X = Xin(:,indNonConst) ;

% NOTE: 
% In the algorithm that follows the symbols from (Lebrun and Dutfoy, 2009)
% are mostly adopted

%% Do isoprobabilistic transformation from X to V
% V denotes the samples in the standard elliptical space where covariance is still existent
switch lower(copula.Type)
    case 'gaussian'
        [vMarginals(1:length(marginals(indNonConst))).Type] = deal('gaussian') ;
        [vMarginals(1:length(marginals(indNonConst))).Parameters] = deal([0 1]) ;
        V = uq_IsopTransform(X, marginals(indNonConst), vMarginals);
    case 'student'
        [vMarginals(1:length(marginals(indNonConst))).Type] = deal('student') ;
        [vMarginals(1:length(marginals(indNonConst))).Parameters] = deal(1) ;
        V = uq_IsopTransform(X, marginals(indNonConst), vMarginals);
    otherwise
        err('Error: Unknown type of elliptical copula!')
        return;
end
%% Perform step 3: mapping from V to U
% if the Cholesky decomposition of the correlation matrix already exists
% use it, otherwise calculate it here (and store it)
if isfield(copula,'cholR') && ~isempty(copula.cholR)
    L = copula.cholR ;
else
    try
        L = chol(copula.Parameters(indNonConst,indNonConst));
        copula.cholR = L;
    catch
        error('Error: The copula correlation matrix is not positive definite or incorrectly defined!')
    end
end

% Produce U (of non - constant marginals)
U = zeros(size(Xin) ) ;
U(:,indNonConst) = V / L;
% the U's of constant marginals (if any) are zero