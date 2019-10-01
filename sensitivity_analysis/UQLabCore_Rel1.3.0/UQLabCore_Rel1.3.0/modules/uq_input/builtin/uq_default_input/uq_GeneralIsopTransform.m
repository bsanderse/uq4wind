function Y = uq_GeneralIsopTransform(X, X_Marginals, X_Copula, Y_Marginals, Y_Copula)
% Y = UQ_GENERALISOPTRANSFORM(X, X_Marginals, X_Copula, Y_Marginals, Y_Copula):
%     Maps a set X of samples from a random vector with an arbitrary 
%     distribution to a sample Y in another probability space. Both spaces  
%     need to be specified via the Marginals and Copula formalism, as 
%     described in the UQLab Input module user manual.
%
%     Note: non-constant marginal distributions can be transformed into
%     constant ones, but the converse is not true and will raise an error.
%     Constants can be mapped to different constants though.
%
% See also:  UQ_INVNATAFTRANSFORM, UQ_ROSENBLATTTRANSFORM, UQ_ISOPTRANSFORM

if any(isnan(X(:)))
    error('Requested Generalized probabilistic transformation for array X containing nans.')
end

[n, M] = size(X);

idNonConst_X = uq_find_nonconstant_marginals(X_Marginals);
idNonConst_Y = uq_find_nonconstant_marginals(Y_Marginals);
idConst_X = setdiff(1:M, idNonConst_X);
idConst_Y = setdiff(1:M, idNonConst_Y);

% Raise error if constants should be mapped to non-constants
idConst_XnotY = setdiff(idConst_X, idConst_Y);
if ~isempty(idConst_XnotY)
    msg = 'Constant marginal X_i cannot be mapped to non-constant marginal Y_i';
    error('%s, i=%s', msg, mat2str(idConst_XnotY))
end

% Identify XX, the sub-vectors of X containing the non-constant variables only, 
% and define the non-constant marginals of X and Y and their copula
idNonConst = idNonConst_Y; 
idConst = idConst_Y; 

% Initialize Y
Y = nan(n, M);

% Transform constant marginals to new constants
for ii = idConst
    Y(:,ii) = Y_Marginals(ii).Parameters;
end

% If all variables are constant, coupled by a copula other than the 
% independence copula, raise error
if length(idConst) == M
    msg='Requested isoprobabilistic transformation between constant variables';
    if ~uq_isIndependenceCopula(X_Copula)
        error('%s X coupled by non-independence copula', msg)
    end
    if ~uq_isIndependenceCopula(Y_Copula)
        error('%s Y coupled by non-independence copula', msg)
    end
end

% Assign X_Copula and Y_Copula the field .Variables, if not existing yet
% and if they have length 1
if length(X_Copula) == 1 && ~uq_isnonemptyfield(X_Copula, 'Variables')
    X_Copula.Variables = 1:M;
end
if length(Y_Copula) == 1 && ~uq_isnonemptyfield(Y_Copula, 'Variables')
    Y_Copula.Variables = 1:M;
end

% Initialize Z, the Rosenblatt transform of X, to nans
Z = nan(n, M);

% Fill the columns of Z associated with non-constant variables with correct 
% values obtained by Rosenblatt transformation (the other columns will
% retain value -1)
for cc = 1:length(X_Copula)
    Cop = X_Copula(cc);
    Vars = Cop.Variables;
    % If some variables couples by Cop are constant, take them out; also
    % raise errors if they were coupled by non-independence copula
    if any(ismember(Vars, idConst))
        if uq_isIndependenceCopula(Cop)
            Vars = Vars(ismember(Vars, idNonConst)); % take non-constant only
            Cop = uq_IndepCopula(length(Vars)); % redefine Cop in lower dim
            Cop.Variables = 1:length(Vars);
        else
            error('X variables %s are constant but coupled by %s copula', ...
                mat2str(Vars), Cop.Type)
        end
    end
    
    if ~isempty(Vars)
        Z(:, Vars) = uq_RosenblattTransform(X(:, Vars), X_Marginals(Vars), Cop);
    end
end


% Perform inverse Roseblatt transform of the non-constant columns in Z
for cc = 1:length(Y_Copula)
    Cop = Y_Copula(cc);
    Vars = Cop.Variables;
    % If some variables coupled by Cop are constant, take them out; also
    % raise errors if they were coupled by non-independence copula
    if any(ismember(Vars, idConst))
        if uq_isIndependenceCopula(Cop)
            Vars = Vars(ismember(Vars, idNonConst));
            Cop = uq_IndepCopula(length(Vars)); % redefine Cop in lower dim
            Cop.Variables = 1:length(Vars);
        else
            error('Constant variables X_%s cannot be coupled by %s copula', ...
                mat2str(Vars), Cop.Type)
        end
    end
    
    if ~isempty(Vars)
        Y(:, Vars) = uq_invRosenblattTransform(Z(:, Vars), Y_Marginals(Vars), Cop);
    end
end

% Final check that all columns of Y have been dealt with (raise error
% otherwise)
if any(isnan(Y(:)))
    error('Generalized Isoprobabilistic Transform Y of X contains nans')
end
