function F = uq_Kriging_eval_F( X, current_model )
% F = UQ_KRIGING_EVAL_F(X, KRGMODEL): calculates the trend matrix (a.k.a. 
% information matrix) F of the trend of a Kriging metamodel KRGMODEL on the 
% points specified in X
%
% See also UQ_KRIGING_CALCULATE, UQ_KRIGING_EVAL, UQ_KRIGING_INITIALIZE, 
% UQ_KRIGING_INITIALIZETREND

M = current_model.Internal.Runtime.MnonConst;
nonConstIdx = current_model.Internal.Runtime.nonConstIdx;
N = size(X,1);

% Obtain the current output
current_output = current_model.Internal.Runtime.current_output ;

%% IF the F function is given use it, otherwise proceed with polynomial F
if isfield(current_model.Internal.Kriging(current_output).Trend,'CustomF') &&...
        ~isempty(current_model.Internal.Kriging(current_output).Trend.CustomF)
    % F is either numeric or cell array of handles
    if isnumeric(current_model.Internal.Kriging(current_output).Trend.CustomF)
        % F is numeric
        F = repmat(current_model.Internal.Kriging(current_output).Trend.CustomF,N,1) ;
    else
        % This is the general case of custom trend
        P = length(current_model.Internal.Kriging(current_output).Trend.CustomF);
        F = zeros(N,P) ;
        for ii = 1 : P
            f_ii =  current_model.Internal.Kriging(current_output).Trend.CustomF{ii};
            F(:,ii) = f_ii(X);
        end
    end
else % the trend is not Custom
    
    PolyType = {current_model.Internal.Kriging(current_output).Trend.PolyTypes{nonConstIdx}};
    P = current_model.Internal.Kriging(current_output).Trend.Degree;
    
    %% first get the polynomial basis index matrix
    
    % current maximum degree. It is set to the first value of Degree array, even for a non
    % basis adaptive scheme
    MaxDegree = current_model.Internal.Kriging(current_output).Trend.Degree(1) ;
    
    % if truncation is set to manual use the user-supplied polynomial indices, otherwise
    % calculate them
    
    if isfield(current_model.Internal.Kriging(current_output).Trend.TruncOptions, 'Custom')% user-specified truncation
        current_model.Internal.Kriging(current_output).Trend.Indices = current_model.Internal.Kriging(current_output).Trend.TruncOptions.Custom ;
    else
        current_model.Internal.Kriging(current_output).Trend.Indices = uq_generate_basis_Apmj(0 : MaxDegree, M, ...
            current_model.Internal.Kriging(current_output).Trend.TruncOptions);
    end
    

    %% ok, now on to calculating the univariate polynomials
    univ_p_val = zeros(N,M, P+1);
    % evaluating the univariate polynomials in the experimental design
    for ii = 1 : M
        currPolyType = PolyType{ii};
        if iscell(currPolyType)
            currPolyType = cell2string(currPolyType);
        end
        switch lower(currPolyType)
            case 'simple_poly'
                univ_p_val(:,ii,:) = uq_eval_simple_poly (P, X(:,ii));
        end
    end
    % get F
    F = uq_PCE_create_Psi(current_model.Internal.Kriging(current_output).Trend.Indices, univ_p_val);
end