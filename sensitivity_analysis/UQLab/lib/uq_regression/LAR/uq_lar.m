function results = uq_lar(Psi, Y, options)
% RESULTS = UQ_LAR(PSI, Y, OPTIONS): sparse regression with Least Angle
%   Regression. PSI is the matrix of evaluations of the regressors in the
%   experimental design points, while Y are the corresponding model
%   responses. The OPTIONS structure is optional and can contain any of the
%   following fields: 
%   
%    'normalize': 0 or 1 (default 1), normalize the both PSI and Y to have 0
%       mean and variance 1. Should in general be enabled
%    'early_stop': 0 or 1 (default 1), stop the LAR iterations if the accuracy starts
%       decreasing
%    'hybrid_lars': 0 or 1 (default 1), perform an OLS regression at the end
%       of the basis selection to obtain more accurate results. Note that
%       if this option is set to 0, it will disable normalization
%
%   The RESULTS structure contains the results as a structure with the following fields:
%     'coefficients':     the array of coefficients
%     'best_basis_index': the index of the iteration LAR has converged to
%     'max_score':        the maximum score of the best iteration (1 - LOO_k)
%     'LOO':              the Leave One Out error estimate for the best iteration
%     'normEmpErr':       the estimated normalizedEmpiricalError 
%     'nz_idx':           the index of non-zero regressors (w.r.t. the original
%                         PSI matrix)
%     'a_scores':         the vector of scores for each iteration of LAR
%     'coeff_array':      the matrix of the coefficients for each iteration
%                         of LAR
%
%   Example: calculate the regression of a given design matrix PSI on the
%            model evaluations Y with normalized lars, and plot the
%            evolution of the scores of the LAR iterations:
%
%   lar_options.early_stop = 1;
%   lar_options.normalize = 1;
%   lar_results = uq_lar(Psi, Y, lar_options)
%
% See also: UQ_PCE_OLS_REGRESSION, UQ_PCE_LOO_ERROR

%% Initialization of the default options
normflag = 1;
early_stop = 1;
hybrid_lars = 1;
modified_loo = 1;
hybrid_loo = 1;
no_selection = 0;
generalized_ls = 0;
display = 0;

%% parsing the options vector (if any)
if exist('options', 'var')
    % normalize prior to lars
    if isfield(options, 'normalize')
        normflag = options.normalize;
    end
   
    % early stop option
    % is this basis adaptive?
    if isfield(options, 'early_stop')
        early_stop = options.early_stop;
    end

    
    % hybrid lars option
    if isfield(options, 'hybrid_lars')
        hybrid_lars = options.hybrid_lars;
    end
    
    % disable basis selection (mostly for debug purposes only)
    if isfield(options, 'no_selection')
        no_selection = options.no_selection;
    end
    
    if isfield(options, 'loo_modified')
        modified_loo = options.loo_modified;
    end
    
    if isfield(options, 'loo_hybrid')
        hybrid_loo = options.loo_hybrid;
    end
    
    % verbosity level
    if isfield(options, 'display')
        DisplayLevel = options.display;
    else
        DisplayLevel = 0;
    end
    % generalized least squares in the presence of a covariance matrix for
    % Y
    if isfield(options, 'CY')
        CY = options.CY;
        generalized_ls = 1;
    end
end

% OLS options
olsoptions.modified_loo = modified_loo;

% get the number of data points and polynomials
[N ,P] = size(Psi);

%% trivial case: only 1 regressor, let's directly return the OLS solution
if P == 1
    % ols solution
    if generalized_ls
        olsoptions.CY = CY;
        ols_results = uq_PCE_OLS_regression(Psi, Y, olsoptions);
    else
        ols_results = uq_PCE_OLS_regression(Psi, Y, olsoptions);
    end
    
    % and now on to assigning the outputs
    results.coefficients = ols_results.coefficients;
    results.LOO = ols_results.LOO;
    results.normEmpErr = ols_results.normEmpErr;
    results.optErrorParams = ols_results.optErrorParams;
    results.coeff_array = results.coefficients;
    results.max_score   = 1-results.LOO;
    results.a_scores    = results.max_score;
    results.loo_scores  = results.LOO;
       
    results.best_basis_index = 1;
    results.nz_idx  = 1;
        
    % and the LARs estimate of the LOO
    results.LOO_lars = results.LOO;
    results.lars_idx = 1;
    
    
    % and exit from the algorithm
    return;
end

% first check for constant regressors (this is necessary before normalization)
constidx = ~any(diff(Psi, 1));
constval = Psi(1,constidx);

% normalize the regressors and the experimental design if necessary
if normflag
    % center and normalize the regressors. Note: this will remove the
    % constant regressors, as they are not used in LAR
    [Psi, mu_Psi, sigma_Psi] = zscore(Psi);
        
    % normalize the data as well
    mu_Y = mean(Y);
    Y = Y - mu_Y;
    Yvar = var(Y);
    modi_diag = 1/N*ones(N,1);
else
    modi_diag = zeros(N,1);
end

if generalized_ls
    % decorrelate the outputs
    CYinv = CY\eye(size(CY));
    L = chol(CYinv);
    Psi = L*Psi;
    Y = L*Y;
end

%% Initialization of the LAR iterations

nvars = min(N-2,P); % maximum number of active predictors: either the full set of basis elements, or N-1 (the - 2 is due to how k is incremented in the loop)
mu = zeros(size(Y)); % initial direction of the LAR 
a_coeff = []; % set of active coefficients
i_coeff = 1:P; % maximal set of coefficients
M = []; % initial information matrix (PsiT_j Psi_j)
maxk = 8*nvars; % maximum number of LARs iterations
coeff_array = zeros(nvars+1, P);

a_scores =-inf(1, nvars+1);
loo_scores = inf*ones(1, nvars+1);


% initialize the best lars leave one out score
refscore = inf;


%% iterative LAR
k = 0;
if DisplayLevel > 1
    fprintf('Maximum LARS candidate basis size: %d\n', P);
end


% get the scores for the constant term only
% initial score: just do OLS with the constant term, if it exists
ols_results = uq_PCE_OLS_regression(ones(size(Y,1),1), Y, olsoptions);
loo_scores(1) = ols_results.LOO;
a_scores(1) = 1-loo_scores(1);

maxiter = min(maxk,nvars);

while k < maxiter && Yvar
    k = k + 1;
    if DisplayLevel > 3
        fprintf('Computing LAR iteration %d\n',k);
    end
    
    % correlation with the residual
    cj = Psi'*(Y - mu);
    % getting the most correlated with the current inactive set
    [C, idx] = max(abs(cj(i_coeff)));
    idx = i_coeff(idx); % translate the index to an absolute index
    
    
    % invert the information matrix at the first iteration, later only update its value on
    % the basis of the previous inverted one
    if k == 1
        M = pinv(Psi(:,idx)'*Psi(:,idx));
    else
        x = (Psi(:,a_coeff))'*Psi(:,idx);
        r = (Psi(:,idx))'*Psi(:,idx) ;
        % update the information matrix based on the last calculated
        % inverse
        M = uq_blockwise_inverse(M,x,x',r) ; 
        
        % if the resulting matrix is singular, throw out a warning and use pinv
        if any(~isfinite(M(:)))
            try
                M = pinv(Psi(:,[a_coeff idx])'*Psi(:,[a_coeff idx]));
            catch me
                if DisplayLevel > 1
                    warning('singular design matrix. Skipping the current basis element and removing it from the candidates to improve stability');
                end
                i_coeff(i_coeff == idx) = [];
                continue;
            end
            
        end
    end
    
   % update the set of active predictors with the newfound idx
    a_coeff = [a_coeff idx];
    
    
    % now we get the vector of correlation signs (cf. pg 207 of Blatman's Thesis)
    s = sign(cj(a_coeff));
    % set the null signs to positive
    s(~s) = 1;
    
    %% now trying to calculate gamma, w and u based on Blatman's Thesis and Efron et al. 2004
    % full reference Blatman Thesis: Blatman, G, 2009, Adaptive sparse polynomial chaos
    % expansions for uncertainty propagation and sensitivity analysis, PhD Thesis,
    % Universit?? Blaise Pascal - Clermont II
    
    % full reference Efron et al. 2004 (gamma calculation): EFRON, HASTIE,JOHNSTONE and
    % TIBSHIRANI, Least Angle Regression, The Annals of Statistics 2004, Vol. 32, No. 2,
    % 407???499 

    % Variable naming after Blatman et al. 2004 (Simpler)
    c = 1/sqrt((s'*M)*s);
    
    % descent direction
    w = c*M*s;
    % descent versor in the residual space
    u = Psi(:,a_coeff)*w;
    % 
    aj = Psi'*u;
    
    % calculating gamma, based on Efron et al. 2004. 
    % Please note the variable naming change (Efron => this code):
    % AA => C, cj => cj(i_coeff). The rest is the same
    if k < nvars
        tmp = [(C - cj(i_coeff))./(c - aj(i_coeff)); (C + cj(i_coeff))./(c + aj(i_coeff))];
        % careful, this array may have all zeros!!!!
        gamma = min(tmp(tmp > 0));
        if isempty(gamma)
            gamma = 0;
            warning('Warning: numerical instability!! Gamma for LAR iteration %d was set to 0 to prevent crashes.', k);
        end
    else % we are at the OLS solution, so update the vectors accordingly
        gamma = C/c;
    end
    
    % remove the coefficient from the candidate predictors
    i_coeff(i_coeff == idx) = [];
    
    % now update the residual 
    mu = mu + gamma*u;
    
    % and finally update the coefficients in the correct direction to yield equicorrelated
    % vectors
    coeff_array(k+1,a_coeff) = coeff_array(k,a_coeff) + gamma*w' ;
    
    
    %% adaptive LARS: Modified Leave One Out error estimate Q^2:
    % based on Blatman, 2009 (PhD Thesis), pg. 115-116
        
    % modified Leave One Out Estimate should be (eq. 5.10, pg 116)
    % Err_loo = mean((M(x^(i))- metaM(x^(i))/(1-h_i)).^2), with
    % h_i = tr(Psi*(PsiTPsi)^-1 *PsiT)_i and
    % divided by var(Y) (eq. 5.11) and multiplied by the T coefficient (eq 5.13):
    % T(P,NCoeff) = NCoeff/(NCoeff-P) * (1 + tr(PsiTPsi)^-1)
    
    % corrected leave-one-out error:
    if ~hybrid_loo
        loo = uq_PCE_loo_error(Psi(:,a_coeff), M, Y, coeff_array(k+1,a_coeff)', modified_loo, modi_diag);
    else
        loo = uq_PCE_loo_error(Psi(:,a_coeff), M, Y, [], modified_loo, modi_diag);
    end
    if loo < 0
        warning('leave one out error negative!!')
    end
        
    loo_scores(k+1) = loo;
    a_scores(k+1) = 1 - loo;
    
        
    %  Stop the iterations if the error increases again:
    mm = round(nvars*0.1) ;
    mm = max(mm,100);
    mm = min(mm, nvars);
    
    % update the current best score
    if loo < refscore
        refscore = loo;
    end
    
    if k > mm
        % simply stop if the loo error is consistently above the reference loo for
        % at least 10% of the iterations
        if (loo_scores(k-mm) <= refscore) && early_stop
            if DisplayLevel > 1
                fprintf('Early stop at coefficient %d/%d \n', k-mm, P);
            end
            break;
        end
    end
end


% get the best score in the current array
if ~no_selection
    [maxScore, k] = max(a_scores);
else
    maxScore = 1-loo;
    k = length(a_scores);
end


%% Assigning the coefficients with the best candidate basis via OLS (Hybrid LARS)
% recompute the coefficients with the correct basis, by rescaling back the
% Psi and Y matrices

nz_idx = abs(coeff_array(k,:)) > 0;




% let's first scale back the Psi matrix to the original shape:
if normflag
    totcoeff = length(sigma_Psi);
    Psi = bsxfun(@plus, Psi * spdiags(sigma_Psi', 0, totcoeff, totcoeff),mu_Psi);
    
    % set back the constant column of Psi to the original value for the
    % first constant element only
    if sum(constidx)
        Psi(:,constidx(1)) = constval(1);
        
        % don't forget to calculate the first element and to add it to the accepted
        % ones!
        nz_idx(constidx(1)) = 1;
    end
    % and the model evaluations as well:
    Y = Y + mu_Y;
end

% now we assign the coefficients, either through an extra hybrid_lars
% iteration, or directly from OLS

coefficients = zeros(1,P);
if hybrid_lars
    % and now let's recalculate the coefficients via standard least squares
    ols_results = uq_PCE_OLS_regression(Psi(:,nz_idx), Y,olsoptions); % do not use the covariance option, as Y and Psi are already decorrelated
    coefficients(nz_idx) = ols_results.coefficients;
    results.LOO = ols_results.LOO;
    results.normEmpErr = ols_results.normEmpErr;
    results.optErrorParams = ols_results.optErrorParams;
else
    coefficients(nz_idx) = coeff_array(k,nz_idx);
    coefficients = coefficients';
    % if we normalized, we have to rescale now
    if normflag 
        coefficients(sigma_Psi ~= 0) = coefficients(sigma_Psi ~= 0)./sigma_Psi(sigma_Psi ~= 0)';
        % and add the constant term, if it exists
        if sum(constidx)
            coefficients(constidx(1)) =  mean(Y-Psi*coefficients);
        end
    end
    [results.LOO, results.normEmpErr, results.optErrorParams] = uq_PCE_loo_error(Psi(:,nz_idx), pinv(Psi(:,nz_idx).'*Psi(:,nz_idx)), Y, coefficients(nz_idx), 1);
end


if DisplayLevel > 1
    fprintf('LAR basis size: %d/%d\n', sum(nz_idx), P);
end


%% Assign the remaining outputs
results.coeff_array = coeff_array;
results.max_score   = maxScore;

% now let's check that the coefficients array has the correct dimensions, 
% otherwise clear it (it means it is outdated)

results.coefficients = coefficients;
results.a_scores     = a_scores;
results.loo_scores     = loo_scores;

results.best_basis_index = k;
results.nz_idx  = nz_idx;

% useful to reorder the matrix
results.lars_idx = [constidx(1) a_coeff(1:(k-1))];

% now get the actual error from the results of the hybrid LARS
% and the LARs estimate of the LOO
results.LOO_lars = 1 - maxScore;
