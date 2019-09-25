function success = uq_initialize_uq_default_input(module)
% UQ_INITIALIZE_UQ_DEFAULT_INPUT initializes an Input object of type
% uq_default_input
success = 0; 

%% parameters
DefaultCopula.Type = 'Independent';
defaultSampling.Method = 'MC';
defaultSampling.DefaultMethod = 'MC';
% number of points for the lookup tables (for kernel-smoothing)
ks.NP = 1e4;
% number of plus/minus standard deviations for the bounds of
% kernel-smoothing
ks.Nstd = 8;
%% retrieve the Input object
if exist('module', 'var')
    current_input = uq_getInput(module);
else
    current_input = uq_getInput;
end
% The module must be of "uq_default_input" type
if ~strcmp(current_input.Type, 'uq_default_input')
    error('uq_initialize_default_input error: you must initialize a uq_default_input type object, not a % one!!', current_input.Type);
end
% retrieve the options
Options = current_input.Options;

%% Initialize marginals
[Marginals, Options] = uq_process_option(Options, 'Marginals');
if Marginals.Missing|| Marginals.Invalid
    error('Input Marginals must be defined!');
end
marginals = Marginals.Value;
% Get the number of Inputs
M = length(marginals);

for ii = 1 : M
    % Fill in the variable names
    if ~isfield(marginals(ii),'Name') || isempty(marginals(ii).Name)
        % set a default name if none is specified 
        marginals(ii).Name = sprintf('X%d', ii);
    end
    %% Initialize kernel smoothing-based marginals: caching and lookup tables
    if strcmpi(marginals(ii).Type, 'ks')
        % the data to be used for KS
        KSPar = marginals(ii).Parameters;
        % kernel smoothing options
        ksopts = {};
        if isfield(marginals(ii), 'Options') && ~isempty(marginals(ii).Options)
            ff = fieldnames(marginals(ii).Options);
            for kk = 1:length(ff)
                ksopts = [ksopts ff{kk} marginals(ii).Options.(ff{kk})];
            end
        else
            marginals(ii).Options = [];
        end
        
        % set the bounds to reasonable values
        if isfield(marginals(ii), 'Bounds')
            uBounds = marginals(ii).Bounds;
        else
            stdKSPar = std(KSPar);
            uBounds(1) = min(KSPar) - ks.Nstd*stdKSPar;
            uBounds(2) = max(KSPar) + ks.Nstd*stdKSPar;
        end
        
        % create the u vector that will be used to cache the values of KS
        marginals(ii).KS.u = linspace(uBounds(1), uBounds(2), ks.NP);
        
        % cache the pdf
        marginals(ii).KS.pdf = ksdensity(KSPar, marginals(ii).KS.u, ksopts{:});
        [marginals(ii).KS.cdf, cdfidx] = unique(ksdensity(KSPar, marginals(ii).KS.u, ksopts{:}, 'function', 'cdf'));
        marginals(ii).KS.ucdf = marginals(ii).KS.u(cdfidx);
    end
end
% Fill-in Parameters (Moments) field if Moments (Parameters) is given
marginals = uq_MarginalFields(marginals);
% Get the indices of the marginals that are non-constant and store them in
% the input object
nonConstIdx = find(~ismember(lower({marginals.Type}),'constant'));
% Create the property that is needed
uq_addprop(current_input, 'nonConst', nonConstIdx);

%% Initialize copula
% If no copula is specified the copula type is assigned the 'independent' value 
[Copula, Options] = uq_process_option(Options, 'Copula', DefaultCopula);
copula = Copula.Value;
% Make sure that copula type is known
typeKnown = strcmpi(copula.Type, 'independent')||strcmpi(copula.Type, 'gaussian');
if ~typeKnown
    error('Error: Unknown copula type: %s !',copula.Type )
end

% Copula correlation matrix: consistency checks
% Naming convention: 
%   copula.Parameters : linear correlation matrix R(for Gaussian copula)
%   copula.RankCorr   : Spearman correlation matrix Rc (for Gaussian copula)

if ~strcmpi(copula.Type, 'independent')
    COPULA_R_EXISTS = isfield(copula,'Parameters') && ~isempty(copula.Parameters);
    COPULA_RC_EXISTS = isfield(copula,'RankCorr') && ~isempty(copula.RankCorr);
    if COPULA_RC_EXISTS && ~strcmpi(copula.Type, 'gaussian')
        error('The Spearman correlation matrix can only be used when the copula is Gaussian!')
    end
    switch (COPULA_R_EXISTS + COPULA_RC_EXISTS)
        case 1
            % do nothing
        case 2
            error('You can only define one of the possible copula correlation matrices!');
        case 0
            %No copula correlation matrix has been defined
            error('No copula correlation matrix defined!');
    end
    
    if COPULA_RC_EXISTS % The Spearman correlation matrix is given
        % the linear correlation matrix can be directly obtained
        copula.Parameters = 2*sin(pi/6*copula.RankCorr);
    end
    % R should be square
    if size(copula.Parameters,1) ~= size(copula.Parameters,2)
        error('Error: The copula matrix should be square!')
    end
    % R size should match the number of marginals
    if size(copula.Parameters,1) ~= length(marginals)
        error('Error: The copula matrix size is not consistent with the number of marginals!')
    end
    % R should be positive definite
    if nnz(eig(copula.Parameters(nonConstIdx,nonConstIdx)) > 0 ) < size(copula.Parameters(nonConstIdx,nonConstIdx),1)
        error('Error: The copula matrix is not positive definite!')
    end
else % independent copula
    % Create the copula's linear and Spearman correlation matrices and
    % store them
    copula.Parameters = eye(M);
    copula.RankCorr = copula.Parameters ;
end


%% Initialize sampling
% Set some defaults to the sampling method that is going to be used when
% getting samples from the random vector specified by this input object
[Sampling, Options] = uq_process_option(Options, 'Sampling', defaultSampling);
sampling = Sampling.Value;

%% Store the (initialized) main ingredients of the Input object
uq_addprop(current_input, 'Marginals', marginals);
uq_addprop(current_input, 'Copula', copula);
uq_addprop(current_input, 'Sampling', sampling);

success = 1;
