function [X, u] = uq_default_input(current_input,varargin)
% UQ_DEFAULT_INPUT generates samples of random vectors. The elements of the random vector
% may be dependent. The dependency is specified via the copula formalism. At the current stage 
% only Gaussian copulas are supported. This function is called by uq_getSample 
% when the currently selected input object is of type uq_default_input.
%
% UQ_DEFAULT_INPUT(N) generates N samples from the random vector specified 
% by the currently selected Input object using the default sampling strategy
% 
% UQ_DEFAULT_INPUT(N, Options) allows to specify additional
% options (e.g. the sampling strategy). For more information please refer
% to the UQLab user manual: The Input module
% 
% X = UQ_DEFAULT_INPUT(...) returns an NxM matrix of the samples in the
% physical space
%
% [X, u] = UQ_DEFAULT_INPUT(...) additionally returns an NxM matrix of
% the samples in the uniform space
%



%% Parameters
LHSiterations_default = 5;

%% Parse sampling options
nargs = length(varargin);
if ~nargs
   error('Number of samples not specified. Usage: uq_getSample(N)');
end
N = varargin{1};

if nargs > 1 
    sampling = varargin{2};
end

if nargs > 1 && ~isempty(sampling)
    % Obtain sampling method if specified
    current_input.Sampling.Method = sampling;
else 
    % otherwise use the default one
    current_input.Sampling.Method = current_input.Sampling.DefaultMethod;
end
samplingOptions = current_input.Sampling;

%% Parse additional options
if nargs > 2
    p = inputParser;
    % 1) LHS number of iterations
    p.addOptional('iterations', LHSiterations_default, @isscalar);
   
    % Now parse the additional options
    p.parse(varargin{3:end});
   
    % Assign the parsed values
    samplingOptions.LHSiterations = p.Results.iterations ;
else
    % Assign the default values
    samplingOptions.LHSiterations = LHSiterations_default;
end

%% Calculate/Retrieve samples
if strcmpi(current_input.Sampling.Method,'Data')
    % Since the sampling method is 'data', [X,u] are going to be retrieved from an external file 
    % It is assumed that the input data file correctly contains them, i.e.
    % no additional checks are taking place after the values are loaded IF
    % loading succeeds
    try
        load(current_input.DataFile);
    catch
        error('Error: The data file %s could not be loaded!',current_input.DataFile)
    end
else
    % In this case the samples need to be generated
    
    %% Keep track of constant/non-constant variables
    % Get the indices of the marginals of type constant (if any)
    Types = {current_input.Marginals(:).Type};
    
    % Get the indices of non-constant marginals
    indConst = find(strcmpi(Types, 'constant'));
    indNonConst    = 1:length(Types);
    indNonConst(indConst) = [];
    M = length(indNonConst) ;
    
    %% Generate samples in uniform space (u)
    % At least one non-constant element should exist in the random vector
    % to proceed in generating samples
    if M
        u(:,indNonConst) = uq_sampleU(N, M,samplingOptions);
        if strcmpi(current_input.Sampling.Method, 'Sobol')
            if ~isfield(current_input.Sampling, 'SobolGen') || ...
                    isempty(current_input.Sampling.SobolGen)
                % create the Sobol-set generator if it doesn't exist
                current_input.Sampling.SobolGen = sobolset(M);
                current_input.Sampling.SobolGen.Skip = N+1;
            else
            % If generator exists, update Sobol index seed
            current_input.Sampling.SobolGen.Skip = ...
                current_input.Sampling.SobolGen.Skip + N;
            end
            
        elseif strcmpi(current_input.Sampling.Method, 'Halton')
            
            if ~isfield(current_input.Sampling, 'HaltonGen') || ...
                    isempty(current_input.Sampling.HaltonGen)
                % create the Halton-set generator if it doesn't exist
                current_input.Sampling.HaltonGen = haltonset(M);
                % scramble the sequence
                current_input.Sampling.HaltonGen = ...
                    scramble( current_input.Sampling.HaltonGen,'RR2');
                current_input.Sampling.HaltonGen.Skip = N+1;
            else
                % If generator exists, update Halton index seed
                current_input.Sampling.HaltonGen.Skip = ...
                    current_input.Sampling.HaltonGen.Skip + N;
            end
        end
    end
    % Prepare the description of the u vector. This is going to be used in
    % order to perform an Isoprobabilistic transform of the samples from uniform
    % to physical space
    switch lower(current_input.Sampling.Method)
        case 'data'
            uMarginals = current_input.Marginals;
        otherwise
            % Assign marginal Types of u marginals (unit hypercube) that
            % correspond to non-constant marginals
            [uMarginals(indNonConst).Type] = deal('uniform') ;
            % Assign marginal Parameters of u marginals that
            % correspond to non-constant marginals
            [uMarginals(indNonConst).Parameters] = deal([0 1]) ;
            % Assign marginal Types of  u marginals that
            % correspond to constant marginals
            [uMarginals(indConst).Type] = deal('constant') ;
            % Assign marginal Parameters of u marginals that
            % correspond to constant marginals
            [uMarginals(indConst).Parameters] = deal(0.5) ;
            % Fix the values of all u samples that correspond to constant
            % marginal to 0.5 (mean of uniform distribution in [0,1])
            u(:,indConst) = 0.5 * ones(N, numel(indConst)) ;
    end
    
    %% Obtain X samples from u samples
    if strcmpi(current_input.Copula.Type, 'independent')
        %  If marginals are independent: directly perform the isop. transformation
        %  from u to X
        X = uq_IsopTransform(u, uMarginals, current_input.Marginals);
    elseif strcmpi(current_input.Copula.Type, 'gaussian')
        % Assign marginal Types of  U marginals (standard normal) that correspond
        % to non-constant marginals
        [UMarginals(indNonConst).Type] = deal('gaussian') ;
        % Assign marginal Parameters of  U marginals that correspond
        % to non-constant marginals
        [UMarginals(indNonConst).Parameters] = deal([0 1]) ;
        % Assign marginal Types of U marginals that correspond to constant
        % marginals
        [UMarginals(indConst).Type] = deal('constant') ;
        % Assign marginal Parameters of U marginals that correspond to constant
        % marginals
        [UMarginals(indConst).Parameters] = deal(0) ;
        %  If marginals are dependent, with elliptical copula:
        % First map: u -> U(standard normal space) using isoprobabilistic transformation
        U = uq_IsopTransform(u, uMarginals, UMarginals);
        % Then map: U -> X(physical space) using inverse Nataf transformation
        X = uq_invNatafTransform(U, current_input.Marginals, current_input.Copula);
    else
        %  If the dependency of the marginals is unknown: throw error
        error('Error: The input copula type is not known!')
    end
    
end




