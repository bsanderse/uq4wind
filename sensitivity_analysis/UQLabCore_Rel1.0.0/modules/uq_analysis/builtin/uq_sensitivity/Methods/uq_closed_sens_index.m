function [Sclo,Cost,Options] = uq_closed_sens_index(Options,VariableSet)
% [SCLO,COST,OPTIONS] = UQ_CLOSED_SENS_INDEX(OPTIONS,VARIABLESET)
% produces the closed index of X_VARIABLESET, which corresponds to
% Kucherenko's first order effect, and the associated cost.
%   The OPTIONS are essentially the current analysis' options including a
%   field named IndexOpts that is of type structure and contains
%   .SampleSize, .Sampling, .Estimator and possibly samples .X1,.Y1,.X2,.Y2

% See also: UQ_KUCHERENKO_INDICES, UQ_SHAPLEY_INDICES, UQ_TOTAL_SENS_INDEX

%% GET THE INPUT ARGUMENTS

% Verbosity
Display = Options.Display;

% Input object
myInput = Options.Input;
% Just another check to make sure the copula is Gaussian or Independent
if  all(~strcmpi(myInput.Copula.Type,{'Independent','Gaussian'}))
    error('Cannot compute sensitivity indices for copula of type "%s"', myInput.Copula.Type)
end

M = length(myInput.Marginals);

% Get sampling specs
N = Options.IndexOpts.SampleSize;
Method = Options.IndexOpts.Sampling;

% Kucherenko estimator
Estimator = Options.IndexOpts.Estimator;

% Determine which variables are not in the VariableSet
noSet = myInput.nonConst(~ismember(myInput.nonConst,VariableSet));

% Get the sample location
Samples = Options.IndexOpts;

% Turn variable indices into logicals
% First check the format of the indices. We want them to be logical.
if islogical(VariableSet)
    % there should be M entries
    if length(VariableSet) ~= M
        fprintf('\n\nError: The conditioning indices are provided as logicals, \n but the length of the array is not equal to M!\n');
        error('While initiating the conditional sampling')
    end
elseif isnumeric(VariableSet)
    % maybe it's 1's and 0's and meant to be logical
    if all(ismember(VariableSet,[0 1])) && length(VariableSet)==M
        VariableSet = logical(VariableSet);
    
    % but maybe it's variable indices $\subset (1,...,M)$, turn into logical
    elseif all(VariableSet < M+1) && length(unique(VariableSet)) == length(VariableSet)
        logidx = false(1,M);
        logidx(VariableSet) = true;
        VariableSet = logidx;
        
    else
        fprintf('\n\nError: The provided conditioning indices are neither logical nor numeric!\n');
        error('While initiating the conditional sampling')
    end
else
    fprintf('\n\nError: The provided conditioning indices are neither logical nor numeric!\n');
    error('While initiating the conditional sampling')
end


%% CONDITIONAL SAMPLING & SAMPLE ASSEMBLING

switch Estimator
    case {'standard','modified'}
        [MixedSample,Options.IndexOpts] = uq_getKucherenkoSamples(myInput,N,Method,VariableSet,Options.IndexOpts,Estimator);
                
    case 'alternative'
        [MixedSample,Options.IndexOpts] = uq_getKucherenkoSamples(myInput,N,Method,VariableSet,Options.IndexOpts,Estimator);
        
    otherwise
        fprintf('\nError: Estimator "%s" not available for Kucherenko indices. Please\n',Estimator)
        fprintf('choose one of the available: Standard, Modified or Alternative.\n\n')
        error('While calculating the indices');
end

%% EVALUATION AND ESTIMATION OF THE INDICES

% Initiate Cost & add everytime uq_evalModel gets used
Cost = 0;

switch Estimator
    case 'standard'
                
        % Needed evaluations
        if ~isfield(Options.IndexOpts,'y1')
            % Save it for later use
            Options.IndexOpts.y1 = uq_evalModel(Options.Model,Options.IndexOpts.x1);
            Cost = Cost+N;
        end
        EvalMix = uq_evalModel(Options.Model,MixedSample);
        Cost = Cost+N;
        
        % Variance of each output
        D = var(Options.IndexOpts.y1,1);
        
        % Index estimation
        term1 = Options.IndexOpts.y1.*EvalMix;
        term2 = sum(Options.IndexOpts.y1,1)/N;
        Sclo = ( (sum(term1,1)) /N - (term2).^2 ) ./D;
        
    case 'modified'
        
        % Needed evaluations
        if ~isfield(Options.IndexOpts,'y1')
            % Save it for later use
            Options.IndexOpts.y1 = uq_evalModel(Options.Model,Options.IndexOpts.x1);
            Cost = Cost+N;
        end
        % Needed evaluations
        if ~isfield(Options.IndexOpts,'y2')
            % Save it for later use
            Options.IndexOpts.y2 = uq_evalModel(Options.Model,Options.IndexOpts.x2);
            Cost = Cost+N;
        end
        EvalMix = uq_evalModel(Options.Model,MixedSample);
        Cost = Cost+N;
        
        % Variance of each output
        D = var(Options.IndexOpts.y1,1);        
        
        % Index estimation
        term = Options.IndexOpts.y1.*(EvalMix-Options.IndexOpts.y2);
        Sclo = sum(term,1)./(N*D);
        
    case 'alternative'
        % Needed evaluations
        if ~isfield(Options.IndexOpts,'y1')
            % Save it for later use
            Options.IndexOpts.y1 = uq_evalModel(Options.Model,Options.IndexOpts.x1);
            Cost = Cost+N;
        end
        % Needed evaluations
        if ~isfield(Options.IndexOpts,'y2')
            % Save it for later use
            Options.IndexOpts.y2 = uq_evalModel(Options.Model,Options.IndexOpts.x2);
            Cost = Cost+N;
        end
        EvalMix = uq_evalModel(Options.Model,MixedSample);
        Cost = Cost+N;
        
        % Variance of each output
        D = var(Options.IndexOpts.y1,1);        
        
        % Index estimation
        term = Options.IndexOpts.y2.*(EvalMix-Options.IndexOpts.y1);
        Sclo = sum(term,1)./(N*D);        
end

Options.IndexOpts.TotalVariance = D;

