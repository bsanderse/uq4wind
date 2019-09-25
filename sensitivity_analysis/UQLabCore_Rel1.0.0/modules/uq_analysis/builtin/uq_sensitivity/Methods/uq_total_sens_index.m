function [Stot,Cost,Options] = uq_total_sens_index(Options,VariableSet)
% [STOT,COST] = UQ_TOTAL_SENS_INDEX(OPTIONS,VARIABLESET)
% produces the total index of X_VARIABLESET, which corresponds to
% Kucherenko's total effect, and the associated cost.
%   The OPTIONS are essentially the current analysis' options including a
%   field named IndexOpts that is of type structure and contains
%   .SampleSize, .Sampling and possibly samples .X1, .Y1, .X2, .Y2
%
% See also: UQ_KUCHERENKO_INDICES, UQ_SHAPLEY_INDICES, UQ_CLOSED_SENS_INDEX

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

% Determine which variables are not in the VariableSet
noSet = myInput.nonConst(~ismember(myInput.nonConst,VariableSet));

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
% we need the standard/modified estimators cross-conditioned sample
estim = 'standard';
[MixedSample,Options.IndexOpts] = uq_getKucherenkoSamples(myInput,N,Method,~VariableSet,Options.IndexOpts,estim);

%% EVALUATION AND ESTIMATION OF THE INDICES

Cost = 0;

% Needed evaluations
if ~isfield(Options.IndexOpts,'y1')
    % Save it for later use
    Options.IndexOpts.y1 = uq_evalModel(Options.Model,Options.IndexOpts.x1);
    Cost = Cost+N;
end
EvalMix = uq_evalModel(Options.Model,MixedSample);
Cost = Cost+N;

% Variance of each output
Options.IndexOpts.TotalVariance = var(Options.IndexOpts.y1,1);

% Estimation of index
term = (Options.IndexOpts.y1 - EvalMix).^2;
Stot = sum(term,1)./(2*N*Options.IndexOpts.TotalVariance);

