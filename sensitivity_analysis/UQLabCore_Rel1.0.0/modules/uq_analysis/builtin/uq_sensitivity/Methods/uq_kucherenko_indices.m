function Results = uq_kucherenko_indices(current_analysis)
% RESULTS = UQ_KUCHERENKO_INDICES(ANALYSISOBJ) produces the first order and
% total Kucherenko importance indices. They form a generalisation of the
% Sobol' indices for dependent input variables.

% See also: UQ_SENSITIVITY,UQ_INITIALIZE_UQ_SENSITIVITY

%% GET THE OPTIONS
Options = current_analysis.Internal;

% Current model
CurrentModel = Options.Model;

% Verbosity
Display = Options.Display;

% Number of variables:
M = Options.M;
% Analysed variables
Mprime = sum(Options.FactorIndex);
MprimeIdx = find(Options.FactorIndex);

% Input object
myInput = Options.Input;
if  all(~strcmpi(myInput.Copula.Type,{'Independent','Gaussian'}))
    error('Cannot compute sensitivity indices for copula of type "%s"', myInput.Copula.Type)
end


%% GET THE INDICES
% Pre allocation of indices as cell-array to avoid "growing" as a result of
% multiple outputs
First = cell(Mprime,1);
Cost1 = zeros(Mprime,1);

Total = cell(Mprime,1);
CostTot = zeros(Mprime,1);

% Cycle through the input variables find(FactorIdx)
for ii = 1:Mprime
    % First order indices and associated cost plus Options with new samples
    [First{ii}, Cost1(ii), Options] = uq_closed_sens_index(Options,MprimeIdx(ii));
        
    % Total indices and associated cost
    [Total{ii}, CostTot(ii), ~] = uq_total_sens_index(Options,MprimeIdx(ii));
end


%% COLLECT THE INDICES IN THE RESULT STRUCTURE
First = vertcat(First{:});
Total = vertcat(Total{:});

Results.FirstOrder = zeros(M,size(First,2));
Results.FirstOrder(MprimeIdx,:) = First;
Results.Total = zeros(M,size(Total,2));
Results.Total(MprimeIdx,:) = Total;

Results.TotalVariance = Options.IndexOpts.TotalVariance;

% Cost
if ~isprop(CurrentModel, 'MetaType')
    Results.Cost = sum(Cost1) + sum(CostTot);
else
    Results.Cost = 0;
end

if Display > 0
    fprintf('\nKucherenko: finished.\n');
end