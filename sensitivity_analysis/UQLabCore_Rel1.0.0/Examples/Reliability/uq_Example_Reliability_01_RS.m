%% RELIABILITY: RESISTANCE VS. STRESS (R-S)
%
% This example showcases the application of different reliability analysis
% methods available in UQLab to the simple R-S example.

%% 1 - INITIALIZE UQLAB
%
% Clear all variables from the workspace, set the random number generator
% for reproducible results, and initialize the UQLab framework:
clearvars
rng(100,'twister')
uqlab

%% 2 - COMPUTATIONAL MODEL
%
% The R-S function is defined as:
%
% $$g(\mathbf{x}) = R - S$$
%
% where $\mathbf{x} = \{R,S\}$.
% $R$ and $S$ are the resistance and stress variables, respectively.
%
% Create a limit state function model using a string (vectorized):
ModelOpts.mString = 'X(:,1) - X(:,2)';
ModelOpts.isVectorized = true; 

myModel = uq_createModel(ModelOpts);

%% 3 - PROBABILISTIC INPUT MODEL
%
% The probabilistic input model consists of two independent Gaussian 
% random variables:
%
% $$R \sim \mathcal{N}(5, 0.8), \; S \sim \mathcal{N}(2, 0.6)$$
%
% Specify the probabilistic input model for the $R$ and $S$ variables:
InputOpts.Marginals(1).Name = 'R'; % resistance variable
InputOpts.Marginals(1).Type = 'Gaussian';
InputOpts.Marginals(1).Moments = [5 0.8];

InputOpts.Marginals(2).Name = 'S'; % stress variable
InputOpts.Marginals(2).Type = 'Gaussian';
InputOpts.Marginals(2).Moments = [2 0.6];

%%
% Create the INPUT object:
myInput = uq_createInput(InputOpts);

%% 4 - RELIABILITY ANALYSIS
%
% Failure event is defined as $g(\mathbf{x}) \leq 0$.
% The failure probability is then defined
% as $P_f = P[g(\mathbf{x}) \leq 0]$.
%
% Reliability analysis is performed with the following methods:
%
%  * Monte Carlo simulation
%  * First-order reliability method (FORM)
%  * Importance sampling (IS)
%  * Subset simulation
%  * Adaptive Kriging-Monte Carlo Simulation (AK-MCS)

%% 4.1 Monte Carlo simulation
%
% Specify the Reliability module and select the Monte Carlo simulation
% (MCS) method:
MCSOptions.Type = 'Reliability';
MCSOptions.Method = 'MCS';

%%
% Specify the maximum sample size:
MCSOptions.Simulation.MaxSampleSize = 1e6;

%%
% Run reliability analysis with MCS:
MCSAnalysis = uq_createAnalysis(MCSOptions);

%%
% Print out a report of the results:
uq_print(MCSAnalysis)

%%
% Create a graphical representation of the results:
uq_display(MCSAnalysis)

%% 4.2 FORM
%
% Select FORM as the reliability method:
FORMOptions.Type = 'Reliability';
FORMOptions.Method = 'FORM';

%%
% Run the FORM analysis:
FORMAnalysis = uq_createAnalysis(FORMOptions);

%%
% Print out a report of the results:
uq_print(FORMAnalysis)

%%
% Create a graphical representation of the results:
uq_display(FORMAnalysis)

%% 4.3 Importance sampling
%
% Select Importance sampling (IS) as the reliability analysis tool:
ISOptions.Type = 'Reliability';
ISOptions.Method = 'IS';

%%
% Run the IS reliability analysis:
ISAnalysis = uq_createAnalysis(ISOptions);

%%
% Print out a report of the results:
uq_print(ISAnalysis)

%%
% Create a graphical representation of the results:
uq_display(ISAnalysis)

%% 4.4 Subset simulation
%
% Select subset simulation as the reliability analysis tool:
SSimOptions.Type = 'Reliability';
SSimOptions.Method = 'Subset';

%%
% Run the subset simulation:
SSimAnalysis = uq_createAnalysis(SSimOptions);

%%
% Print out a report of the results:
uq_print(SSimAnalysis)

%%
% Create a graphical representation of the results:
uq_display(SSimAnalysis)

%% 4.5 AK-MCS
%
% Select AK-MCS as the reliability analysis tool:
AKMCSOptions.Type = 'Reliability';
AKMCSOptions.Method = 'AKMCS';

%%
% Run the AK-MCS analysis
AKMCSAnalysis = uq_createAnalysis(AKMCSOptions);

%%
% Print out a report of the results:
uq_print(AKMCSAnalysis)

%%
% Create a graphical representation of the results:
uq_display(AKMCSAnalysis)