%% RELIABILITY: TWO-DIMENSIONAL HAT FUNCTION
%
% This example showcases the application of various reliability analysis
% methods in UQLab to a two-dimensional hat function.

%% 1 - INITIALIZE UQLAB
%
% Clear all variables from the workspace,
% set the random number generator for reproducible results,
% and initialize the UQLab framework:
clearvars
rng(100,'twister')
uqlab

%% 2 - COMPUTATIONAL MODEL
%
% The two-dimensional hat function is defined as follows:
%
% $$g(x_1, x_2) = 20 - (x_1 - x_2)^2 - 8 (x_1 + x_2 - 4)^3$$
%
% Create a limit state function model based on the hat function
% using a string (vectorized):
ModelOpts.mString = '20 - (X(:,1)-X(:,2)).^2 - 8*(X(:,1)+X(:,2)-4).^3';
ModelOpts.isVectorized = true;
myModel = uq_createModel(ModelOpts);

%% 3 - PROBABILISTIC INPUT MODEL
%
% The probabilistic input model consists of two independent
% and identically-distributed Gaussian random variables:
%
% $X_i \sim \mathcal{N}(0.25, 1), \quad i = 1, 2$
%
% Specify the probabilistic model for the two input random variables:
InputOpts.Marginals(1).Name = 'X1'; 
InputOpts.Marginals(1).Type = 'Gaussian';
InputOpts.Marginals(1).Parameters = [0.25 1];

InputOpts.Marginals(2).Name = 'X2';  
InputOpts.Marginals(2).Type = 'Gaussian';
InputOpts.Marginals(2).Parameters = [0.25 1];

%%
% Create the INPUT object:
myInput = uq_createInput(InputOpts);

%% 4 - STRUCTURAL RELIABILITY
%
% Failure event is defined as $g(\mathbf{x}) \leq 0$.
% The failure probability is then defined as
% $P_f = P[g(\mathbf{x})\leq 0]$.
%
% Reliability analysis is performed with the following methods:
%
% * Monte Carlo simulation (MCS)
% * Subset simulation
% * Adaptive-Kriging-Monte-Carlo-Simulation (AK-MCS)
% * Adaptive-Polynomial-Chaos-Kriging-Monte-Carlo-Simulation (APCK-MCS)

%% 4.1 Monte Carlo simulation
%
% Select the Reliability module and the Monte Carlo simulation (MCS)
% method:
MCSOptions.Type = 'Reliability';
MCSOptions.Method = 'MCS';

%%
% Specify the maximum sample size, the size of the batch,
% and the target coefficient of variation:
MCSOptions.Simulation.MaxSampleSize = 1e6;
MCSOptions.Simulation.BatchSize = 1e5;
MCSOptions.Simulation.TargetCoV = 1e-2;

%%
% Run the Monte Carlo simulation:
myMCSAnalysis = uq_createAnalysis(MCSOptions);

%%
% Print out a report of the results:
uq_print(myMCSAnalysis)

%% 
% Visualize the results of the analysis:
uq_display(myMCSAnalysis)

%% 4.2 Subset simulation
%
% Select the Reliability module and the Subset simulation method:
SubsetSimOptions.Type = 'Reliability';
SubsetSimOptions.Method = 'Subset';

%%
% Specify the sample size in each subset:
SubsetSimOptions.Simulation.BatchSize = 1e4;

%%
% Run the subset simulation analysis
mySubsetSimAnalysis = uq_createAnalysis(SubsetSimOptions);

%%
% Print out a report of the results:
uq_print(mySubsetSimAnalysis)

%%
% Visualize the results of the analysis:
uq_display(mySubsetSimAnalysis)

%% 4.3 Adaptive-Kriging-Monte-Carlo-Simulation (AK-MCS)
%
% Select the Reliability module and the AK-MCS method:
AKOptions.Type = 'Reliability';
AKOptions.Method = 'AKMCS';

%%
% Specify the size of the Monte Carlo sample set used for 
% the Monte Carlo simulation
% and as the candidate set for the learning function:
AKOptions.Simulation.MaxSampleSize = 1e6;

%% 
% Specify the maximum number of sample points added
% to the experimental design:
AKOptions.AKMCS.MaxAddedED = 20;

%%
% Specify the initial experimental design:
AKOptions.AKMCS.IExpDesign.N = 20;
AKOptions.AKMCS.IExpDesign.Sampling = 'LHS';

%%
% Specify the options for the Kriging metamodel
% (note that all Kriging options are supported):
AKOptions.AKMCS.Kriging.Corr.Family = 'Gaussian';

%%
% Specify the convergence criterion
% for the adaptive experimental design algorithm
% (here, it is based on the failure probability estimate):
AKOptions.AKMCS.Convergence = 'stopPf';

%%
% Specify the learning function
% (here, it is the _expected feasibility function (EFF)_):
AKOptions.AKMCS.LearningFunction = 'EFF';

%%
% Run the AK-MCS analysis:
myAKAnalysis = uq_createAnalysis(AKOptions);

%%
% Print out a report of the results:
uq_print(myAKAnalysis)

%%
% Visualize the results of the analysis:
uq_display(myAKAnalysis)

%% 4.4 Adaptive-Polynomial-Chaos-Kriging-Monte-Carlo-Simulation (APCK-MCS)
%
% APCK-MCS is a variation of AK-MCS in which the Kriging model
% is replaced by a PC-Kriging (PCK) model.
%
% Select the Reliability module and the AK-MCS method:
APCKOptions.Type = 'Reliability';
APCKOptions.Method = 'AKMCS';

%%
% Select PCK as the metamodel type:
APCKOptions.AKMCS.MetaModel = 'PCK';

%%
% Specify the correlation function for the Kriging metamodel 
% (here, it is the Gaussian correlation function):
APCKOptions.AKMCS.PCK.Kriging.Corr.Family = 'Gaussian';

%% 
% Set the number of sample points in the initial experimental design:
APCKOptions.AKMCS.IExpDesign.N = 5;

%%
% Specify the size of the Monte Carlo set used for
% the Monte Carlo simulation
% and as the candidate set for the learning function:
APCKOptions.Simulation.MaxSampleSize = 1e6;

%%
% Run the APCK-MCS analysis:
myAPCKAnalysis = uq_createAnalysis(APCKOptions);

%%
% Print out a report of the results:
uq_print(myAPCKAnalysis)

%% 
% Visualize the results of the analysis:
uq_display(myAPCKAnalysis)