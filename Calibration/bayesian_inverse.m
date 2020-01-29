% Calibration of a simply supported beam with multiple model inputs
% Quarter and midspan deflection
clc
close all
clearvars

%% Initialize UQLab 
uqlab
rng(100)

%% Prior distribution (Input parameters)
PriorOpts.Name = 'Model parameters prior';
PriorOpts.Marginals(1).Name = 'b';
PriorOpts.Marginals(1).Type = 'Constant';
PriorOpts.Marginals(1).Parameters = 0.15; % (m)
PriorOpts.Marginals(2).Name = 'h';
PriorOpts.Marginals(2).Type = 'Constant';
PriorOpts.Marginals(2).Parameters = 0.3; % (m)
PriorOpts.Marginals(3).Name = 'L';
PriorOpts.Marginals(3).Type = 'Constant';
PriorOpts.Marginals(3).Parameters = 5; % (m)
PriorOpts.Marginals(4).Name = 'E';
PriorOpts.Marginals(4).Type = 'Lognormal';
PriorOpts.Marginals(4).Moments = [30000 4500] ; % (MPa)
PriorOpts.Marginals(5).Name = 'p';
PriorOpts.Marginals(5).Type = 'Constant';
PriorOpts.Marginals(5).Parameters = 0.012; % (kN/m)

myPriorDist = uq_createInput(PriorOpts);


%% Forward model (Black-box model)
ModelOpts.Name = 'Forward model';
ModelOpts.mFile = 'uq_SimplySupportedBeamTwo';
myForwardModel = uq_createModel(ModelOpts);


%% Experimental data
V = [[ 8.98; 8.66; 8.85; 9.19; 8.64],... % L /4 (m)
[12.84; 13.12; 12.13; 12.19; 12.67]]/1000; % L /2 (m)
myData.Name = 'Beam quarter and midspan deflection'
myData.y = V;


%% Perform Bayesian inverse
BayesOpts.Type = 'Inversion';
BayesOpts.Data = myData;
myBayesianAnalysis = uq_createAnalysis(BayesOpts);

%% Discrepancy
%DiscrepancyOpts.Type = 'Gaussian';
%DiscrepancyOpts.Parameters = 1e-6; % single scalar
%BayesOpts.Discrepancy = DiscrepancyOpts;

%% Post-processing
uq_postProcessInversion(myBayesianAnalysis,'priorPredictive',1000);
uq_print(myBayesianAnalysis)
uq_display(myBayesianAnalysis)