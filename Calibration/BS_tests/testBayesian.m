%% Sobol indices computed with Monte Carlo and PCE-type methods

clc
close all
clearvars


%% initialize UQlab

% add path
addpath(genpath('../../UQLabCore_Rel1.3.0/'));
% start uqlab
uqlab


%% mdoel

ModelOpts.mFile = 'uq_SimplySupportedBeam';
ModelOpts.isVectorized = true;

myForwardModel = uq_createModel(ModelOpts);

%% prior
PriorOpts.Marginals(1).Name = 'b';               % beam width
PriorOpts.Marginals(1).Type = 'Constant';
PriorOpts.Marginals(1).Parameters = [0.15];      % (m)

PriorOpts.Marginals(2).Name = 'h';               % beam height
PriorOpts.Marginals(2).Type = 'Constant';
PriorOpts.Marginals(2).Parameters = [0.3];       % (m)

PriorOpts.Marginals(3).Name = 'L';               % beam length
PriorOpts.Marginals(3).Type = 'Constant';
PriorOpts.Marginals(3).Parameters = 5;           % (m)

PriorOpts.Marginals(4).Name = 'E';               % Young's modulus
PriorOpts.Marginals(4).Type = 'LogNormal';
PriorOpts.Marginals(4).Moments = [30000 4500];   % (MPa)

PriorOpts.Marginals(5).Name = 'p';               % uniform load
PriorOpts.Marginals(5).Type = 'Gaussian';
PriorOpts.Marginals(5).Moments = [0.012 0.012*0.05]; % (kN/m)

myPriorDist = uq_createInput(PriorOpts);

% display input properties
uq_print(myPriorDist);
uq_display(myPriorDist);

%% data
myData.y = [12.84; 13.12; 12.13; 12.19; 12.67]/1000; % (m)
myData.Name = 'Mid-span deflection';

%% Bayes options and run
BayesOpts.Type = 'Inversion';
BayesOpts.Data = myData;
myBayesianAnalysis = uq_createAnalysis(BayesOpts);

%% postprocessing
uq_display(myBayesianAnalysis)
