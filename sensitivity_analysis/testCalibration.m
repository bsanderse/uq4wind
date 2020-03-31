%% Bayesian calibration with UQLab

clc
close all
clearvars

rng default

caseName = 'aero_module'; % 'airfoil_lift','aero_module', etc;
input_file = caseName; % specify directory which contains test case settings and model

%% Add paths for dependent routines located in the directories 'NURBS','AEROmoduleWrapper' and 'Geometry'
addpath([pwd,'/AEROmoduleWrapper/']);
addpath([pwd,'/NURBS/']);
addpath([pwd,'/Geometry/']);

%% initialize UQlab

% add path
run('config.m');
addpath(genpath(UQLab_path));

% start uqlab
uqlab


%% process input files
run(['cases/' input_file '/initialize_calibration.m']);

%% set prior distribution
myPrior = uq_createInput(Prior);

% display input properties
uq_print(myPrior);
uq_display(myPrior);

%% set forward model
ForwardModel = uq_createModel(Model);

%% start Bayesian calibration
BayesianAnalysis = uq_createAnalysis(BayesOpts);

%% postprocessing
uq_display(BayesianAnalysis)
