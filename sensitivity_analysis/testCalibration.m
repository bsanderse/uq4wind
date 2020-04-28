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
addpath([pwd,'/AEROmodule/NM80_calibrate/output']);
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
uq_print(BayesianAnalysis)
uq_display(BayesianAnalysis)
uq_display(BayesianAnalysis, 'meanConvergence', 'all')
uq_display(BayesianAnalysis, 'trace', 'all')
uq_display(BayesianAnalysis, 'acceptance', 'true')
uq_postProcessInversion(BayesianAnalysis,'pointEstimate', 'MAP')
uq_postProcessInversion(BayesianAnalysis,'gelmanRubin', 'true')
R_hat = BayesianAnalysis.Results.PostProc.MPSRF;

if R_hat <= 1.5
    disp('The MCMC simulation has converged')
else
    disp('The MCMC simulation has not converged. Increase the number of samples or fine tune the algorithm.')
end