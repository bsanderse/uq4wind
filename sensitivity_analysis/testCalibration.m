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

%% Surrogate
MetaOpts.Input = myPrior;
MetaOpts.FullModel = ForwardModel;
mySurrogateModel = uq_createModel(MetaOpts);
% |mySurrogateModel| in lieu of the original |myForwardModel|:
BayesOpts.ForwardModel.Model = mySurrogateModel;

%% Bayesian analysis for surrogate model
% Run the Bayesian inversion analysis:
myBayesianAnalysis_surrogateModel = uq_createAnalysis(BayesOpts);

%% Post-processing
% Print out a report of the results:
uq_print(myBayesianAnalysis_surrogateModel)
uq_display(myBayesianAnalysis_surrogateModel)
uq_display(myBayesianAnalysis_surrogateModel, 'meanConvergence', 'all')
uq_display(myBayesianAnalysis_surrogateModel, 'trace', 'all')
uq_display(myBayesianAnalysis_surrogateModel, 'acceptance', 'true')
uq_postProcessInversion(myBayesianAnalysis_surrogateModel,'pointEstimate', 'MAP')
uq_postProcessInversion(myBayesianAnalysis_surrogateModel,'gelmanRubin', 'true')
R_hat = myBayesianAnalysis_surrogateModel.Results.PostProc.MPSRF;

if R_hat <= 1.5
    disp('The MCMC simulation has converged')
else
    disp('The MCMC simulation has not converged. Increase the number of samples or fine tune the algorithm.')
end


