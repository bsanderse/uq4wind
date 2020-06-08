%% Bayesian calibration with UQLab
clc
close all
clearvars
rng default

%% Case study
caseName = 'aero_module'; % 'airfoil_lift','aero_module', etc;
input_file = caseName; % specify directory which contains test case settings and model

%% Add paths for dependent routines located in the directories 'NURBS','AEROmoduleWrapper' and 'Geometry'
addpath([pwd,'/AEROmoduleWrapper/']);
addpath([pwd,'/NURBS/']);
addpath([pwd,'/Geometry/']);
addpath([pwd,'/AEROmodule/NM80_calibrate/output']);

%% Initialize UQlab
% add path
run('config.m');
addpath(genpath(UQLab_path));
% start uqlab
uqlab

%% Initialization
run(['cases/' input_file '/initialize_calibration.m']);

%% Set prior distribution
myPrior = uq_createInput(Prior);
BayesOpts.Input = myPrior;

% display input properties
uq_print(myPrior);
uq_display(myPrior);

%% set forward model
myForwardModel = uq_createModel(Model);

%% Surrogate
Bayes_full = 0;
Surrogate_model_type = 0; % 0 or 1

if (Bayes_full == 0) % create a PCE surrogate model to be used
    if (Surrogate_model_type == 0)
        load('surrogate.mat');
        disp('loading surrogate model from file: surrogate.mat');
        
        
    elseif (Surrogate_model_type == 1)
        disp('creating surrogate model');
        % use prior also as input uncertainties
        MetaOpts.Input     = myPrior;
        MetaOpts.FullModel = myForwardModel;
        mySurrogateModel   = uq_createModel(MetaOpts);
        save surrogate.mat mySurrogateModel % Saves the surrogate model
    end
    % |mySurrogateModel| in lieu of the original |myForwardModel|:
    BayesOpts.ForwardModel.Model = mySurrogateModel;
    
else % do Bayes with full model
    BayesOpts.ForwardModel.Model = myForwardModel;
end


%% Bayesian analysis options
run(['cases/' input_file '/PostProcessingCalibration.m']);

