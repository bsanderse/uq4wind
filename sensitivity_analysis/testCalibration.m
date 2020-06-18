%% Bayesian calibration: ECN AeroModule
% In this example, the ECN Aero-Module (aerodynamic model) is calibrated
% using DANAERO measurements of axial force along the blade radial locations.
% The calibration is solved using an MCMC sampler using the full forward
% model and using a polynomial chaos expansion (PCE) surrogate model.
clc
close all
clearvars
rng default

%% Case study
caseName = 'aero_module'; % 'airfoil_lift','aero_module', etc;
input_file = caseName; % specify directory which contains test case settings and model

%% Add paths for dependent routines located in the directories:'NURBS','AEROmoduleWrapper' and 'Geometry'
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
% display input properties
uq_print(myPrior);
uq_display(myPrior);

%% Set forward model
% The Aero-Module model has been implemented in the function
% |uq_createModel(Model)| supplied with UQLab. The function
% evaluates the model using the input parameters |P| given 
% in the structure |getParameterAeroModule(turbineName)|.
myForwardModel = uq_createModel(Model);

%% Loading full/surrogate model for Bayesian analysis
if (Bayes_full == 0) % create a PCE surrogate model to be used
    if (Surrogate_model_type == 0)
        disp(['loading surrogate model from file: ' Surrogate_model_filename]);
        loaded_surrogate_model = load(Surrogate_model_filename);
        BayesOpts.ForwardModel.Model = loaded_surrogate_model.mySurrogateModel;
    elseif (Surrogate_model_type == 1)
        disp('creating surrogate model');
        % use prior also as input uncertainties
        MetaOpts.Input     = myPrior;
        MetaOpts.FullModel = myForwardModel;
        mySurrogateModel   = uq_createModel(MetaOpts);
        BayesOpts.ForwardModel.Model = mySurrogateModel;
    end
else % do Bayes with full model
    BayesOpts.ForwardModel.Model = myForwardModel;
end

%% Set Bayes options
BayesOpts.Prior = myPrior;                  % Prior
BayesOpts.Data  = Data;                     % Measurement data
BayesOpts.Type  = 'Inversion';              % Calibration
BayesOpts.Discrepancy = DiscrepancyOpts;    % Likelihood
BayesOpts.Solver = Solver;                  % MCMC

%% Run the Bayesian inversion analysis
BayesianAnalysis = uq_createAnalysis(BayesOpts);

%% Post-processing
run(['cases/' input_file '/PostProcessingCalibration.m']);

