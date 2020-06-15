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

% display input properties
uq_print(myPrior);
uq_display(myPrior);


%% set forward model
myForwardModel = uq_createModel(Model);

%% Surrogate
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
        %         s =  MetaOpts.ExpDesign.NSamples;
        %         %fname = sprintf('surrogate_%PCE.mat', s);
        %         save surrogate/surrogatePCE_2.mat mySurrogateModel
    end
    % |mySurrogateModel| in lieu of the original |myForwardModel|:
    %BayesOpts.ForwardModel.Model = loaded_surrogate_model.mySurrogateModel;
    %BayesOpts.ForwardModel.Model = mySurrogateModel;
else % do Bayes with full model
    BayesOpts.ForwardModel.Model = myForwardModel;
end

%% set Bayes options
BayesOpts.Prior = myPrior;           % Prior
BayesOpts.Data  = Data;              % measurement data
BayesOpts.Type  = 'Inversion';       
BayesOpts.Discrepancy = DiscrepancyOpts;    % likelihood
BayesOpts.Solver = Solver;           % MCMC

%% Run the Bayesian inversion analysis
BayesianAnalysis = uq_createAnalysis(BayesOpts);

%% postprocess Bayesian analysis options
run(['cases/' input_file '/PostProcessingCalibration.m']);

