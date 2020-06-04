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
BayesOpts.Input = myPrior;

% display input properties
uq_print(myPrior);
uq_display(myPrior);

%% set forward model
myForwardModel = uq_createModel(Model);    

%% Surrogate
if (Bayes_full == 0) % create a PCE surrogate model to be used
    if (Surrogate_model_type == 0)
        disp(['loading surrogate model from file: ' Surrogate_model_filename]);        
        mySurrogateModel = load(Surrogate_model_filename);

    elseif (Surrogate_model_type == 1)
        disp('creating surrogate model');        
        % use prior also as input uncertainties
        MetaOpts.Input     = myPrior;
        MetaOpts.FullModel = myForwardModel;   
        mySurrogateModel   = uq_createModel(MetaOpts);
    end
    % |mySurrogateModel| in lieu of the original |myForwardModel|:   
    BayesOpts.ForwardModel.Model = mySurrogateModel;

else % do Bayes with full model
    BayesOpts.ForwardModel.Model = myForwardModel;
end


%% Bayesian analysis options

% Run the Bayesian inversion analysis
BayesianAnalysis = uq_createAnalysis(BayesOpts);
% Print out a report of the results:
uq_print(BayesianAnalysis)
uq_display(BayesianAnalysis)
uq_display(BayesianAnalysis, 'meanConvergence', 'all')
uq_display(BayesianAnalysis, 'trace', 'all')
uq_display(BayesianAnalysis, 'acceptance', 'true')
uq_postProcessInversion(BayesianAnalysis,'pointEstimate', 'MAP')
uq_postProcessInversion(BayesianAnalysis,'gelmanRubin', 'true')
R_hat_full = BayesianAnalysis.Results.PostProc.MPSRF;

if R_hat_full <= 1.5
    disp('The MCMC simulation has converged')
else
    disp('The MCMC simulation has not converged. Increase the number of samples or fine tune the algorithm.')
end


