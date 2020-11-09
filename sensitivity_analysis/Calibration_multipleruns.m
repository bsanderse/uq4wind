%% Bayesian calibration: ECN AeroModule
% In this example, the ECN Aero-Module (aerodynamic model) is calibrated
% using DANAERO measurements of axial force along the blade radial locations.
% The calibration is solved using an MCMC sampler using the full forward
% model and using a polynomial chaos expansion (PCE) surrogate model.
clc
close all
clearvars
rng default

root_folder = pwd;
    
%% Case study
caseName = 'NewMexico_calibrate_allruns'; % 'airfoil_lift','NM80', etc;
% specify directory which contains test case settings and model
% often this is simply the caseName
input_file = caseName; 

%% Initialize UQlab
% add path
run('config.m');
addpath(genpath(UQLab_path));
% start uqlab
uqlab

%% Initialization
run(['cases/' input_file '/initialize_calibration.m']);

%% Add paths for dependent routines located in the directories:'NURBS','AEROmoduleWrapper' and 'Geometry'
addpath([root_folder,'/AEROmoduleWrapper/']);
addpath([root_folder,'/NURBS/']);
addpath([root_folder,'/Geometry/']);
%addpath(strcat(pwd,'/AEROmodule/',turbineName,'/output'));

%% empty the contents of the output folder of the AeroModule to prevent that old information is being loaded
delete(strcat(root_folder,'/AEROmodule/',turbineName,'/current/output/*'));


%% display prior distribution
uq_print(myPrior);
uq_display(myPrior);
pause(0.01)


%% Set Bayes options

% Loading full/surrogate model for Bayesian analysis
if (Bayes_full == 0) % create a PCE surrogate model to be used
    if (Surrogate_model_type == 0)

        BayesOpts.ForwardModel = loaded_surrogate_models.mySurrogateModels;
        
    elseif (Surrogate_model_type == 1)
        disp('creating surrogate model');
        % use prior also as input uncertainties for surrogate model
        % other MetaOpts should have been set in the initialize_calibration
        % file
        BayesOpts.ForwardModel = mySurrogateModels;
    end
else % do Bayes with full model
    BayesOpts.ForwardModel = myForwardModels;
end

BayesOpts.Prior = myPrior;                  % Prior
BayesOpts.Data  = Data;                     % Measurement data
BayesOpts.Type  = 'Inversion';              % Calibration
BayesOpts.Discrepancy = DiscrepancyOpts;    % Likelihood
BayesOpts.Solver = Solver;                  % MCMC options


%% Run the Bayesian inversion analysis
disp('performing Bayesian analysis');
BayesianAnalysis = uq_createAnalysis(BayesOpts);

%% Post-processing
pp_file = ['cases/' input_file '/postProcessing_calibration.m'];
if (exist(pp_file,'file')==2)
    run(pp_file);
else
    warning('postprocessing file not available');
end

