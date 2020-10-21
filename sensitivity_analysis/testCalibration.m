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
caseName = 'NewMexico_calibrate'; % 'airfoil_lift','NM80', etc;
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
addpath([pwd,'/AEROmoduleWrapper/']);
addpath([pwd,'/NURBS/']);
addpath([pwd,'/Geometry/']);
%addpath(strcat(pwd,'/AEROmodule/',turbineName,'/output'));

%% empty the contents of the output folder of the AeroModule to prevent that old information is being loaded
delete(strcat(pwd,'/AEROmodule/',turbineName,'/current/output/*'));

%% Set prior distribution
myPrior = uq_createInput(Prior);
% display input properties
uq_print(myPrior);
uq_display(myPrior);
pause(0.01)

%% Set forward model
% The Aero-Module model has been implemented in the function
% |uq_createModel(Model)| supplied with UQLab. The function
% evaluates the model using the input parameters |P| given
% in the structure |getParameterAeroModule(turbineName)|.
myForwardModel = uq_createModel(Model);

% do a test run with the forward model at unperturbed settings
if (exist('test_run','var'))
    if (test_run == 1)
        disp('Performing test run at unperturbed settings');
        uq_evalModel(zeros(1,ndim));
    end
end

%% Loading full/surrogate model for Bayesian analysis
if (Bayes_full == 0) % create a PCE surrogate model to be used
    if (Surrogate_model_type == 0)
        disp(['loading surrogate model from file: ' Surrogate_model_filename]);
        loaded_surrogate_model = load(Surrogate_model_filename);
        % check whether loaded surrogate model is having same features as
        % the uncertainties given in the prior (those we are trying to calibrate)
        if (~isequaln(myPrior.Marginals,loaded_surrogate_model.mySurrogateModel.Options.Input.Marginals))
            error('Marginals specified for Prior do not correspond with marginals used in surrogate model');
        end
        BayesOpts.ForwardModel.Model = loaded_surrogate_model.mySurrogateModel;
        
    elseif (Surrogate_model_type == 1)
        disp('creating surrogate model');
        % use prior also as input uncertainties for surrogate model
        % other MetaOpts should have been set in the initialize_calibration
        % file
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
disp('performing Bayesian analysis');
BayesianAnalysis = uq_createAnalysis(BayesOpts);

%% Post-processing
pp_file = ['cases/' input_file '/postProcessing_calibration.m'];
if (exist(pp_file,'file')==2)
    run(pp_file);
else
    warning('postprocessing file not available');
end

