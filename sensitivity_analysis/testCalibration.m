%% Bayesian calibration: ECN AeroModule
% In this example, the ECN Aero-Module (aerodynamic model) is calibrated
% using DANAERO measurements of axial force along the blade radial locations.
% The calibration is solved using an MCMC sampler using the full forward
% model and using a polynomial chaos expansion (PCE) surrogate model.
clc
close all
clearvars
% fix random  number generator (this makes results reproducible, but is
% not always wanted):
rng default

root_folder = pwd;

% check current path for existence of folders from this working directory
path_all = strsplit(path,';');
ind = startsWith(path_all,root_folder);
% remove those
rmpath(strjoin(string(path_all(ind)),';'))

%% Case study
caseName = 'NM80_calibrate'; % 'airfoil_lift','NM80', etc;
% specify directory which contains test case settings and model
% often this is simply the caseName
input_file = caseName; 

%% Initialize UQlab
% add path
run('config.m');
addpath(genpath(UQLab_path));
% start uqlab
uqlab

%% Add paths for dependent routines located in the directories:'NURBS','AEROmoduleWrapper' and 'Geometry'
% remove folders from the path to prevent that files from wrong folders are called
addpath(fullfile(root_folder,'AEROmoduleWrapper'));
addpath(fullfile(root_folder,'NURBS'));
addpath(fullfile(root_folder,'Geometry'));
addpath(fullfile(root_folder,'cases',caseName));
addpath(fullfile(root_folder,'Other'));

%% Initialization
run(['cases/' input_file '/initialize_calibration.m']);


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
ndim = length(myPrior.Marginals);
% set unperturbed vector:
for i=1:ndim
    % we take the mean of each parameter as the unperturbed
    % condition
    X_unperturbed(1,i) = myPrior.Marginals(i).Moments(1);
end
if (exist('test_run','var'))
    if (test_run == 1)
        disp('Performing test run at unperturbed (mean value) settings');
        
        Y_unperturbed = uq_evalModel(X_unperturbed);
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

% This is a bug fix for UQLab
% UQLab gives an error if the number of field names in the prior are not
% the same as in the discrepancy options
% see https://uqworld.org/t/bayesian-inference-error-when-fields-of-prior-are-different-from-field-of-discrepancy-options/937

% Here we add the fields of the discrepancy to be the same as those of the prior
for i=1:ndim
    fieldnames_prior = fieldnames(myPrior.Marginals(i));
    fieldnames_disc  = fieldnames(DiscrepancyOpts(i).Prior.Marginals);
    match_fieldnames = isfield(DiscrepancyOpts(i).Prior.Marginals,fieldnames_prior);
    ind = find(match_fieldnames==0);
    for j=1:length(ind)
        add_field = fieldnames_prior(ind(j));
        % add field and set to be empty
        DiscrepancyOpts(i).Prior.Marginals.(add_field{1}) = [];
    end
end

%%
BayesianAnalysis = uq_createAnalysis(BayesOpts);

%% Post-processing
pp_file = ['cases/' input_file '/postProcessing_calibration.m'];
if (exist(pp_file,'file')==2)
    run(pp_file);
else
    warning('postprocessing file not available');
end

