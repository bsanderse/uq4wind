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
delete(strcat(pwd,'/AEROmodule/',turbineName,'/output/*'));

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
run(['cases/' input_file '/PostProcessingCalibration.m']);

%% Write calibration files
run(['cases/' input_file '/write_calibration.m']);
run(['AEROmodule/NM80_calibrate/cal_files.m']);

% %% Cross-validation (Optional)
% figure()
% R = [13,19,30,37];
% exp = [474.7 817.9 1210.8 1254.3];
% aero = [0.5378    0.8691    1.3804    1.4778]*10^3;
% calibrated = [0.4747    0.8179    1.2108    1.2543]*10^3;
% plot(R, exp,'g-*', 'LineWidth', 1)
% hold on
% plot(R, aero, 'r-o', 'LineWidth', 1)
% hold on
% plot(R, calibrated, 'b-o', 'LineWidth', 1)
% xlabel('R [m]')
% ylabel('\mu_{n} [N/m]')
% grid on
% legend('Experimental', 'Non-calibrated AeroModule', 'Calibrated AeroModule','Location','southeast')

