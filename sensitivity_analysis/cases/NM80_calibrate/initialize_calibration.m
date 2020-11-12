%% This file is used to initialize all settings in order to perform the calibration of the AeroModule 
% for the case that experimental data corresponding to different operating
% conditions is available

%% Turbine name and folders
% Name of Matlab file that describes the uncertain parameters used for calibration
% check this file for definition of uncertainties
turbineName = 'NM80_calibrate'; 

% folder used to get reference input files
ref_folder      = 'AEROmodule/NM80_calibrate/reference/';
% folder used to write new input files
current_folder  = 'AEROmodule/NM80_calibrate/current/';
% name of input file is assumed to be input.txt


%% Forward model description
% Name of Matlab file representing the model
Model.mHandle = @aero_module;
% Quantity of interest
% note that we look at the mean (in time) of the QoI, as the AeroModule predictions
% for this test case are basically steady state.
QoI = 'Sectional_normal_force';


%% Experimental data
% Marco's (ECN) script for reading the data in N/m
% Currently, two options can be passed on:
% 1. raw_normal.dat for 10 minute normal force measurements
% 2. raw_tangential.dat for 10 minute tangential force measurements
filename_exp = fullfile(root_folder,'..','Experimental','WINDTRUE','raw_normal.dat');
output_raw   = readNM80(filename_exp, 2);

% the position of the sections of the experimental data which are used for
% interpolation of the aeromodule results: see NM80_calibrate_readoutput.m
r_exp_data = [11.87, 17.82, 28.97, 35.53]; % Measurement radial stations

 % 'mean', 'full', 'synthetic'; other options also possible but need to be implemented
 % below
Data_type = 'full';

% Because the model has different discrepancy options at different radial locations,
% the measurement data is stored in four different data structures:
% these data are column vectors; each entry corresponding to a different
% measurement
Data(1).y = output_raw.Fy03; % [N/m]
Data(1).Name = 'Fy03';
Data(1).MOMap = 1; % Model Output Map 1

Data(2).y = output_raw.Fy05; % [N/m]
Data(2).Name = 'Fy05';
Data(2).MOMap = 2; % Model Output Map 2

Data(3).y = output_raw.Fy08; % [N/m]
Data(3).Name = 'Fy08';
Data(3).MOMap = 3; % Model Output Map 3

Data(4).y = output_raw.Fy10; % [N/m]
Data(4).Name = 'Fy10';
Data(4).MOMap = 4; % Model Output Map 4


switch Data_type
    
    case 'mean'
        
        % Because the model has different discrepancy options at different radial locations,
        % the measurement data is stored in four different data structures:
        Data(1).y = mean(Data(1).y); % [N/m]        
        Data(2).y = mean(Data(2).y); % [N/m] 
        Data(3).y = mean(Data(3).y); % [N/m]        
        Data(4).y = mean(Data(4).y); % [N/m]
        
    case 'full'
        Data(1).y = Data(1).y(1:100); % [N/m]        
        Data(2).y = Data(2).y(1:100); % [N/m] 
        Data(3).y = Data(3).y(1:100); % [N/m]        
        Data(4).y = Data(4).y(1:100); % [N/m]
        
        
    case 'synthetic'
        % Y0 is obtained by running the code at unperturbed conditions
        Y0 = 1.0e+03 *[0.5385    0.8705    1.3804    1.4787];
        N_synth = 500;
        sigma   = [10;10;10;10];
        
        for i=1:4
            Data(i).y = Y0(i) + sigma(i)*randn(N_synth,1);
        end
end


%% Get and store uncertain and fixed parameters

% Pass parameters to model via the cell array FixedInputs
[FixedParameters,UncertainInputs] = getParameterAeroModule(turbineName);

FixedParameters.root_folder    = root_folder;
FixedParameters.ref_folder     = ref_folder;
FixedParameters.current_folder = current_folder;
FixedParameters.QoI            = QoI;
FixedParameters.r_exp          = r_exp_data;

P.FixedParameters = FixedParameters;
P.UncertainInputs = UncertainInputs;

Model.Parameters = P;
%  note that the model output is (in general) a vector, containing e.g. the
%  force at different radial sections; this does not require any additional
%  specification in UQLab
%  the model is however not vectorized in the sense that we cannot give all
%  the possible parameter values at once as input to AeroModule, but
%  instead we need to do this sequentially
Model.isVectorized = false;

%% Prior
% Set the Prior equal to the Input
ndim  = length(UncertainInputs.Marginals);
Prior = UncertainInputs;


%% Likelihood description
% Here, the discrepancy for each data point y are chosen to be
% independent and identically distributed Gaussian random variables with
% mean zero and unknown variance. A prior is needed for the variance.
% we assume the following prior for sigma^2
% discrepancy bounds for prior (specify here the expected variance between model
% and data, and multiply by some factor, e.g. 5 to have sufficiently broad prior
prior_sigma2_bounds = 5*[0 100; 0 100; 0 100; 0 100];


for i=1:size(prior_sigma2_bounds,1)

    DiscrepancyPriorOpts.Name = strcat('Prior of Sigma2 ',num2str(i));
    DiscrepancyPriorOpts.Marginals(1).Name = strcat('Sigma2-',num2str(i));
    DiscrepancyPriorOpts.Marginals(1).Type = 'Uniform';
    DiscrepancyPriorOpts.Marginals(1).Parameters = prior_sigma2_bounds(i,:);
    DiscrepancyPrior = uq_createInput(DiscrepancyPriorOpts);

    DiscrepancyOpts(i).Type = 'Gaussian';
    DiscrepancyOpts(i).Prior = DiscrepancyPrior;
        
end


%% Surrogate model options
test_run = 0; % perform test run with Forward Model without uncertainties

% Switch for Bayesian analysis with the AeroModule or with the surrogate model
Bayes_full = 0; % 0: use and/or set-up surrogate model (PCE); 1: run full model for Bayes (Computationally expensive!)

% If Bayes_full = 0, we need to specify options for loading a surrogate model
Surrogate_model_type = 0; % 0: Uses a stored PCE surrogate model, 1: create surrogate model

% Options for loading a surrogate model
Surrogate_model_filename = 'StoredSurrogates/NM80_calibrate/PCE_LARS_CL_LHS40.mat'; % Specify the surrogate model file to be used

% Options for creating a surrogate model
% These are used if Bayes_full = 0 and Surrogate_model_type = 1
MetaOpts.Type = 'Metamodel';
MetaOpts.MetaType = 'PCE';
MetaOpts.Method = 'LARS'; % Quadrature, OLS, LARS

MetaOpts.ExpDesign.Sampling = 'LHS';
MetaOpts.ExpDesign.NSamples = 40;
MetaOpts.Degree = 1:4;
MetaOpts.TruncOptions.qNorm = 0.75;

%% MCMC options
% MCMC parameters
Solver.Type = 'MCMC';
% MCMC algorithms available in UQLab
MCMC_type = 'AIES'; % 'MH': Metropolis-Hastings, 'AM': Adaptive Metropolis
                    % 'AIES': Affine invariant ensemble, 'HMC': Hamilton Monte Carlo

switch MCMC_type
    
    case 'MH'
        
        Solver.MCMC.Sampler = 'MH';
        Solver.MCMC.Steps = 1e2;
        Solver.MCMC.NChains = 1e2;
        Solver.MCMC.T0 = 1e1;
        
    case 'AM'
        Solver.MCMC.Sampler = 'AM';
        Solver.MCMC.Steps = 1e2;
        Solver.MCMC.NChains = 1e2;
        Solver.MCMC.T0 = 1e1;
        Solver.MCMC.Epsilon = 1e-2;
        
    case 'AIES'
        Solver.MCMC.Sampler = 'AIES';
        Solver.MCMC.Steps = 5e2;
        Solver.MCMC.NChains = 1e2;
        Solver.MCMC.a = 5;
        
    case 'HMC'
        Solver.MCMC.Sampler = 'HMC';
        Solver.MCMC.LeapfrogSteps = 1e3;
        Solver.MCMC.LeapfrogSize = 0.01;
        Solver.MCMC.Mass = 100;
        
    otherwise
        error('wrong MCMC type provided');
        
end
