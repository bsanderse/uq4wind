%% This file is used to initialize all settings in order to perform sensitivity study
% of the NM80 (DanAero) case

% Name of Matlab file representing the uncertain parameters
turbineName = 'NM80'; % 'NM80', 'AVATAR'

% folder used to get reference input files
ref_folder      = 'AEROmodule/NM80/reference_adapted_thickness/';
% folder used to write new input files
current_folder  = 'AEROmodule/NM80/current/';

%% Forward model description
% Name of Matlab file representing the model
Model.mHandle = @aero_module;

% Quantity of interest
% QoI = 'Axial_Force';
QoI = 'Sectional_normal_force';

% for sectional normal force, we look at sensitivity with respect to the force
% at the following points, corresponding to the experimental data
% Measurement radial stations, distance from blade root (not rotor center)
r_exp_data = [11.876, 17.820, 28.976, 35.535];
% Alternatively, distance from rotor center:
% r_exp_data = [13.0, 19.0, 30.0, 37.0];  % see e.g. the DanAero MW final report
% r_exp_data = [13.116, 19.06, 30.216, 36.775]; % according to Koen

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


%% Input uncertainties
ndim  = length(UncertainInputs.Marginals);
Input = UncertainInputs;

%% list of UQ methods to be used for analysis

% specify a list of options from the following list:

% methods = {'PCE_OLS'}; % {'MC','PCE_Quad','PCE_OLS','PCE_LARS'};

methods = {'PCE_LARS'}; % {'MC','PCE_Quad','PCE_OLS','PCE_LARS'};


% for MC, specify number of times to repeat MC-based methods to obtain 'nice' convergence
% graphs
MC_repeat = 1;
% number of samples with MC
NsamplesMC = [8 16 32];

% for PCE_Quad, specify the polynomial degrees to be tested
DegreesQuad = 1:3; %[1 2 3 4 5 6];

% % for PCE-OLS:
NsamplesOLS = [8 16 32 64 128]; % if not specified, the number of samples from Quad is taken
OLS_repeat = 1; % like MC_repeat
 
% for PCE-LARS:

NsamplesLARS = [4;8;16;32]; %32]; % if not specified, the number of samples from Quad is taken

LARS_repeat = 5; % like MC_repeat

%% check location of ECNAeroModule
%path_found  = findAeroModulePath();
