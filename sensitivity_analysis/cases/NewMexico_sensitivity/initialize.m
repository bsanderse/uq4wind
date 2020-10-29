%% Turbine input parameters
% Name of Matlab file representing the turbine data for sensitivity
turbineName = 'NewMexico'; %
% check NewMexico.m for definition of uncertainties

% folder used to get reference input files
ref_folder      = 'AEROmodule/NewMexico/reference/';
% folder used to write new input files
current_folder  = 'AEROmodule/NewMexico/current/';

%% Experimental data

% the position of the sections of the experimental data which are used for
% interpolation of the aeromodule results: see NewMexico_readoutput.m
r_exp_data = [0.25 0.35 0.6 0.82 0.92]*2.25; 

%% Forward model description
% Name of Matlab file representing the model
Model.mHandle = @aero_module;
% Quantity of interest
QoI = 'Axial_Force';

% Pass parameters to model via the cell array FixedInputs
[FixedParameters,UncertainInputs] = getParameterAeroModule(turbineName);

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
NsamplesMC = 4; %[8 16 32];

% for PCE_Quad, specify the polynomial degrees to be tested
DegreesQuad = 1:3; %[1 2 3 4 5 6];

% for PCE-OLS:
NsamplesOLS = [8 16 32 64 128]; % if not specified, the number of samples from Quad is taken
OLS_repeat = 1; % like MC_repeat
 
% for PCE-LARS:
NsamplesLARS = [16]; % if not specified, the number of samples from Quad is taken

LARS_repeat = 1; % like MC_repeat


%% check location of ECNAeroModule
path_found  = findAeroModulePath();


