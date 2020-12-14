%% This file is used to initialize all settings in order to perform sensitivity study
% of the NewMexico case

% Name of Matlab file representing the turbine data for sensitivity
turbineName = 'NewMexico'; %
% check NewMexico.m for definition of uncertainties

% folder used to get reference input files
ref_folder      = 'AEROmodule/NewMexico/reference/';
% folder used to write new input files
current_folder  = 'AEROmodule/NewMexico/current/';

%% Description of experimental dataset and choice of operating conditions

% NewMexicoData as obtained from Koen Boorsma (TNO)
% folder that contains the experimental data
folder_exp    = fullfile(root_folder,'..','Experimental','NewMexicoData');
% filename with description of datasets
filename_runs = fullfile(folder_exp,'DPN_overview.csv');
% choose the conditions (columns) from the table that are used in the input files
% these should be strings that match the column headers
changing_conditions = {'AIRDENSITY','PITCHANGLE','YAWANGLE','WINDSPEED'}; 
% choose the run for which the sensitivity analysis is to be performed
% note: only a single run is possible at this moment
select_runs = 1002; 

% the position of the sections of the experimental data which are used for
% interpolation of the aeromodule results: see NewMexico_calibrate_readoutput.m
r_exp_data = [0.25 0.35 0.6 0.82 0.92]*2.25;


%% Forward model description
% Name of Matlab file representing the model
Model.mHandle = @aero_module;
% Quantity of interest
% QoI = 'Axial_Force_Blade';
QoI = 'Sectional_normal_force';
% 'mean' or 'full'; in case of 'full', specify:
% * the radial indices (sections to consider in the analysis)
% * the number of revolutions 
% * number of Fourier coefficients that are required
QoI_type = 'mean';
               
% number of revolutions to consider (counting from end of time series)
n_rev = 4;
% number of Fourier coefficients to keep (including mean)
% note: we get (n_fourier-1)*2 + 1 coefficients since there is both a real and
% imaginary component (stored as amplitude and phase angle) for each
% frequency
n_fourier = 1;
% radial indices (blade sections) to consider:
r_index = 1:5;

% Pass parameters to model via the cell array FixedInputs
[FixedParameters,UncertainInputs_NoOC] = getParameterAeroModule(turbineName);

FixedParameters.root_folder    = root_folder;
FixedParameters.ref_folder     = ref_folder;
FixedParameters.current_folder = current_folder;
FixedParameters.QoI            = QoI;
FixedParameters.QoI_type       = QoI_type;

FixedParameters.r_exp          = r_exp_data;

switch QoI_type
    case 'full'
        % store parameters in struct
        FixedParameters.n_rev          = n_rev;
        FixedParameters.n_fourier      = n_fourier;
        FixedParameters.r_index        = r_index;
end

P.FixedParameters = FixedParameters;


%% find the selected run in the CSV file and add operating condition to the model
% - set up a forward and surrogate model for each operating condition

% get the experimental conditions into a table
data_runs     = readtable(filename_runs,'PreserveVariableNames',true);
% number of columns chosen
n_oper = length(changing_conditions);
% total number of columns in table
n_cols = size(data_runs,2); 

% find the index in the table of the current run
i_run = find(data_runs.PNT == select_runs(1));

if (isempty(i_run))
    error(['run with ID ' num2str(select_runs(1)) ' not found']);
end

% loop over the columns of the table and build an object
% OperatingConditions for this run that contains the selected operating
% conditions
k = 1;
conditions_covered = zeros(n_oper,1);
for col=1:n_cols

    % select the column that matches the selected operating conditions
    id_var = strcmp(data_runs.Properties.VariableNames{col},changing_conditions);

    % store in OperatingConditions
    if (any(id_var))
        ind = find(id_var==1);
        OperatingCondition(k).Name = data_runs.Properties.VariableNames{col};
        OperatingCondition(k).Parameters = data_runs(i_run,col).Variables;
        conditions_covered(k) = 1;
        k=k+1;

    end

end
ind_notcovered = find(conditions_covered==0);
if (~isempty(ind_notcovered))
    for k=1:length(ind_notcovered)
        warning(['Variable ' changing_conditions{k} ' not found in the table']);
    end
end

% add Operatingconditions to the structure that contains the
% uncertainties
UncertainInputs = addOperatingConditions(UncertainInputs_NoOC,OperatingCondition);
clear OperatingConditions; 

% total dimension of the uncertain inputs vector, this includes any
% constants
ndim  = length(UncertainInputs.Marginals);



%% finish set-up of forward model

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

methods = {'MC','PCE_Quad'}; % {'MC','PCE_Quad','PCE_OLS','PCE_LARS'};


% for MC, specify number of times to repeat MC-based methods to obtain 'nice' convergence
% graphs
MC_repeat = 1;
% number of samples with MC
NsamplesMC = 4; %[8 16 32];

% for PCE_Quad, specify the polynomial degrees to be tested
DegreesQuad = 1:6; %[1 2 3 4 5 6];

% for PCE-OLS:
NsamplesOLS = 8;%[8 16 32 64 128]; % if not specified, the number of samples from Quad is taken
OLS_repeat = 1; % like MC_repeat
 
% for PCE-LARS:
NsamplesLARS = [80]; % if not specified, the number of samples from Quad is taken

LARS_repeat = 1; % like MC_repeat


%% check location of ECNAeroModule
%path_found  = findAeroModulePath();


