%% model description 

% name of Matlab file representing the model
Model.mHandle = @aero_module;

% optionally, one can pass parameters to the model
P = getParameterAeroModule(); % Get values of deterministic parameters used in AERO module
Model.Parameters = P;

% Add paths for dependent routines located in the directories 'NURBS',
% 'AEROmoduleWrapper' and 'Geometry'
addpath('C:\Users\pkumar\Dropbox\WindTrue\windtrue\AEROmoduleWrapper\')
addpath('C:\Users\pkumar\Dropbox\WindTrue\windtrue\NURBS\')
addpath('C:\Users\pkumar\Dropbox\WindTrue\windtrue\Geometry\')

%% list of UQ methods to be used for analysis

% specify a list of options from the following list:
% methods = {'MC','PCE_Quad','PCE_OLS','PCE_LARS'};
methods = {'PCE_Quad','PCE_OLS','PCE_LARS'};

% % for Monte Carlo, specify number of times to repeat MC-based methods to obtain 'nice' convergence
% % graphs
% MC_repeat = 1;
% % number of samples with MC
% NsamplesMC = [1e1 1e2 1e3 1e4];

% for PCE-Quad, specify the polynomial degrees to be tested
DegreesQuad = 1:8; %[1 2 3 4 5 6];

% for PCE-OLS:
NsamplesOLS = [5 10 20 30 40]; % if not specified, the number of samples from Quad is taken
OLS_repeat = 1; % like MC_repeat

% for PCE-LARS:
NsamplesLARS = [5 10 20 30 40]; % if not specified, the number of samples from Quad is taken
LARS_repeat = 1; % like MC_repeat


%% input description
% X1: CL
% X2: V

% number of random variables
ndim = 2;

% marginal distribution X1
Input.Marginals(1).Name = 'CL';
Input.Marginals(1).Type = 'Uniform'; 
Input.Marginals(1).Parameters = [0.1 1];
Input.Marginals(1).Bounds = [0.1 1]; 

% marginal distribution X2
Input.Marginals(2).Name = 'V';
Input.Marginals(2).Type = 'Weibull'; 
Input.Marginals(2).Parameters = [10 2]; % scale and shape parameter

%% exact mean/std for error computation, if known:
% mean_exact = ;
% var_exact  = ;
% std_exact  = sqrt(var_exact);

% mean_ref and std_ref should be taken different from mean_exact and
% std_exact if it equals zero
% mean_ref   = mean_exact;
% std_ref    = std_exact;
 
%% exact Sobol indices
% variances:
% D1 = ;
% D2 = ;
% D3 = ;
% D12 = ;
% D13 = ;
% D23 = ;
% D123 = ;
 
% first order indices
% S1 = D1/var_exact;
% S2 = D2/var_exact;
% S3 = D3/var_exact;
% S12 = D12/var_exact;
% S13 = D13/var_exact;
% S23 = D23/var_exact;
% S123 = D123/var_exact;