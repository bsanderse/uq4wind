%% model description 

% name of Matlab file representing the model
Model.mHandle = @aero_module;

% optionally, one can pass parameters to the model
P = getParameterAeroModule(); % Get values of deterministic parameters used in AERO module
Model.Parameters = P;
Model.isVectorized = false;

%% Add paths for dependent routines located in the directories 'NURBS','AEROmoduleWrapper' and 'Geometry'
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

% % for PCE-Quad, specify the polynomial degrees to be tested
DegreesQuad = 1:2; %[1 2 3 4 5 6];

% % for PCE-OLS:
NsamplesOLS = [8 16]; % if not specified, the number of samples from Quad is taken
OLS_repeat = 1; % like MC_repeat
% 
% % % for PCE-LARS:
NsamplesLARS = [8 16]; % if not specified, the number of samples from Quad is taken
LARS_repeat = 1; % like MC_repeat

%% input description
% X1: CL
% X2: V

% number of random variables
ndim = length(P{31}) + length(P{30}); % Number of random control points for chord and twist

% marginal distributions for each control points
for i = 1:ndim
    Input.Marginals(i).Name = ['CP',num2str(i)];
    Input.Marginals(i).Type = 'Uniform'; 
    Input.Marginals(i).Parameters = [-0.5 0.5];
    Input.Marginals(i).Bounds = [-0.5 0.5]; 
end