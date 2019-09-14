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
methods = {'PCE_LARS'};
% % for Monte Carlo, specify number of times to repeat MC-based methods to obtain 'nice' convergence
% % graphs
% MC_repeat = 1;
% % number of samples with MC
% NsamplesMC = [1e1 1e2 1e3 1e4];

% % for PCE-Quad, specify the polynomial degrees to be tested
% DegreesQuad = 1:4; %[1 2 3 4 5 6];

% % % for PCE-OLS:
% NsamplesOLS = [16]; % if not specified, the number of samples from Quad is taken
% OLS_repeat = 1; % like MC_repeat
% 
% % % for PCE-LARS:
NsamplesLARS = [256 512 1024 2048]; % if not specified, the number of samples from Quad is taken
LARS_repeat = 1; % like MC_repeat

%% number of random variables
ndim = length(P{32}) + length(P{33}) + length(P{34}) + 1; % Number of random control points for twist, chord, thickness

% marginal distributions for each control points
for i = 1:ndim-1
    Input.Marginals(i).Name = ['CP',num2str(i)];
    Input.Marginals(i).Type = 'Uniform'; 
    Input.Marginals(i).Parameters = [-0.5 0.5];
    Input.Marginals(i).Bounds = [-0.5 0.5]; 
end

Input.Marginals(i).Name = 'YAW';
Input.Marginals(ndim).Type = 'Gaussian'; 
Input.Marginals(ndim).Parameters = [P{20}, P{35}];
Input.Marginals(ndim).Bounds = [-10 10]; 