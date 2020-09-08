%% Linear model calibration with vector output and
% calibration of discrepancy parameter sigma

clc
close all
clearvars

rng default

%% input settings
p = 2; % dimension of parameters to be calibrated
d = 2; % dimension of model output
N_data = 20; % dimension of data points per output, total data d*N

dom = [0,2]; % domain of x

%% plot settings
fontsize = 14;
fontname = 'Helvetica';
set(0,'defaultlinelinewidth',2)

%% generate data and exact solution

% generate artificial measurement data by using the linear model
% with Gaussian noise on top of it

% stdev in measurement error - should also be used in likelihood in Bayes
sigma1 = 0.5;
sigma2 = 0.2;

% exact value for beta is only used to generate artificial measurement data, and to
% plot exact solution
beta_exact = [0.3 0.7]; % row vector of length p

% measurement data, size N_data*d,
z_data  = linearmodel_vectoroutput(beta_exact) + [sigma1*randn(N_data,1) sigma2*randn(N_data,1)];

% generate exact solution for plotting purposes
y_exact = linearmodel_vectoroutput(beta_exact);


%% Bayesian solution with UQLab

%% initialize UQlab

% add path
addpath(genpath('../../../UQLabCore_Rel1.3.0/'));
% start uqlab
uqlab

% model settings
ModelOpts.mFile = 'linearmodel_vectoroutput';
ModelOpts.isVectorized = true;
myForwardModel = uq_createModel(ModelOpts);

% prior
PriorOpts.Marginals(1).Name = 'beta0';
PriorOpts.Marginals(1).Type = 'Uniform';
PriorOpts.Marginals(1).Parameters = [0,1];
%
PriorOpts.Marginals(2).Name = 'beta1';
PriorOpts.Marginals(2).Type = 'Uniform';
PriorOpts.Marginals(2).Parameters = [0,1];
%
myPriorDist = uq_createInput(PriorOpts);


% display input properties
% uq_print(myPriorDist);
% uq_display(myPriorDist);


%% data
Data(1).y = z_data(:,1); % specify column vectors
Data(1).Name = 'Measurements of y1';
Data(1).MOMap = 1; % Model Output Map

Data(2).y = z_data(:,2);
Data(2).Name = 'Measurements of y2';
Data(2).MOMap = 2; % Model Output Map


%% likelihood and prior for discrepancy
DiscrepancyPriorOpts1.Name = 'Prior of sigma 1';
DiscrepancyPriorOpts1.Marginals(1).Name = 'Sigma21';
DiscrepancyPriorOpts1.Marginals(1).Type = 'Uniform';
DiscrepancyPriorOpts1.Marginals(1).Parameters = [0, 4*sigma1^2];
DiscrepancyPrior1 = uq_createInput(DiscrepancyPriorOpts1);

DiscrepancyOpts(1).Type = 'Gaussian';
DiscrepancyOpts(1).Prior = DiscrepancyPrior1;

DiscrepancyPriorOpts2.Name = 'Prior of sigma 2';
DiscrepancyPriorOpts2.Marginals(1).Name = 'Sigma22';
DiscrepancyPriorOpts2.Marginals(1).Type = 'Uniform';
DiscrepancyPriorOpts2.Marginals(1).Parameters = [0, 4*sigma2^2];
DiscrepancyPrior2 = uq_createInput(DiscrepancyPriorOpts2);

DiscrepancyOpts(2).Type = 'Gaussian';
DiscrepancyOpts(2).Prior = DiscrepancyPrior2;

%% Bayes options
Solver.Type = 'MCMC';
% Adaptive Metropolis:
Solver.MCMC.Sampler = 'AM';
Solver.MCMC.Steps = 1e4;
Solver.MCMC.NChains = 1e2;
Solver.MCMC.T0 = 1e2;
%     Solver.MCMC.Proposal.PriorScale = 0.1;
% AIES:
%     Solver.MCMC.Sampler = 'AIES';

% show MCMC chain convergence:
%     Solver.MCMC.Visualize.Parameters = 1;
%     Solver.MCMC.Visualize.Interval = 100;

BayesOpts.Data  = Data;
BayesOpts.Type  = 'Inversion';
BayesOpts.Discrepancy = DiscrepancyOpts;
BayesOpts.Solver = Solver;
BayesOpts.Prior = myPriorDist;

%% perform MCMC
myBayesianAnalysis = uq_createAnalysis(BayesOpts);

%% postprocessing
% text output of Bayesian analysis
uq_print(myBayesianAnalysis)
% graphical display of posterior
uq_display(myBayesianAnalysis)

% All post-processing results generated by the uq_postProcessInversion function are stored
% in the myBayesianAnalysis.Results.PostProc structure.

% get the Maximum a Posteriori value
uq_postProcessInversion(myBayesianAnalysis,'pointestimate','Mean');
beta_UQLab_mean = myBayesianAnalysis.Results.PostProc.PointEstimate.X
% get the mean posterior value
uq_postProcessInversion(myBayesianAnalysis,'pointEstimate', 'MAP');
beta_UQLab_MAP = myBayesianAnalysis.Results.PostProc.PointEstimate.X

y_UQLab_MAP = linearmodel_vectoroutput(beta_UQLab_MAP(1:2));
y_UQLab_mu  = myBayesianAnalysis.Results.PostProc.PostPred.model.pointEstimateRun; % posterior predictive