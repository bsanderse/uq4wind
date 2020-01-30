%% Sobol indices computed with Monte Carlo and PCE-type methods

clc
close all
clearvars


%% initialize UQlab

% add path
addpath(genpath('../../../UQLabCore_Rel1.3.0/'));
% start uqlab
uqlab


%% generate data and exact solution

% generate artificial measurement data by using the linear model
% with Gaussian noise on top of it

% number of data points
N_data = 20;
% interval [0,1]
x_data = linspace(0,1,N_data)';

% exact value for beta used to generate artificial measurement data, and to
% plot exact solution
beta_exact = [2;3];
% measurement data
y_data  = linearmodel(beta_exact,x_data) + 0.5*randn(N_data,1);

% generate exact solution
N_exact = 100;
x_exact = linspace(0,1,N_exact)';
y_exact = linearmodel(beta_exact,x_exact);

figure(1)
plot(x_data,y_data,'x');
hold on
plot(x_exact,y_exact,'-');


%% ordinary least squares solution
% we now create the design matrix, built with the x_data points
% we can create it automatically by calling linearmodel several times,
% each time with only one of the betas active
p = length(beta_exact); % number of terms to calibrate
A = zeros(N_data,p); % size of design matrix, normally N_data>p
beta_test = eye(p,p);
for i=1:p
    A(:,i) = linearmodel(beta_test(i,:),x_data);
end
% solve the (overdetermined) least-squares problem with the regress command
% this gives the estimate for the beta parameters
beta_OLS = regress(y_data,A)
% note that this is the same as 
% beta_OLS = (X' * X)\(X' * y_data)
% and also the same as
% beta_OLS = pinv(X)*y_data

% with this estimate of beta we get the following solution 
y_OLS = linearmodel(beta_OLS,x_exact);
% alternatively, we can use linearmodel_vectorized, which will return the
% OLS solution at the data points 
% y_OLS = linearmodel_vectorized(beta_OLS,A);
figure(1)
plot(x_exact,y_OLS,'--');


%% Bayesian solution

% model settings
ModelOpts.mFile = 'linearmodel_vectorized';
% pass design matrix as parameter to the M-file
ModelOpts.Parameters   = A; 
ModelOpts.isVectorized = true;
myForwardModel = uq_createModel(ModelOpts);
 
% prior
PriorOpts.Marginals(1).Name = 'beta_0';
PriorOpts.Marginals(1).Type = 'Uniform';
PriorOpts.Marginals(1).Parameters = [1,3];
% 
PriorOpts.Marginals(2).Name = 'beta_1';
PriorOpts.Marginals(2).Type = 'Uniform';
PriorOpts.Marginals(2).Parameters = [2,4];
% 
myPriorDist = uq_createInput(PriorOpts);


% display input properties
uq_print(myPriorDist);
uq_display(myPriorDist);
 

%% likelihood
DiscrepancyOptsKnown.Type = 'Gaussian';
DiscrepancyOptsKnown.Parameters = 1;

%% Bayes options 
Solver.Type = 'MCMC';
Solver.MCMC.Sampler = 'AM';
Solver.MCMC.Steps = 1e3;
Solver.MCMC.NChains = 1e2;
Solver.MCMC.T0 = 1e2;
% Solver.MCMC.Proposal.PriorScale = 0.1;

myData.y       = y_data'; % note: UQLab uses a row vector here
myData.Name    = 'measurement data';
BayesOpts.Data = myData;
BayesOpts.Type = 'Inversion';
BayesOpts.Discrepancy = DiscrepancyOptsKnown;
BayesOpts.Solver = Solver;

%% perform MCMC
myBayesianAnalysis = uq_createAnalysis(BayesOpts);
 
%% postprocessing
% text output of Bayesian analysis
uq_print(myBayesianAnalysis)
% graphical display of posterior
uq_display(myBayesianAnalysis)

% uq_postProcessInversion(myBayesianAnalysis,'pointEstimate', 'MAP')

% get the Maximum a Posteriori value
beta_MAP = myBayesianAnalysis.Results.PostProc.PointEstimate.X

y_MAP = linearmodel(beta_MAP,x_exact);

figure(1)
plot(x_exact,y_MAP,'-.');
legend('measurement data','exact solution','least-squares fit','Bayesian MAP');
