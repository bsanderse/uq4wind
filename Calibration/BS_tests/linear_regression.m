%% Sobol indices computed with Monte Carlo and PCE-type methods

clc
close all
clearvars


%% initialize UQlab

% add path
addpath(genpath('../../UQLabCore_Rel1.3.0/'));
% start uqlab
uqlab

%% 

% generate data
N_data = 50;
x_data = linspace(0,1,N_data)';

beta_exact = [2;3];
y_data  = linearmodel(x_data,beta_exact) + 0.5*randn(N_data,1);

N_exact = 100;
x_exact = linspace(0,1,N_exact)';
y_exact = linearmodel(x_exact,beta_exact);

plot(x_data,y_data,'x');
hold on
plot(x_exact,y_exact,'-');


%% ordinary least squares solution
% design matrix
% X = [ones(N_data,1) x_data];
% automatically create design matrix:
p = 2;
X = zeros(N_data,p);
X(:,1) = linearmodel(x_data,[1 0]);
X(:,2) = linearmodel(x_data,[0 1]);
   
beta_OLS = regress(y_data,X)
% this is the same as 
% beta_OLS = (X' * X)\(X' * y_data)
% or
% beta_OLS = pinv(X)*y_data

y_OLS = linearmodel(x_exact,beta_OLS);
plot(x_exact,y_OLS,'--');

legend('data','exact','OLS');

%% Bayesian solution

% 
ModelOpts.mFile = 'linearmodel';
ModelOpts.isVectorized = true;
% 
myForwardModel = uq_createModel(ModelOpts);
% 
% prior
PriorOpts.Marginals(1).Name = 'beta_0';
PriorOpts.Marginals(1).Type = 'Uniform';
PriorOpts.Marginals(1).Parameters = [1,3];
% 
PriorOpts.Marginals(2).Name = 'beta_1';
PriorOpts.Marginals(2).Type = 'Uniform';
PriorOpts.Marginals(2).Parameters = [2,4];

% 
% PriorOpts.Marginals(3).Name = 'L';               % beam length
% PriorOpts.Marginals(3).Type = 'Constant';
% PriorOpts.Marginals(3).Parameters = 5;           % (m)
% 
% PriorOpts.Marginals(4).Name = 'E';               % Young's modulus
% PriorOpts.Marginals(4).Type = 'LogNormal';
% PriorOpts.Marginals(4).Moments = [30000 4500];   % (MPa)
% 
% PriorOpts.Marginals(5).Name = 'p';               % uniform load
% PriorOpts.Marginals(5).Type = 'Gaussian';
% PriorOpts.Marginals(5).Moments = [0.012 0.012*0.05]; % (kN/m)
% 
myPriorDist = uq_createInput(PriorOpts);
% 
% display input properties
uq_print(myPriorDist);
uq_display(myPriorDist);
% 
% %% data
% myData.y = [12.84; 13.12; 12.13; 12.19; 12.67]/1000; % (m)
% myData.Name = 'Mid-span deflection';
% 
% Bayes options and run
BayesOpts.Type = 'Inversion';
BayesOpts.Data = myData;
myBayesianAnalysis = uq_createAnalysis(BayesOpts);
% 
% %% postprocessing
% uq_display(myBayesianAnalysis)
