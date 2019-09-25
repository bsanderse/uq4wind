%% INPUT MODULE: MARGINALS AND COPULA
%
% This example showcases how to define a probabilistic input model
% with or without a copula dependency.

%% 1 - INITIALIZE UQLAB
%
% Clear all variables from the workspace, set the random number generator
% for reproducible results, and initialize the UQLab framework:
clearvars
rng(1,'twister')
uqlab

%% 2 - PROBABILISTIC INPUT MODEL
%
% The probabilistic input model consists of two random variables:
%
% $$X_1 \sim \mathcal{N}(\mu = 0, \sigma = 1)$$
%
% $$X_2 \sim \mathcal{B}(r = 1, s = 3)$$

% Specify the marginals:
InputOpts.Marginals(1).Type = 'Gaussian';
InputOpts.Marginals(1).Parameters = [0 1];
InputOpts.Marginals(2).Type = 'Beta';
InputOpts.Marginals(2).Parameters = [1 3];

%% 2.1 Without dependency
%
% By default, the variables are considered independent.

%%
% Create an INPUT object:
myInput1 = uq_createInput(InputOpts);

%%
% Print a summary of the INPUT object:
uq_print(myInput1)

%% 2.2 With dependency: Gaussian copula
%
% The marginal distributions of the random variables variables
% of the probabilistic input model are already defined inside
% the structure |Input|.
% Dependency using a Gaussian copula is added as follows:
InputOpts.Copula.Type = 'Gaussian';
InputOpts.Copula.RankCorr = [1 0.8; 0.8 1]; % Spearman correlation matrix

%%
% Create a dependent INPUT object under the specified Gaussian copula:
myInput2 = uq_createInput(InputOpts) ;

%%
% Print a report on the INPUT object:
uq_print(myInput2)

%% 3 - COMPARISON OF THE GENERATED INPUT MODELS
%
% Each of the create INPUT objects can be quickly visualized
% using the function |uq_display|:
uq_display(myInput1)
uq_display(myInput2)