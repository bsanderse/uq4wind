%% KRIGING METAMODELING: TREND TYPES
%
% This example showcases how to perform Kriging for a 
% simple one-dimensional function using various trend types. 

%% 1 - INITIALIZE UQLAB
%
% Clear all variables from the workspace,
% set the random number generator for reproducible results,
% and initialize the UQLab framework:
clearvars
rng(100,'twister')
uqlab

%% 2 - COMPUTATIONAL MODEL
%
% The computational model is a simple analytical function defined by:
%
% $$Y = x \sin(x), \; x \in [0, 15]$$
% 
% Specify this model in UQLab using a string with vectorized operation:
ModelOpts.Name = 'XsinX';
ModelOpts.mString = 'X.*sin(X)';
ModelOpts.isVectorized = true;

myModel = uq_createModel(ModelOpts);

%% 3 - PROBABILISTIC INPUT MODEL
%
% The probabilistic input model consists of one uniform random variable:
% 
% $X \sim \mathcal{U}(0, 15)$
%
% Specify the marginal and create a UQLab INPUT object:
InputOpts.Marginals.Type = 'Uniform';
InputOpts.Marginals.Parameters = [0 15];

myInput = uq_createInput(InputOpts);

%% 4 - EXPERIMENTAL DESIGN AND MODEL RESPONSES
%
% Generate $10$ samples of $X$ using the Sobol' sequence sampling:
X = uq_getSample(10,'Sobol');

%%
% Evaluate the corresponding model responses:
Y = uq_evalModel(X);

%% 5 - KRIGING METAMODELS
%
% Three trend types are considered in this example:
%
% * Constant (so-called _ordinary Kriging_)
% * 3rd-degree polynomial
% * Custom type 

%% 5.1 Ordinary Kriging
%
% Select the metamodeling tool and the Kriging module:
MetaOpts.Type = 'Metamodel';
MetaOpts.MetaType = 'Kriging';

%%
% Assign the experimental design and the corresponding model responses:
MetaOpts.ExpDesign.X = X;
MetaOpts.ExpDesign.Y = Y;

%%
% Set the Kriging trend type to ordinary Kriging:
Metaopts.Trend.Type = 'ordinary';

%%
% Set the correlaton family to Matérn 3/2:
MetaOpts.Corr.Family = 'matern-3_2';

%%
% Use cross-validation to estimate the hyperparameters:
MetaOpts.EstimMethod = 'CV';

%%
% Use a hybrid genetic algorithm (GA)-based optimization to solve the
% cross-validation problem:
MetaOpts.Optim.Method = 'HGA';

%%
% In the hybrid genetic algorithm, the final solution of the genetic
% algorithm is used as the starting point of a gradient-based optimization
% method (BFGS).

%% 
% Specify the maximum number of generations, convergence tolerance,
% and population size of the GA-based optimization: 
MetaOpts.Optim.MaxIter = 100;
MetaOpts.Optim.Tol = 1e-6;
MetaOpts.Optim.HGA.nPop = 40;

%%
% Create a Kriging metamodel:
myKrigingOrdinary = uq_createModel(MetaOpts);

%% 5.2 Polynomial Trend
%
% A second Kriging metamodel is created using a 3rd-degree
% polynomial as the trend function,
% while the other options remain the same.
%
% Set the Kriging trend type to polynomial and specify the degree:
MetaOpts.Trend.Type = 'polynomial';
MetaOpts.Trend.Degree = 3;

%%
% Create the second Kriging metamodel:
myKrigingPolynomial = uq_createModel(MetaOpts);

%% 5.3 Custom Trend
%
% A third Kriging metamodel is created with a custom functional
% basis in the trend:
%
% $$f(x) = x^2 + \sqrt{|x|}$$
%
% Set the Kriging trend type to custom and specify the custom function
% using a function handle:
MetaOpts.Trend.Type = 'custom';
MetaOpts.Trend.CustomF = @(X) X.^2 + sqrt(abs(X));

%%
% It is advised to remove fields from the options that are
% not relevant to the metamodel to be created.
% For example, the |.Trend.Degree| option is not relevant here:
MetaOpts.Trend.Degree = [];

%%
% Create a Kriging MODEL object with the specified options:
myKrigingCustom = uq_createModel(MetaOpts);

%% 6 - VALIDATION
%
% Create a validation set with $10^3$ samples over a regular grid:
Xval = uq_getSample(10^3,'grid');

%%
% Evaluate the true model responses on the validation set:
Yval = uq_evalModel(myModel,Xval);

%%
% Evaluate the corresponding responses of each of the three 
% Kriging metamodels:
[YMeanOrdinary,YVarOrdinary] = uq_evalModel(myKrigingOrdinary,Xval);
[YMeanPolynomial,YVarPolynomial] = uq_evalModel(myKrigingPolynomial,Xval);
[YMeanCustom,YVarCustom] = uq_evalModel(myKrigingCustom,Xval);

%%
% For each metamodel, the mean and the variance of the
% Kriging predictor are calculated.

%%
% Compare the mean prediction of each Kriging metamodel on the validation
% set (also taking into account the true model responses):
myColors = uq_cmap(4);
uq_figure('Position', [50 50 500 400])

uq_plot(Xval, Yval, 'k', 'Color', myColors(1,:))
hold on
uq_plot(Xval, YMeanOrdinary, 'Color', myColors(2,:))
uq_plot(Xval, YMeanPolynomial, 'Color', myColors(3,:))
uq_plot(Xval, YMeanCustom, 'Color', myColors(4,:))
uq_plot(X, Y, 'ko', 'MarkerFaceColor', 'k')
hold off

xlim([0 15])
ylim([-15 30])

legend({'True model', 'Ordinary Kriging', '3rd deg. polynomial trend',...
    'Custom trend', 'Observations'},...
    'Location', 'northwest', 'FontSize', 14, 'Interpreter', 'latex')

uq_setInterpreters(gca)
xlabel('$\mathrm{X}$')
ylabel('$\mathrm{\mu_{\widehat{Y}}}$')
set(gca, 'FontSize', 14)

%%
% Finally, compare the variance that is predicted by each Kriging
% metamodel:
uq_figure('Position', [50 50 500 400])

uq_plot(Xval, YVarOrdinary, 'Color', myColors(2,:))
hold on
uq_plot(Xval, YVarPolynomial, 'Color', myColors(3,:))
uq_plot(Xval, YVarCustom, 'Color', myColors(4,:))
uq_plot(X, zeros(size(X)), 'ko', 'MarkerFaceColor', 'k')
hold off

xlim([0 15])
ylim([0 40])

legend({'Ordinary Kriging', '3rd-deg. polynomial trend', 'Custom trend',...
    'Observations'},...
    'Location', 'northwest', 'FontSize', 14, 'Interpreter', 'latex')

uq_setInterpreters(gca)
xlabel('$\mathrm{X}$')
ylabel('$\mathrm{\sigma^2_{\widehat{Y}}}$')
set(gca, 'FontSize', 14)