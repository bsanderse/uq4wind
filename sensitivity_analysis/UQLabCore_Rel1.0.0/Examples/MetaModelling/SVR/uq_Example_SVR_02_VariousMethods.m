%% SVR METAMODELING: VARIOUS METHODS
%
% This example showcases how to perform Support Vector Machine for 
% Regression (SVR) metamodeling for a simple one-dimensional function,
% using various hyperparameter estimation and optimization methods.

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
% $$Y(x) = x \sin(x), \; x \in [0, 15]$$
%
% In UQLab, the model can be specified directly as a string
% (with vectorized operation):
ModelOpts.mString = 'X.*sin(X)';
ModelOpts.isVectorized = true;

myModel = uq_createModel(ModelOpts);

%% 3 - PROBABILISTIC INPUT MODEL
%
% The probabilistic input model consists of a uniform random variable:
%
% $$X \sim \mathcal{U}(0, 15)$$

%%
% Specify the probabilistic model of the input variable:
InputOpts.Marginals.Type = 'Uniform';
InputOpts.Marginals.Parameters = [0 15];

%%
% Create an INPUT object:
myInput = uq_createInput(InputOpts);

%% 4 - EXPERIMENTAL DESIGN AND MODEL RESPONSES
%
% The experimental design and the model responses are calculated 
% and are later used for creating the SVR metamodel.
%
% Generate $10$ sample points from the input distribution
% using the latin hypercube sampling:
X = uq_getSample(10,'LHS');

%%
% Evaluate the corresponding model responses:
Y = uq_evalModel(X);

%% 5 - SVR METAMODELS
%
% Select SVR as the metamodeling tool:
MetaOpts.Type = 'Metamodel';
MetaOpts.MetaType = 'SVR';

%% 
% Use the experimental design and corresponding model responses 
% generated earlier:
MetaOpts.ExpDesign.X = X;
MetaOpts.ExpDesign.Y = Y;

%%
% Select the 'Matern 5/2' kernel family:
MetaOpts.Kernel.Family = 'Matern-5_2';

%%
% The previous set of options are fixed and will be used for all 
% the SVR metamodels created below.
% The optimization and hyperparameter estimation methods varies 
% for each SVR metamodel.

%%
% * Create an SVR model using the Span leave-one-out (LOO) error estimator
%   and BFGS optimization method:
disp(['> Estimation Method: LOO span estimate, ',...
    'Optimization method: BFGS'])
MetaOpts.EstimMethod = 'SpanLOO';
MetaOpts.Optim.Method = 'BFGS';
mySVR_Span_BFGS = uq_createModel(MetaOpts);

%%
% * Create an SVR model using the Span LOO error estimator
%   and HCMAES optimization method:
disp(['> Estimation Method: LOO span estimate, ',...
    'Optimization method: HCMAES'])
MetaOpts.Optim.Method = 'HCMAES';
mySVR_Span_HCMAES = uq_createModel(MetaOpts);

%%
% * Create an SVR model using the Smooth LOO error estimator
%   and BFGS optimization method:
disp(['> Estimation Method: Smoothed LOO span estimate, ',...
    'Optimization method: BFGS'])
MetaOpts.EstimMethod = 'SmoothLOO';
MetaOpts.Optim.Method = 'BFGS';
mySVR_Smooth_BFGS = uq_createModel(MetaOpts);

%%
% * Create an SVR model using the Smooth LOO error estimator 
%   and HCE optimization method:
disp(['> Estimation Method: Smoothed LOO span estimate, ',...
    'Optimization method: HCE'])
MetaOpts.Optim.Method = 'HCE';
mySVR_Smooth_HCE = uq_createModel(MetaOpts);

%%
% * Create an SVR model using the Cross Validation error estimator 
%   and BFGS optimization method:
disp(['> Estimation Method: Cross-Validation, ',...
    'Optimization method: BFGS'])
MetaOpts.EstimMethod = 'CV';
MetaOpts.Optim.Method = 'BFGS';
MetaOpts.Optim.maxIter = 30;
mySVR_CV_BFGS = uq_createModel(MetaOpts);

%% 
% * Create a SVR model using the Cross-Validation (CV) error estimator
%   and HGA optimization method:
disp(['> Estimation Method: Cross-Validation, ',...
    'Optimization method: HGA'])
MetaOpts.Optim.Method = 'HGA';
mySVR_CV_HGA = uq_createModel(MetaOpts) ;

%% 6 - COMPARISON OF THE METAMODELS
%
% Create a validation set of size $10^3$:
Xval = linspace(0, 15, 1e3)';

%%
% Evaluate the full model responses at the validation set points:
Yval = uq_evalModel(myModel,Xval);

%%
% Evaluate the corresponding responses for each of the generated SVR
% metamodel.
Y_Span_BFGS  = uq_evalModel(mySVR_Span_BFGS,Xval);
Y_Span_HCMAES = uq_evalModel(mySVR_Span_HCMAES,Xval);
Y_Smooth_BFGS = uq_evalModel(mySVR_Smooth_BFGS,Xval);
Y_Smooth_HCE = uq_evalModel(mySVR_Smooth_HCE,Xval);
Y_CV_BFGS = uq_evalModel(mySVR_CV_BFGS,Xval);
Y_CV_HGA = uq_evalModel(mySVR_CV_HGA,Xval);

%%
% Comparative plots of the SVR predictors are generated.
% They are divided into three groups based on the hyperparameter
% estimation methods.
%
% * LOO Span estimate:
myColors = uq_cmap(2); % use UQLab colormap
uq_figure('Position', [50 50 500 400])

uq_plot(Xval, Yval, 'k')
hold on
uq_plot(Xval, Y_Span_BFGS, 'Color', myColors(1,:))
uq_plot(Xval, Y_Span_HCMAES, '--', 'Color', myColors(2,:), 'LineWidth', 3)
uq_plot(X, Y, 'ro')
hold off

axis([0 15 -15 25])
legend({'Original model', 'SVR, optim. method: BFGS',...
    'SVR, optim. method: HCMAES', 'Observations'},...
    'Location', 'north', 'Interpreter', 'latex')
uq_setInterpreters(gca)
xlabel('$\mathrm{X}$')
ylabel('$\mathrm{\widehat{Y}(x)}$')
title('LOO span estimate')
set(gca, 'Fontsize', 14, 'Linewidth', 2)

%%
% * Smooth LOO span estimate:
uq_figure('Position',[50 50 500 400])

uq_plot(Xval, Yval, 'k')
hold on
uq_plot(Xval, Y_Smooth_BFGS, 'Color', myColors(1,:));
uq_plot(Xval, Y_Smooth_HCE, '--', 'Color', myColors(2,:), 'LineWidth', 3)
uq_plot(X, Y, 'ro')
hold off

axis([0 15 -15 25])
legend({'Original model', 'SVR, optim. method: BFGS',...
    'SVR, optim. method: HCE', 'Observations'},...
    'Location', 'north', 'Interpreter', 'latex')
uq_setInterpreters(gca)
xlabel('$\mathrm{X}$')
ylabel('$\mathrm{\widehat{Y}(x)}$')
title('Smoothed LOO span estimate')
set(gca, 'Fontsize', 14, 'Linewidth', 2)

%%
% * Cross-validation (CV)-based estimation:
uq_figure('Position', [50 50 500 400])

uq_plot(Xval, Yval, 'k')
hold on
uq_plot(Xval, Y_CV_BFGS, 'Color', myColors(1,:))
uq_plot(Xval, Y_CV_HGA, '--', 'Color', myColors(2,:), 'LineWidth', 3);
uq_plot(X, Y, 'ro')
hold off

axis([0 15 -15 25])
legend({'Original model', 'SVR, optim. method: BFGS',...
    'SVR, optim. method: HGA', 'Observations'},...
    'Location', 'North', 'Interpreter', 'latex')
uq_setInterpreters(gca)
xlabel('$\mathrm{X}$')
ylabel('$\mathrm{\widehat{Y}(x)}$')
title('(LOO) Cross-Validation')
set(gca, 'Fontsize', 14, 'Linewidth', 2)