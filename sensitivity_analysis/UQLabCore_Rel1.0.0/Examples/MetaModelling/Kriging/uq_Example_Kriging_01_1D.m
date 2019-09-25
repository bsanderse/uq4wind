%% KRIGING METAMODELING: ONE-DIMENSIONAL EXAMPLE
%
% This example showcases how to perform Kriging metamodeling
% on a simple one-dimensional function
% using various types of correlation families.

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
% Generate $10$ samples of $X$ using the Sobol' sequence sampling:
X = uq_getSample(10,'Sobol');

%%
% Evaluate the corresponding model responses:
Y = uq_evalModel(X);

%% 5 - KRIGING METAMODELS
%
% Select the metamodeling tool and the Kriging module:
MetaOpts.Type = 'Metamodel';
MetaOpts.MetaType = 'Kriging';

%% 
% Use the experimental design and corresponding model responses 
% generated earlier:
MetaOpts.ExpDesign.X = X;
MetaOpts.ExpDesign.Y = Y;

%% 5.1 Matérn correlation
%
% Perform ordinary Kriging
% (i.e., Kriging with a constant but unknown trend):
MetaOpts.Trend.Type = 'ordinary';

%%
% Set the correlation family to Matérn 5/2:
MetaOpts.Corr.Family = 'matern-5_2';

%%
% Use the maximum likelihood estimation (MLE) 
% to estimate the hyperparameters:
MetaOpts.EstimMethod = 'ML';

%%
% Create the Kriging metamodel:
myKrigingMat = uq_createModel(MetaOpts);

%% 
% Print out a report on the resulting Kriging object:
uq_print(myKrigingMat)

%%
% Plot a representation of the mean and the 95% confidence bounds of the
% Kriging predictor:
uq_display(myKrigingMat)

%% 4.2 Linear correlation
%
% Create another Kriging metamodel with the same configuration
% options but use a linear correlation family instead: 
MetaOpts.Corr.Family = 'linear';
myKrigingLin = uq_createModel(MetaOpts);

%% 4.3 Exponential correlation
%
% Finally, create a Kriging metamodel using the exponential correlation family:
MetaOpts.Corr.Family = 'exponential';
myKrigingExp = uq_createModel(MetaOpts);

%% 6 - METAMODELS VALIDATION
%
% Create a validation set with $10^3$ samples over a regular grid:
Xval = uq_getSample(10^3,'grid');

%%
% Evaluate the true model responses for the validation set:
Yval = uq_evalModel(myModel,Xval);

%%
% Evaluate the Kriging surrogate predictions on the validation set.
% In each case, both the mean and the variance of the Kriging predictor
% are calculated:
[YMeanMat,YVarMat] = uq_evalModel(myKrigingMat,Xval);
[YMeanLin,YVarLin] = uq_evalModel(myKrigingLin,Xval);
[YMeanExp,YVarExp] = uq_evalModel(myKrigingExp,Xval);

%%
% Compare the mean prediction of each Kriging metamodel on the validation
% set (also taking into account the true model responses):
uq_figure('position',[50 50 500 400])

uq_plot(Xval, Yval, 'k', 'LineWidth', 2)
hold on
uq_plot(Xval, YMeanMat, 'LineWidth', 2)
uq_plot(Xval, YMeanLin, 'LineWidth', 2)
uq_plot(Xval, YMeanExp, 'LineWidth', 2)
uq_plot(X, Y, 'ko','MarkerFaceColor','k')
hold off

xlim([0 15])
ylim([-20 30])

uq_legend({'True model', 'Kriging, R: Matern 5/2',...
    'Kriging, R: Linear', 'Kriging, R: Exponential', 'Observations'},...
    'Location', 'northwest', 'FontSize', 14, 'Interpreter', 'latex')

uq_setInterpreters(gca)
xlabel('$\mathrm{X}$')
ylabel('$$\mathrm{\mu_{\widehat{Y}}}$$')
set(gca, 'FontSize', 14)

%%
% Finally, compare the variance that is predicted by each Kriging
% metamodel:
uq_figure('position',[50 50 500 400])

uq_plot(Xval, YVarMat, 'LineWidth', 2)
hold on
uq_plot(Xval, YVarLin, 'LineWidth', 2)
uq_plot(Xval, YVarExp, 'LineWidth', 2)
uq_plot(X, zeros(size(X)), 'ko', 'MarkerFaceColor','k')
hold off

xlim([0 15])
ylim([0 50])

legend({'Kriging, R: Matern 5/2', 'Kriging, R: Linear',...
    'Kriging, R: Exponential', 'Observations'},...
    'Location', 'north', 'FontSize', 14, 'Interpreter', 'latex')

uq_setInterpreters(gca)
xlabel('$\mathrm{X}$')
ylabel('$$\mathrm{\sigma^2_{\widehat{Y}}}$$')
set(gca, 'FontSize', 14)