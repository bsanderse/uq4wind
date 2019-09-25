%% PCE METAMODELING: BOOTSTRAP PCE FOR LOCAL ESTIMATION
%
% This example showcases the creation and use of a bootstrap polynomial
% chaos expansion (bPCE).
% The full computational model model is an analytical non-linear, 
% non-monotonic function.

%% 1 - INITIALIZE UQLAB
%
% Clear all variables from the workspace, set the random number generator
% for reproducible results, and initialize the UQLab framework:
clearvars
rng(100,'twister')
uqlab

%% 2 - COMPUTATIONAL MODEL
%
% The full computational model is a simple analytical function defined by:
%
% $$\mathcal{M}(X) = X \sin(X)$$
%
% where $X \in [0, 15]$.

%%
% Create a MODEL object using a string:
ModelOpts.mString = 'X.*sin(X)';
myFullModel = uq_createModel(ModelOpts);

%% 3 - PROBABILISTIC INPUT MODEL
%
% The probabilistic input model consists of a uniform random variable:
%
% $$X \sim \mathcal{U}(0, 15)$$

%%
% Specify the probabilistic model of the input variable:
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [0 15];

%%
% Create the INPUT object:
myInput = uq_createInput(InputOpts);

%% 4 - BOOTSTRAP POLYNOMIAL CHAOS EXPANSION METAMODEL
%
% Select the metamodelling tool and the PCE module:
MetaOpts.Type = 'Metamodel';
MetaOpts.MetaType = 'PCE';

%%
% Specify the degree of the expansion:
MetaOpts.Degree = 11;

%%
% Specify the number of points in the experimental design:
MetaOpts.ExpDesign.NSamples = 15;

%%
% Enable bootstrapping by specifying the number of bootstrap replications:
MetaOpts.Bootstrap.Replications = 100;

%%
% Create the PCE metamodel:
myPCE = uq_createModel(MetaOpts);

%% 5 - VALIDATION OF THE METAMODEL
%
% Create a boostrap validation sample:
Xval = linspace(0, 15, 1000)';

%%
% Evaluate the true model on the validation sample:
Yval = uq_evalModel(myFullModel,Xval);

%%
% Evaluate the PCE metamodel and the corresponding bootstrap replications
% on the validation sample:
[YPCval,YPC_var,YPCval_Bootstrap] = uq_evalModel(myPCE,Xval);

%%
% Create plots with confidence bounds (based on empirical quantiles)
% and bootstrap replications:
uq_figure('Position', [50 50 500 400])
myColors = uq_cmap(3);
pb = plot(Xval, YPCval_Bootstrap', 'linewidth', 0.5,...
    'Color', myColors(3,:));
hold on
p(1) = plot(Xval, Yval, 'Color', myColors(1,:));
p(2) = plot(Xval, YPCval, 'Color', myColors(2,:));
p(3) = plot(Xval, quantile(YPCval_Bootstrap,0.025,2), '--',...
    'Color', myColors(2,:));
p(4) = plot(Xval, quantile(YPCval_Bootstrap,0.975,2), '--',...
    'Color', myColors(2,:));
p(5) = plot(myPCE.ExpDesign.X, myPCE.ExpDesign.Y, 'o',...
    'MarkerSize', 5, 'Color', myColors(1,:),...
    'MarkerFaceColor', myColors(1,:));
hold off
grid on
% Format the axes and add labels and a legend
uq_setInterpreters(gca)
legend([p([5 1:3]) pb(1)],...
    {'Exp Design', 'True', 'PCE', '95\% Confidence Bounds',...
    'Replications'}, 'Interpreter', 'latex', 'Location', 'best',...
    'Fontsize', 14)
xlabel('$\mathrm{X}$')
ylabel('$\mathcal{M}(X)$')