%% SVR METAMODELING: ONE-DIMENSIONAL EXAMPLE
%
% This example showcases how to perform Support Vector Machine
% for Regression (SVR) metamodeling on a simple one-dimensional function
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
% Specify the distribution of the input variables:
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
% Use the experimental design generated above:
MetaOpts.ExpDesign.X = X;
MetaOpts.ExpDesign.Y = Y;

%%
% Select the L2-SVR method (mean-square loss function):
MetaOpts.Loss = 'l2-eps';

%%
% Use a hybrid genetic algorithm (GA) for optimizing the hyperparameters:
MetaOpts.Optim.Method = 'HGA';

%% 5.1 Matérn 5/2 kernel
%
% Select the 'Matern 5/2' kernel family:
MetaOpts.Kernel.Family = 'Matern-5_2';

%%
% Create the SVR metamodel with 'Matern 5/2' kernel:
mySVR_matern = uq_createModel(MetaOpts);

%% 
% Print a report on the main features of the resulting L2-SVR model:
uq_print(mySVR_matern)

%%
% Visualize the SVR predictor and the $95\%$ confidence bounds
% (only available for one-dimensional and two-dimensional functions):
uq_display(mySVR_matern)

%% 5.2 Gaussian kernel
%
% Create another SVR metamodel using a Gaussian kernel:
MetaOpts.Kernel.Family = 'Gaussian';
mySVR_exp = uq_createModel(MetaOpts);

%% 5.3 User-defined kernel
%
% A user-defined (custom) kernel is built up as a mixture of a Gaussian
% and a polynomial kernel:
%
% $k(x,x') = \theta_3 \cdot \exp[||x-x'||^2/ ( 2\, \theta_1)^2 ] + (1-\theta_3) \cdot (x \cdot x' + \theta_2 )^5$
%
% It is defined using a handle to a function written in an m-file
% (shipped with UQLab):
MetaOpts.Kernel.Handle = @uq_MixedSVRKernel;

%% 
% Initialize the parameters of the kernel:
MetaOpts.Hyperparameters.theta = [0.5 1 0.5];

%%
% Set bounds on the search space for the hyperparameters calibration:
MetaOpts.Optim.Bounds.C = [10; 1000];
MetaOpts.Optim.Bounds.epsilon = [1e-3; 1];
MetaOpts.Optim.Bounds.theta = [0.1 1e-3 1e-3; 5 10 1];

%%
% Create the SVR metamodel with custom kernel:
mySVR_cusK = uq_createModel(MetaOpts);

%% 7 - VALIDATION OF THE METAMODELS
%
% Create a validation sample of size $10^3$:
Xval = linspace(0, 15, 1000)';

%%
% Evaluate the model at the validation sample points:
Yval = uq_evalModel(myModel,Xval);

%%
% Evaluate the corresponding predictions for each
% of the three SVR metamodels:
Y_mat = uq_evalModel(mySVR_matern,Xval);
Y_exp = uq_evalModel(mySVR_exp,Xval);
Y_cusK = uq_evalModel(mySVR_cusK,Xval);

%%
% Create a comparative plot of the SVR predictor of each metamodel
% as well as the output of the true model:
myColors = uq_cmap(4); % use UQLab colormap
uq_figure('position', [50 50 500 400])

uq_plot(Xval, Yval, 'k', 'LineWidth', 2)
hold on
uq_plot(Xval, Y_mat,'LineWidth', 2)
uq_plot(Xval, Y_exp, 'LineWidth', 2)
uq_plot(Xval, Y_cusK,'lineWidth', 2)
uq_plot(X, Y, 'ko','MarkerFaceColor', 'k')
hold off
axis([0 15 -20 20])

legend({'Full model', 'Mat{\''e}rn 5/2', 'Exponential',...
    'Custom kernel', 'Observations'},...
    'Location', 'southwest', 'Interpreter', 'latex')

uq_setInterpreters(gca)
xlabel('$\mathrm X$')
ylabel('$\mathrm{\widehat{Y}(x)}$')
set(gca, 'FontSize', 14)