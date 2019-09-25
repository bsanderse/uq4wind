%% KRIGING METAMODELING: TRUSS DATA SET
%
% This example showcases how to perform Kriging metamodeling
% using existing data sets.
% The data sets come from a finite element model of a truss structure
% and are retrieved from different MAT-files.
% The files consist of an experimental design of size $200$
% and a validation basis of size $10^4$.

%% 1 - INITIALIZE UQLAB
%
% Clear all variables from the workspace 
% and initialize the UQLab framework:
clearvars
uqlab

%% 2 - RETRIEVE DATA SETS
%
% The experimental design and the validation basis are stored
% in two separate files in the following location:
FILELOCATION = fullfile(uq_rootPath, 'Examples', 'SimpleDataSets',...
    'Truss_Matlab_FEM');

%%
% Read the experimental design data set file and store the contents 
% in matrices:
load(fullfile(FILELOCATION,'Truss_Experimental_Design.mat'), 'X', 'Y');

%%
% Read the validation basis data set file and store the contents
% in matrices:
load(fullfile(FILELOCATION,'Truss_Validation_Basis.mat'), 'Xval', 'Yval');

%% 3 - KRIGING METAMODEL
%
% Select Kriging as the metamodeling tool:
MetaOpts.Type = 'Metamodel';
MetaOpts.MetaType = 'Kriging';

%% 
% Use experimental design loaded from the data files:
MetaOpts.ExpDesign.X = X;
MetaOpts.ExpDesign.Y = Y;

%%
% Use maximum-likelihood to estimate the hyperparameters:
MetaOpts.EstimMethod = 'ML';

%%
% Use covariance-matrix adaptation evolution strategy optimization
% algorithm (offered by UQLab) for calibratiing the hyperparameters:
MetaOpts.Optim.Method = 'CMAES';

%%
% Provide the validation data set to get the validation error:
MetaOpts.ValidationSet.X = Xval;
MetaOpts.ValidationSet.Y = Yval;

%%
% Create the Kriging metamodel:
myKriging = uq_createModel(MetaOpts);

%% 
% Print a summary of the resulting Kriging metamodel:
uq_print(myKriging)

%% 4 - VALIDATION

%%
% Evaluate the Kriging metamodel at the validation set:
YKRG = uq_evalModel(myKriging,Xval);

%%
% Plot histograms of the true output and the Kriging prediction:
uq_figure('position', [50 50 500 400])
myColors = uq_cmap(2);

uq_histogram(Yval, 'FaceColor', myColors(1,:))
hold on
uq_histogram(YKRG, 'FaceColor', myColors(2,:))
hold off

legend({'True model response', 'Kriging prediction'},...
    'Location', 'northwest')

uq_setInterpreters(gca)
xlabel('$\mathrm{Y}$')
ylabel('Counts')
set(gca, 'Linewidth', 2, 'FontSize', 14)

%% 
% Plot the true vs. predicted values:
uq_figure('position', [50 50 500 400])

hold on 
uq_plot(Yval, YKRG, '+')
uq_plot([min(Yval) max(Yval)], [min(Yval) max(Yval)])
hold off

axis equal 
axis([min(Yval) max(Yval) min(Yval) max(Yval)]) 

uq_setInterpreters(gca)
xlabel('$\mathrm{Y_{true}}$')
ylabel('$\mathrm{Y_{KRG}}$')
set(gca, 'FontSize', 14)

%%
% Print the validation and leave-one-out (LOO) cross-validation errors:
fprintf('Kriging metamodel validation error: %5.4e\n', myKriging.Error.Val)
fprintf('Kriging metamodel LOO error:        %5.4e\n', myKriging.Error.LOO)