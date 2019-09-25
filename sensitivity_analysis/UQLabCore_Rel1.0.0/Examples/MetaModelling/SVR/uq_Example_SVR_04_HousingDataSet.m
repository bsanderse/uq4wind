%% SVR METAMODELING: BOSTON HOUSING DATA SET
%
% This example showcases how to use Support Vector Machine for Regression
% (SVR) metamodeling on an existing data set. 
% The data set under consideration is related to the housing prices
% in Boston, a standard machine learning benchmark data set.
% 
% For more information, see: Harrison, D. and D. L. Rubinfeld, D.L. (1978).
% Hedonic prices and the demand for clean air.
% J. Environ. Economics & Management, vol. 5, pp. 81-102.

%% 1 - INITIALIZE UQLAB
%
% Clear all variables from the workspace, set the random number generator
% for reproducible results, and initialize the UQLab framework:
clearvars
rng(1,'twister')
uqlab

% This example will be using SMO QP optimizer. This feature only appears in Matlab R2015b. 
% Check its availability. Otherwise, issue a warning and revert to interior point (IP) ;
QPSolver = 'SMO' ;
if ~exist('fitrsvm','file')
    QPSolver = 'IP' ;
    warning('SMO is not available on this version of Matlab. Setting QP Solver to IP');
    % IP with a large training set may be slow to converge
    fprintf('This example will take some time!') ;  
end
%% 2 - RETRIEVE DATA SETS
%
% The housing data set is stored in a MAT-file in the following location:
FILELOCATION = fullfile(uq_rootPath, 'Examples', 'SimpleDataSets',...
    'Boston_Housing');

%%
% Read the data set and store the contents in matrices:
load(fullfile(FILELOCATION,'housing.mat'), 'X', 'Y')

%% 3 - TRAINING AND VALIDATION SETS
%
% Retrieve the size of the experimental design:
[N,M] = size(X);

%%
% Use $85\%$ of the data for training and the rest for validation:
Ntrain = floor(0.85*N);
Nval = N - Ntrain;

%% 4 - SVR METAMODELS
%
% Three SVR metamodels are created using different
% training sets of the same size (|Ntrain|) randomly sampled
% from the whole housing data set.
% The steps (repeated three times) are as follows:
%
% # Randomly split the available data into a traning and a validation sets
% # Define and create a SVR metamodel based on the training set
% # Evaluate the SVR metamodel at the validation set points

%%
% Initialize the results matrices:
Yval = zeros(Nval,3);
YSVR = zeros(Nval,3);
legend_text = {};

%%
% Loop over different sample of training data:
for iter = 1:3

    % Randomly split the data into a training and a validation sets 
    idx_train = randperm(N, Ntrain);
    idx_val = setdiff(1:N,idx_train);
    Xtrain = X(idx_train,:);
    Ytrain = Y(idx_train);
    Xval = X(idx_val,:);
    Yval(:,iter) = Y(idx_val);
    
    % Select SVR as the metamodeling technique:
    MetaOpts.Type = 'Metamodel';
    MetaOpts.MetaType = 'SVR';
    
    % Assign the training data set as the experimental design
    MetaOpts.ExpDesign.X = Xtrain;
    MetaOpts.ExpDesign.Y = Ytrain;
    
    % Use the L1-penalization scheme:
    MetaOpts.Loss = 'l1-eps';
    
    % Use an anisotropic kernel:
    MetaOpts.Kernel.Isotropic = 0;
    
    % Use the BFGS method for the hyperparameters calibration:
    MetaOpts.Optim.Method = 'BFGS';
    
    % Use the Sequential Minimal Optimization (SMO) algorithm to solve the
    % quadratic programming problem:
    MetaOpts.QPSolver = QPSolver ;

    % Create the SVR metamodel:
    mySVR = uq_createModel(MetaOpts);
    
    % Evaluate the SVR metamodel at the validation set points
    YSVR(:,iter) = uq_evalModel(mySVR,Xval);
    
    % Create a text for legend
    legend_text = [legend_text, sprintf('Iteration %i', iter)];
end

%% 5 - VALIDATION
%
% Compare the SVR predictions at the validation set points
% against the true values over the $3$ repetitions:
uq_figure('Position', [50 50 500 400])

plot(Yval, YSVR, '+')
hold on
plot([0 60], [0 60], 'Color', 'k')
hold off
box on
grid on

axis equal
axis([0 60 0 60])

colormap(uq_cmap(5))

legend(legend_text, 'Location', 'southeast', 'Interpreter', 'latex')
uq_setInterpreters(gca)
xlabel('$\mathrm{Y_{true}}$')
ylabel('$\mathrm{Y_{SVR}}$')
set(gca, 'FontSize', 14, 'Linewidth', 2)