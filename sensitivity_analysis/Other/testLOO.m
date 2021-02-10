clc
close all
clearvars
rng default


%% initialize UQlab
% % add path
% run('config.m');
% addpath(genpath(UQLab_path));

% start uqlab
uqlab


%% input
Input.Marginals(1).Type = 'Uniform';
Input.Marginals(1).Parameters = [0, 1];
Input.Marginals(1).Bounds = [0, 1]; % not necessary for uniform but useful for plotting


%% model
% Model.mHandle = @aero_module;
Model.mString = 'X(:,1).^2';
Model.isVectorized = true;

%% set-up UQLab model and input
% create and add the model to UQLab
myModel = uq_createModel(Model);
% create input object with UQLab
myInput = uq_createInput(Input) ;


%% metamodel (PCE)
metamodelLARS.FullModel = myModel;
metamodelLARS.Input     = myInput;
metamodelLARS.Type      = 'Metamodel';
metamodelLARS.MetaType  = 'PCE';
metamodelLARS.Method    = 'LARS';
metamodelLARS.Degree    = 1:4; % this automatically switches on degree adaptive PCE
metamodelLARS.TruncOptions.qNorm = 0.75;

metamodelLARS.ExpDesign.Sampling = 'LHS'; % or 'LHS' or 'Sobol' or 'Halton'
metamodelLARS.ExpDesign.NSamples = 4;

myPCE_LARS  = uq_createModel(metamodelLARS);

myPCE_LARS.Error.LOO
myPCE_LARS.Error.ModifiedLOO

%% plot results

% points where the code has been evaluated:
% note that Xsamples includes the constants
Xsamples  = myPCE_LARS.ExpDesign.X;
Ysamples  = myPCE_LARS.ExpDesign.Y;

% evaluate surrogate model at many points:
Ntest  = 100;
domain = getMarginalBounds(myInput.Marginals(1));
X_PCE = linspace(domain(1),domain(2),Ntest)';
Y_PCE = uq_evalModel(myPCE_LARS,X_PCE);

figure
plot(X_PCE,Y_PCE,'-');
hold on
plot(Xsamples,Ysamples,'s');
