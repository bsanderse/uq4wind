%% INVERSION: BAYESIAN CALIBRATION USING PCE SURROGATE MODEL
% In this example, it is shown how a surrogate model can be constructed
% and then used instead of the original forward model in a Bayesian inversion
% analysis. For a computationally expensive forward model, this approach
% can yield considerable time savings in the analysis.
% The problem considered here is for simple beam model calibration.
%% 1 - INITIALIZE UQLAB
% Clear all variables from the workspace, set the random number generator
% for reproducible results, and initialize the UQLab framework:
clc
close all
clearvars
rng(100,'twister')
uqlab

%% 2.1 - FORWARD MODEL
% The simply supported beam problem is shown in the following figure:
uq_figure
[I,~] = imread('SimplySupportedBeam.png');
image(I)
axis equal
set(gca, 'visible', 'off')

%% 2.2 - FORWARD MODEL (VECTORIZED)
% Define the forward model as a MODEL object using the function
% |uq_SimplySupportedBeam(X)|:
ModelOpts.mFile = 'uq_SimplySupportedBeam';
ModelOpts.isVectorized = true;
myForwardModel = uq_createModel(ModelOpts);

%% 3 - PRIOR DISTRIBUTION OF THE MODEL PARAMETERS
% The prior information about the model parameters is gathered in a
% probabilistic model that includes both known (constant) and unknown
% parameters.

PriorOpts.Marginals(1).Name = 'b';               % beam width
PriorOpts.Marginals(1).Type = 'Constant';
PriorOpts.Marginals(1).Parameters = [0.15];      % (m)

PriorOpts.Marginals(2).Name = 'h';               % beam height
PriorOpts.Marginals(2).Type = 'Constant';
PriorOpts.Marginals(2).Parameters = [0.3];       % (m)

PriorOpts.Marginals(3).Name = 'L';               % beam length
PriorOpts.Marginals(3).Type = 'Constant';
PriorOpts.Marginals(3).Parameters = 5;           % (m)

PriorOpts.Marginals(4).Name = 'E';               % Young's modulus
PriorOpts.Marginals(4).Type = 'LogNormal';
PriorOpts.Marginals(4).Moments = [30000 4500];   % (MPa)

PriorOpts.Marginals(5).Name = 'p';               % uniform load
PriorOpts.Marginals(5).Type = 'Gaussian';
PriorOpts.Marginals(5).Moments = [0.012 0.012*0.05]; % (kN/m)

myPriorDist = uq_createInput(PriorOpts);

%% 4 - MEASUREMENT DATA
% The measurement data consists of $N = 5$ independent measurements of
% the beam mid-span deflection.
% The data is stored in the column vector |y|:
myData.y = [12.84; 13.12; 12.13; 12.19; 12.67]/1000;  % (m)
myData.Name = 'Mid-span deflection';

%% 5 - BAYESIAN ANALYSIS
% The options of the Bayesian inversion analysis are specified with
% the following structure:
BayesOpts.Type = 'Inversion';
BayesOpts.Data = myData;
BayesOpts.Solver.Type = 'MCMC';
BayesOpts.Solver.MCMC.Sampler = 'MH';
BayesOpts.Solver.MCMC.Steps = 1e3;
BayesOpts.Solver.MCMC.NChains = 1e2;
BayesOpts.Solver.MCMC.T0 = 1e1;
BayesOpts.Solver.MCMC.a = 2;

Full = 1; % Full model
Surrogate = 0; % Surrogate model

if (Full==1)
    BayesOpts.ForwardModel.Model = myForwardModel; % Full model
    myBayesianAnalysis_fullModel = uq_createAnalysis(BayesOpts);
    uq_print(myBayesianAnalysis_fullModel)
    uq_display(myBayesianAnalysis_fullModel)
end


if (Surrogate==1)
    
    MetaOpts.Type = 'Metamodel';
    MetaOpts.MetaType = 'PCE'; % Polynomial chaos expansion
    MetaOpts.Method = 'LARS'; % Least-angle regression
    
    MetaOpts.ExpDesign.Sampling = 'LHS'; % Latin hypercube sampling
    MetaOpts.ExpDesign.NSamples = 50; % Number of samples
    MetaOpts.Degree = 1:4; % Polynomial degree
    
    MetaOpts.Input = myPriorDist;
    MetaOpts.FullModel = myForwardModel;
    mySurrogateModel = uq_createModel(MetaOpts);
    
    BayesOpts.ForwardModel.Model = mySurrogateModel;
    myBayesianAnalysis_surrogateModel = uq_createAnalysis(BayesOpts);
    uq_print(myBayesianAnalysis_surrogateModel)
    uq_display(myBayesianAnalysis_surrogateModel)
end

