%% INVERSION: CALIBRATION OF MULTIPLE FORWARD MODELS
%
% This example shows how multiple computational models can be calibrated 
% simultaneously. In this example, two different sets of deformation 
% measurements of the same specimen are used to calibrate its Young's  
% modulus $E$, by subjecting it to an uncertain distributed load $p$ and 
% to an uncertain point load $P$.
%
% The first set of measurements refers to the mid-span deflection of the
% simply supported specimen under the distributed load $p$, while the
% second set refers to its elongation under the constant load $P$.

clc
close all
clearvars

%% Initialize UQLab 
uqlab
rng(100)

%% Prior distribution (Input parameters)

% The prior information about the model parameters is gathered in a
% probabilistic model that includes both known and unknown parameters.

PriorOpts.Name = 'Model parameters prior';
PriorOpts.Marginals(1).Name = 'b';
PriorOpts.Marginals(1).Type = 'Constant';
PriorOpts.Marginals(1).Parameters = 0.15; % (m)

PriorOpts.Marginals(2).Name = 'h';
PriorOpts.Marginals(2).Type = 'Constant';
PriorOpts.Marginals(2).Parameters = 0.3; % (m)

PriorOpts.Marginals(3).Name = 'L';
PriorOpts.Marginals(3).Type = 'Constant';
PriorOpts.Marginals(3).Parameters = 5; % (m)

PriorOpts.Marginals(4).Name = 'E';
PriorOpts.Marginals(4).Type = 'Lognormal';
PriorOpts.Marginals(4).Moments = [30000 4500] ; % (MPa)

PriorOpts.Marginals(5).Name = 'p';
PriorOpts.Marginals(5).Type = 'Constant';
PriorOpts.Marginals(5).Parameters = 0.012; % (kN/m)

PriorOpts.Marginals(6).Name = 'P';
PriorOpts.Marginals(6).Type = 'Constant';
PriorOpts.Marginals(6).Parameters = 0.05; % (kN)

myPriorDist = uq_createInput(PriorOpts);


%% Forward models (Black-box models)

% This computation of the first forward model is carried out by the function
% |uq_SimplySupportedBeam.m| supplied with UQLab.
% In addition, specify a |PMap| vector to specify which prior model
% parameters are to be passed to the first forward model:

ModelOpts1.Name = 'Beam bending deflection';
ModelOpts1.mFile = 'uq_SimplySupportedBeam';
ForwardModels(1).Model = uq_createModel(ModelOpts1);
ForwardModels(1).PMap = [1 2 3 4 5];


% The second forward model in this case is directly defined with a string:

ModelOpts2.Name = 'Beam elongation';
ModelOpts2.mString = 'X(:,5).*X(:,3)./(X(:,1).*X(:,2).*X(:,4))';
ForwardModels(2).Model = uq_createModel(ModelOpts2);
ForwardModels(2).PMap = [1 2 3 4 6];


%% Problem figures
% The simply supported beam problem is shown in the following figure:

uq_figure
[I,~] = imread('SimplySupportedBeam.png');
image(I)
axis equal
set(gca, 'visible', 'off')

% The beam elongation problem is shown in the following figure:

uq_figure
[I,~] = imread('ElongationBeam.png');
image(I)
axis equal
set(gca, 'visible', 'off')

%% Experimental data
V_mid = [12.84; 13.12; 12.13; 12.19; 12.67]/1000; % (m)
myData.Name = 'Beam mid-span deflection';
myData.y = V_mid;

% In the case of multiple forward models, and therefore different types of
% data, it is necessary to define the full |MOMap| array to identify which
% model output needs to be compared with which data set:

U = [0.485; 0.466; 0.486]/1000; % (m)
% Data group 1
myData(1).y = V_mid;
myData(1).Name = 'Beam mid-span deflection';
myData(1).MOMap = [ 1;... % Model ID
1]; % Output ID

% Data group 2
myData(2).y = U;
myData(2).Name = 'Beam elongation';
myData(2).MOMap = [ 2;... % Model ID
1]; % Output ID



%% Discrepancy
% % Discrepancy group 1
DiscrepancyPriorOpts1.Name = 'Prior of sigma_1ˆ2';
DiscrepancyPriorOpts1.Marginals(1).Name = 'Sigma2';
DiscrepancyPriorOpts1.Marginals(1).Type = 'Uniform';
DiscrepancyPriorOpts1.Marginals(1).Parameters = [0, 1e-4];
myDiscrepancyPrior1 = uq_createInput(DiscrepancyPriorOpts1);
DiscrepancyOpts(1).Type = 'Gaussian';
DiscrepancyOpts(1).Prior = myDiscrepancyPrior1;

% Discrepancy group 2
DiscrepancyOpts(2).Type = 'Gaussian';
DiscrepancyOpts(2).Parameters = 1e-6; % known discrepancy variance


%% Perform Bayesian inverse
BayesOpts.Type = 'Inversion';
BayesOpts.Name = 'Bayesian multiple models';
BayesOpts.Prior = myPriorDist;
BayesOpts.ForwardModel = ForwardModels;
BayesOpts.Data = myData;
BayesOpts.Discrepancy = DiscrepancyOpts;
% The affine invariant ensemble
% algorithm is employed in this example, with |200| iterations
% and |100| parallel chains:
BayesOpts.Solver.Type = 'MCMC';
BayesOpts.Solver.MCMC.Sampler = 'AIES'; % AM, HMC, AIES
BayesOpts.Solver.MCMC.Steps = 200;
BayesOpts.Solver.MCMC.NChains = 100;

myBayesianAnalysis = uq_createAnalysis(BayesOpts);


%% Post-processing
uq_postProcessInversion(myBayesianAnalysis,'burnIn', 0.3);
uq_print(myBayesianAnalysis)
uq_display(myBayesianAnalysis)

