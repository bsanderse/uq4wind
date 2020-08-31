clc
close all
clearvars

%% Experimental data
% The measurements refers to the mid-span deflection of a
% simply supported beam under the distributed load $p$,
exp_meas = [12.84; 13.12; 12.13; 12.19; 12.67]/1000; % Experimental measurements
x_val = [1:5]'; % Number of counts


%% A frequentist approach!
beta_OLS = regress(exp_meas, x_val)
linear_fit = LinearModel.fit(x_val,exp_meas);
plot(linear_fit, 'marker', 'o', 'MarkerFaceColor', 'b')
title('Linear regression')
xlabel('Number of counts')
ylabel('Deflection (m)')

% Evaluating the goodness of the fit
disp(linear_fit)

% Find a fit for any value of x
% Let's find the value for x = 10
random_estimation =  predict(linear_fit, 10); 
disp('The predicted value is');
disp(random_estimation)


%% A Bayesian approach!

%% Initialize UQLab 
uqlab

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
PriorOpts.Marginals(4).Type = 'Gaussian';
PriorOpts.Marginals(4).Moments = [30000, 4500] ; % (MPa)

PriorOpts.Marginals(5).Name = 'p';
PriorOpts.Marginals(5).Type = 'Constant';
PriorOpts.Marginals(5).Parameters = 0.012; % (kN/m)

PriorOpts.Marginals(6).Name = 'P';
PriorOpts.Marginals(6).Type = 'Constant';
PriorOpts.Marginals(6).Parameters = 0.05; % (kN)

myPriorDist = uq_createInput(PriorOpts);


%% Forward models

ModelOpts1.Name = 'Beam bending deflection';
ModelOpts1.mFile = 'uq_SimplySupportedBeam';
ForwardModels(1).Model = uq_createModel(ModelOpts1);
ForwardModels(1).PMap = [1 2 3 4 5];


%% Calling the experimental data
myData.Name = 'Beam mid-span deflection';
myData.y = exp_meas;
myData(1).y = exp_meas;
myData(1).Name = 'Beam mid-span deflection';
myData(1).MOMap = [ 1;... % Model ID
1]; % Output ID

%% Discrepancy
DiscrepancyPriorOpts1.Name = 'Prior of sigma_1ˆ2';
DiscrepancyPriorOpts1.Marginals(1).Name = 'Sigma2';
DiscrepancyPriorOpts1.Marginals(1).Type = 'Uniform';
DiscrepancyPriorOpts1.Marginals(1).Parameters = [0, 1e-4];
myDiscrepancyPrior1 = uq_createInput(DiscrepancyPriorOpts1);
DiscrepancyOpts(1).Type = 'Gaussian';
DiscrepancyOpts(1).Prior = myDiscrepancyPrior1;


%% Perform Bayesian inverse
BayesOpts.Type = 'Inversion';
BayesOpts.Name = 'Bayesian multiple models';
BayesOpts.Prior = myPriorDist;
BayesOpts.ForwardModel = ForwardModels;
BayesOpts.Data = myData;
BayesOpts.Discrepancy = DiscrepancyOpts(1);
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
