clc
close all
clearvars

%% DANAERO experiment
N_data = 20;  % Number of available dataset
F_max = 0.01; % Maximum azimuthally averaged spanwise force (N)
% Normal force (experimental data) in Newtons
F_n = [0.27; 0.28; 0.29; 0.3; 0.31 ; 0.32; 0.33;  0.53; 0.54; 0.55; ...
        0.56; 0.57; 0.80; 0.81; 0.82; 0.83; 0.84; 0.86; 0.87; 0.88]*F_max;
%From BEM code
% Beta angle distribution
beta_r = [5; 10; 15; 17]*(pi/180); % in radians
beta = repelem(beta_r,[7 5 5 3]); % assuming constant across distribution
% Chord distribution
c_r = [0.9; 0.6; 0.4; 0.3];
c = repelem(c_r,[7 5 5 3]); % assuming constant across distribution
% Relative velocity distribution
v_r = [22; 35; 50; 56];
v = repelem(v_r,[7 5 5 3]); % assuming constant across distribution

% DANAERO plot (modified)
figure
plot(beta,F_n,'.')
xlabel('\beta [radians]')
ylabel('F_{n} [-]')
title('Modified experimental data')

%% Frequentist approach
format compact
format long
p = 2; % Number of parameters to calibrate
A = zeros(N_data,p); % size of design matrix, normally N_data>p
q_test = eye(p,p);
for i=1:p
    A(:,i) = dan_model(q_test(i,:),beta);
end

q_freq = regress(F_n,A);
F_n_OLS = dan_model(q_freq,beta);
%% Bayesian approach
% start uqlab
uqlab

% Model settings
ModelOpts.mFile = 'dan_model_vectorized';
% pass design matrix as parameter to the M-file
ModelOpts.Parameters   = A; 
ModelOpts.isVectorized = true;
myForwardModel = uq_createModel(ModelOpts);
 
% Prior
% Note that the range needs to be a good estimate
PriorOpts.Marginals(1).Name = 'dL*';
PriorOpts.Marginals(1).Type = 'Uniform';
PriorOpts.Marginals(1).Parameters = [0,4e-2];
% 
PriorOpts.Marginals(2).Name = 'dD*';
PriorOpts.Marginals(2).Type = 'Gaussian';
PriorOpts.Marginals(2).Parameters = [4e-4,1e-3];
% 
myPriorDist = uq_createInput(PriorOpts);


% Display input properties
uq_print(myPriorDist);
uq_display(myPriorDist);

% Discrepancy
DiscrepancyOpts.Type = 'Gaussian';
DiscrepancyOpts.Parameters = 1e-6;
BayesOpts.Discrepancy = DiscrepancyOpts;

% Experimental data
myData.y       = F_n'; % note: UQLab uses a row vector here
myData.Name    = 'measurement data';
BayesOpts.Data = myData;
BayesOpts.Type = 'Inversion';


% MCMC parameters
BayesOpts.Solver.Type = 'MCMC';
% MCMC algorithms available in UQLab
MH = 1; % Metropolis-Hastings
AM = 0; % Adaptive Metropolis
AIES = 0; % Affine invariant ensemble
HMC = 0; % Hamilton Monte Carlo  

if (MH==1)
    BayesOpts.Solver.MCMC.Sampler = 'MH';
    BayesOpts.Solver.MCMC.Steps = 1e4;
    BayesOpts.Solver.MCMC.NChains = 1e2;
    BayesOpts.Solver.MCMC.T0 = 1e1;
end

if (AM==1)
    BayesOpts.Solver.MCMC.Sampler = 'AM';
    BayesOpts.Solver.MCMC.Steps = 1e4;
    BayesOpts.Solver.MCMC.NChains = 1e2;
    BayesOpts.Solver.MCMC.T0 = 1e1;
    BayesOpts.Solver.MCMC.Epsilon = 1e-4;
end

if (AIES==1)
    BayesOpts.Solver.MCMC.Sampler = 'AIES';
    BayesOpts.Solver.MCMC.Steps = 1e3;
    BayesOpts.Solver.MCMC.NChains = 1e2;
    BayesOpts.Solver.MCMC.T0 = 1e2;
end

if (HMC==1)
    BayesOpts.Solver.MCMC.Sampler = 'HMC';
    BayesOpts.Solver.MCMC.LeapfrogSteps = 1e3;
    BayesOpts.Solver.MCMC.LeapfrogSize = 0.01;
    BayesOpts.Solver.MCMC.Mass = 1;
end

% Post-processing
myBayesianAnalysis = uq_createAnalysis(BayesOpts);
uq_print(myBayesianAnalysis)
uq_display(myBayesianAnalysis)
uq_display(myBayesianAnalysis, 'meanConvergence', 'all')
uq_display(myBayesianAnalysis, 'trace', 'all')
uq_display(myBayesianAnalysis, 'acceptance', 'true')

%Point value for Q[q1, q2]
uq_postProcessInversion(myBayesianAnalysis,'pointEstimate', 'MAP')
uq_postProcessInversion(myBayesianAnalysis,'gelmanRubin', 'true')
% get the Maximum a Posteriori value (MAP)
q_Bayes = myBayesianAnalysis.Results.PostProc.PointEstimate.X;
F_n_Bayes = dan_model(q_Bayes,beta);


%% Post-processing

cl_freq = (q_freq(1)./0.5*1.225.*v_r.*c_r)
cd_freq = (q_freq(2)./0.5*1.225.*v_r.*c_r)

cl_bayes = (myBayesianAnalysis.Results.PostProc.PointEstimate.X(1)./0.5*1.225.*v_r.*c_r) % MAP
cd_bayes = (myBayesianAnalysis.Results.PostProc.PointEstimate.X(2)./0.5*1.225.*v_r.*c_r) % MAP


figure
plot(beta, F_n, 'kx')
hold on
plot(beta,F_n_OLS, 'g--x')
hold on
plot(beta,F_n_Bayes, 'r--o')
xlabel('\beta [radians]')
ylabel('F_{n} [-]')
title('Experiment vs. Frequentist vs. Bayesian')
legend('measurement data','OLS','Bayesian');

