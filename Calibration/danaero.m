clc
close all
clearvars

%% Data
% DANAERO data for normal aerodynamic force in non-sheared inflow 
% U_0 = 6.1 m/s
% https://backend.orbit.dtu.dk/ws/files/80542014/DanaeroFinalReport.pdf

R = 40; % Blade radius
N_data = 4;  % Number of available dataset
F_n = [ 0.3; 0.55; 0.82; 0.84]; % Normal force
F_n_pos = [ 0.015; 0.025; 0.032; 0.031]; % Standard deviation (+)
F_n_neg = [ 0.015; 0.025; 0.032; 0.031]; % Standard deviation (-)

% Plotting the experimental data
% Note: F_n is normalized
r = [ 13; 19; 30; 37]/R; % Measurements at blade locations
figure(1)
errorbar(r,F_n,F_n_neg,F_n_pos,'.')
xlim([0 1])
ylim([0 1.4])
xlabel('r/R [-]')
ylabel('F_{n} [-]')
title('Experimental data')

% From iterative BEM solution
% Beta values corresponding to the radial stations r (Assumption!)
beta = [66; 71; 80; 86]*(pi/180);  % in radians


figure(2)
errorbar(beta,F_n,F_n_neg,F_n_pos,'.')
xlim([0.9 1.6])
ylim([0 1.4])
xlabel('\beta [radians]')
ylabel('F_{n} [-]')
title('Modified experimental data')


%% Frequentist approach
p = 2; % Number of parameters to calibrate
A = zeros(N_data,p); % size of design matrix, normally N_data>p
q_test = eye(p,p);
for i=1:p
    A(:,i) = dan_model(q_test(i,:),beta);
end

q_freq = regress(F_n,A)
F_n_OLS = dan_model(q_freq,beta)

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
PriorOpts.Marginals(1).Parameters = [0.8,1.2];
% 
PriorOpts.Marginals(2).Name = 'dD*';
PriorOpts.Marginals(2).Type = 'Uniform';
PriorOpts.Marginals(2).Parameters = [-1.6,-0.5];
% 
myPriorDist = uq_createInput(PriorOpts);


% Display input properties
uq_print(myPriorDist);
uq_display(myPriorDist);

% Discrepancy
DiscrepancyOpts.Type = 'Gaussian';
DiscrepancyOpts.Parameters = 1e-3;
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
uq_display(myBayesianAnalysis, 'scatterplot', 'all')


%Point value for Q[q1, q2]
uq_postProcessInversion(myBayesianAnalysis,'pointEstimate', 'MAP')
uq_postProcessInversion(myBayesianAnalysis,'gelmanRubin', 'true')
% get the Maximum a Posteriori value (MAP)
q_Bayes = myBayesianAnalysis.Results.PostProc.PointEstimate.X;
F_n_Bayes = dan_model(q_Bayes,beta);


%%
figure
plot(beta, F_n, 'kx')
errorbar(beta,F_n,F_n_neg,F_n_pos,'.')
hold on
plot(beta,F_n_OLS, 'g--x')
hold on
plot(beta,F_n_Bayes, 'r--o')
xlim([1 1.6])
ylim([0 1.4])
xlabel('\beta [radians]')
ylabel('F_{n} [-]')
title('Experiment vs. Frequentist vs. Bayesian')
legend('measurement data','OLS','Bayesian');