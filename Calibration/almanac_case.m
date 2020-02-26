clc
close all
clearvars

%% Data

% Height and weight data from World Almanac Book Facts
N_data = 15;
height = linspace(58,72,N_data)';
weight = [115; 117; 120; 123; 126; 129; 132; 135; 139; 142; 146; 150; ...
          154; 159; 164];
p = 3; % number of terms to calibrate

%% Frequentist

A = zeros(N_data,p); % size of design matrix, normally N_data>p
q_test = eye(p,p);
for i=1:p
    A(:,i) = almanac(q_test(i,:),height);
end

q_freq = regress(weight,A);

weight_OLS = almanac(q_freq,height);


%% Bayesian
% start uqlab
uqlab

% Model settings
ModelOpts.mFile = 'almanac_vectorized';
% pass design matrix as parameter to the M-file
ModelOpts.Parameters   = A; 
ModelOpts.isVectorized = true;
myForwardModel = uq_createModel(ModelOpts);
 
% Prior
% Note that the range considered is large
PriorOpts.Marginals(1).Name = 'q1';
PriorOpts.Marginals(1).Type = 'Uniform';
PriorOpts.Marginals(1).Parameters = [200,300];
% 
PriorOpts.Marginals(2).Name = 'q2';
PriorOpts.Marginals(2).Type = 'Uniform';
PriorOpts.Marginals(2).Parameters = [-100,-10];
%
PriorOpts.Marginals(3).Name = 'q3';
PriorOpts.Marginals(3).Type = 'Uniform';
PriorOpts.Marginals(3).Parameters = [1,30];
% 
myPriorDist = uq_createInput(PriorOpts);


% Display input properties
uq_print(myPriorDist);
uq_display(myPriorDist);

% Discrepancy
DiscrepancyOpts.Type = 'Gaussian';
DiscrepancyOpts.Parameters = 1e-2;
BayesOpts.Discrepancy = DiscrepancyOpts;

% Experimental data
myData.y       = weight'; % note: UQLab uses a row vector here
myData.Name    = 'measurement data';
BayesOpts.Data = myData;
BayesOpts.Type = 'Inversion';

% MCMC parameters
BayesOpts.Solver.Type = 'MCMC';
BayesOpts.Solver.MCMC.Sampler = 'AM';
BayesOpts.Solver.MCMC.Steps = 1e3;
BayesOpts.Solver.MCMC.NChains = 1e2;
BayesOpts.Solver.MCMC.T0 = 1e2;


% Post-processing
myBayesianAnalysis = uq_createAnalysis(BayesOpts);
uq_print(myBayesianAnalysis)
uq_display(myBayesianAnalysis)
uq_display(myBayesianAnalysis, 'meanConvergence', 'all')
uq_display(myBayesianAnalysis, 'trace', 'all')
    
%Point value for Q[q1, q2, q3]
uq_postProcessInversion(myBayesianAnalysis,'pointEstimate', 'MAP')
 
% get the Maximum a Posteriori value (MAP)
q_Bayes = myBayesianAnalysis.Results.PostProc.PointEstimate.X;
weight_Bayes = almanac(q_Bayes,height);


%% Post-processing (Point estimates)

format bank % Rounding to 2 decimal places

% Frequentist

q_freq = regress(weight,A)

R_freq =  weight - weight_OLS; % Residual estimator 

sigma2_freq = (R_freq.'*R_freq)/(N_data-p)% Variance estimator

cov_freq = sigma2_freq*inv(A'*A)% Covariance matrix estimate

% 2 comes from z value for 0.02275
q_dev_freq = [q_freq - (2*sqrt(sigma2_freq))/sqrt(N_data), q_freq + ...
             (2*sqrt(sigma2_freq))/sqrt(N_data)]%(+-)two standard deviation 


% Bayesian

q_Bayes = myBayesianAnalysis.Results.PostProc.PointEstimate.X'

R_Bayes =  weight - weight_Bayes; % Residual estimator 

sigma2_Bayes = (R_Bayes.'*R_Bayes)/(N_data-p)% Variance estimator

cov_Bayes = sigma2_Bayes*inv(A'*A)% Covariance matrix estimate

% 2 comes from z value for 0.02275
q_dev_Bayes = [q_Bayes - (2*sqrt(sigma2_freq))/sqrt(N_data), q_Bayes + ...
             (2*sqrt(sigma2_freq))/sqrt(N_data)]%(+-) two standard deviation 
         

% Comapre the two approaches to the experimental data        
plot(height, weight, 'kx')
hold on
plot(height, weight_OLS,'b--')
hold on
plot(height,weight_Bayes,'r--o');
legend('measurement data','OLS','Bayesian');





