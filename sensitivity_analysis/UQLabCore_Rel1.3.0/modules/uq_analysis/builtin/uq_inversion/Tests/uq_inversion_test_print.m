function [success] = uq_inversion_test_print(level)
% UQ_INVERSION_TEST_DISPLAY tests the UQ_PRINT_UQ_INVERSION function
%
%   See also: UQ_SELFTEST_UQ_INVERSION

%% START UQLAB
uqlab('-nosplash');

if nargin < 1
    level = 'normal';
end
fprintf(['\nRunning: |' level '| uq_inversion_test_func_AM...\n']);

%% PROBLEM SETUP
[y,~,PriorDist,ForwardModel,~] = uq_inversion_test_MCMC_setup;

%% SOLVER
Solver.MCMC.Sampler = 'AIES';
Solver.MCMC.Steps = 20;
Solver.MCMC.nChains = 3;

%% BAYESIAN MODEL
BayesOpt.Type = 'Inversion';
BayesOpt.Name = 'Bayesian model';
BayesOpt.Prior = PriorDist;
BayesOpt.ForwardModel = ForwardModel;
BayesOpt.Solver = Solver;
BayesOpt.Data.y = y;
BayesianAnalysis = uq_createAnalysis(BayesOpt);
%% TEST print
try
  uq_print(BayesianAnalysis)
  uq_postProcessInversion(BayesianAnalysis,'pointEstimate','map');
  uq_print(BayesianAnalysis)
  success = 1;
catch
  success = 0;
end