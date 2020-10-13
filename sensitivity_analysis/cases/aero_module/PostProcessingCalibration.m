% ================================== PLOTTING and POSTPROCESSING =========================
% Print out a report of the results:
uq_print(BayesianAnalysis)
uq_display(BayesianAnalysis)
%uq_display(BayesianAnalysis, 'meanConvergence', 'all')
uq_display(BayesianAnalysis, 'trace', 'all')
%uq_display(BayesianAnalysis, 'acceptance', 'true')

% check convergence of MCMC
uq_postProcessInversion(BayesianAnalysis,'gelmanRubin', 'true')
R_hat_full = BayesianAnalysis.Results.PostProc.MPSRF;

if R_hat_full <= 5
    disp('The MCMC simulation has converged')
else
    disp('The MCMC simulation has not converged. Increase the number of samples or fine tune the algorithm.')
end

% store the MAP into myBayesianAnalysis.Results.PostProc.PointEstimate:
uq_postProcessInversion(BayesianAnalysis,'pointEstimate', 'MAP')
% Delta interpretation using the MAP of the posterior
Delta_1 = BayesianAnalysis.Results.PostProc.PointEstimate.X(1)*0.2;
Delta_2 = BayesianAnalysis.Results.PostProc.PointEstimate.X(2)*0.2;
Delta_3 = BayesianAnalysis.Results.PostProc.PointEstimate.X(3)*0.3;
Delta_4 = BayesianAnalysis.Results.PostProc.PointEstimate.X(4)*0.3;

