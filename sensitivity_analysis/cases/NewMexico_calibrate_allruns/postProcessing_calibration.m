% ================================== PLOTTING and POSTPROCESSING =========================
% Print out a report of the results:
uq_print(BayesianAnalysis)
uq_display(BayesianAnalysis)

if (exist('test_run','var') && test_run == 1)
    f = gcf;
    k = f.Number - n_runs;
    for j=1:n_runs
        figure(k+j)
        hold on
        plot(1:length(Y_unpert(j,:)),Y_unpert(j,:),'s');
    end
end
% data:
% plot(1:length(Data.y),Data.y,'s')

uq_display(BayesianAnalysis, 'meanConvergence', 'all')
uq_display(BayesianAnalysis, 'trace', 'all')
%uq_display(BayesianAnalysis, 'acceptance', 'true')

% check convergence of MCMC
uq_postProcessInversion(BayesianAnalysis,'gelmanRubin', 'true')
R_hat_full = BayesianAnalysis.Results.PostProc.MPSRF;

if (R_hat_full <= 2)
    disp('The MCMC simulation has converged')
else
    disp('The MCMC simulation has not converged. Increase the number of samples or fine tune the algorithm.')
end

% store the MAP into myBayesianAnalysis.Results.PostProc.PointEstimate:
uq_postProcessInversion(BayesianAnalysis,'pointEstimate', 'MAP')

