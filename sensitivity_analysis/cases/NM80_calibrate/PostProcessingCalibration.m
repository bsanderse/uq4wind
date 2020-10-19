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
% multiply by the scaling value set in NM80_calibrate.m, which is stored in
% P{26}
for i=1:ndim
    Delta(i) = BayesianAnalysis.Results.PostProc.PointEstimate.X(i)*P{26}{i}{3};
end
% Delta_1 = BayesianAnalysis.Results.PostProc.PointEstimate.X(1)*0.2;
% Delta_2 = BayesianAnalysis.Results.PostProc.PointEstimate.X(2)*0.2;
% Delta_3 = BayesianAnalysis.Results.PostProc.PointEstimate.X(3)*0.3;
% Delta_4 = BayesianAnalysis.Results.PostProc.PointEstimate.X(4)*0.3;


%% Write calibrated polars using mean of posterior
run('write_calibration.m');
% run('../../AEROmodule/NM80_calibrate/write_calibrated_polars.m');

% %% Cross-validation (Optional)
% figure()
% R = [13,19,30,37];
% exp = [474.7 817.9 1210.8 1254.3];
% aero = [0.5378    0.8691    1.3804    1.4778]*10^3;
% calibrated = [0.4747    0.8179    1.2108    1.2543]*10^3;
% plot(R, exp,'g-*', 'LineWidth', 1)
% hold on
% plot(R, aero, 'r-o', 'LineWidth', 1)
% hold on
% plot(R, calibrated, 'b-o', 'LineWidth', 1)
% xlabel('R [m]')
% ylabel('\mu_{n} [N/m]')
% grid on
% legend('Experimental', 'Non-calibrated AeroModule', 'Calibrated AeroModule','Location','southeast')


