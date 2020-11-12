% ================================== PLOTTING and POSTPROCESSING =========================

%% Print out a report of the results:
uq_print(BayesianAnalysis)


%% Display posterior
uq_display(BayesianAnalysis,'scatterplot','all')


%% trace plots
%uq_display(BayesianAnalysis, 'meanConvergence', 'all')
% plot trace plot of all parameters:
% uq_display(BayesianAnalysis, 'trace', 'all')
% trace plot of selected parameters:
uq_display(BayesianAnalysis, 'trace', [1;5])
%uq_display(BayesianAnalysis, 'acceptance', 'true')


%% post-process the MCMC results
% note: 
% the data in BayesianAnalysis.Results.PostPro are not updated but 
% (re)created during every call to uq_postProcessInversin using the 
% information in BayesianAnalysis.Results.Sample and in myBayesianAnalysis.Results.ForwardModel.

% UQLab defaults: 
% * burn-in is 50%, i.e. first 50% MCMC samples is discarded
% * 1000 samples of posterior predictive are drawn

% uq_postProcessInversion(BayesianAnalysis,'burnIn',0.5,'posteriorPredictive',0);
uq_postProcessInversion(BayesianAnalysis,'gelmanRubin', 'true', 'burnIn',0.5,'posteriorPredictive',0)
R_hat_full = BayesianAnalysis.Results.PostProc.MPSRF;

if R_hat_full <= 2
    disp('The MCMC simulation has converged')
else
    disp('The MCMC simulation has not converged. Increase the number of samples or fine tune the algorithm.')
end

% store the MAP into BayesianAnalysis.Results.PostProc.PointEstimate:
uq_postProcessInversion(BayesianAnalysis,'pointEstimate', 'MAP', 'burnIn',0.5,'posteriorPredictive',0)
X_MAP = BayesianAnalysis.Results.PostProc.PointEstimate.X;

%%
%evaluate model at MAP
if (Bayes_full == 0) % surrogate model used
    if (Surrogate_model_type == 0)
        Y_MAP = uq_evalModel(loaded_surrogate_model.mySurrogateModel, X_MAP(1:ndim));  
    else
        Y_MAP = uq_evalModel(mySurrogateModel, X_MAP(1:ndim));  
    end
else
    Y_MAP = uq_evalModel(myForwardModel, X_MAP(1:ndim));  
end

figure
hold on
for i=1:4
%      plot(r_exp_data(i)*ones(length(Data(i).y),1),Data(i).y);
    plot(r_exp_data(i),mean(Data(i).y),'x');
end
plot(r_exp_data,Y_MAP,'o');
if (test_run == 1)
    plot(r_exp_data,Y_test,'s');
end

%% Write calibrated polars using mean of posterior
% run('write_calibration.m');
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


