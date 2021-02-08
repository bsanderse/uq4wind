% ================================== PLOTTING and POSTPROCESSING =========================
% Print out a report of the results:
uq_print(BayesianAnalysis)

%% Display posterior
uq_display(BayesianAnalysis,'scatterplot','all')
%
% if (exist('test_run','var') && test_run == 1)
%     f = gcf;
%     k = f.Number - n_runs;
%     for j=1:n_runs
%         figure(k+j)
%         hold on
%         plot(1:length(Y_unpert(j,:)),Y_unpert(j,:),'s');
%     end
% end
% data:
% plot(1:length(Data.y),Data.y,'s')

uq_display(BayesianAnalysis, 'meanConvergence', 'all')
uq_display(BayesianAnalysis, 'trace', 'all')
uq_display(BayesianAnalysis, 'acceptance', 'true')

% check convergence of MCMC
uq_postProcessInversion(BayesianAnalysis,'gelmanRubin', 'true','burnIn',0.5,'posteriorPredictive',0)
R_hat_full = BayesianAnalysis.Results.PostProc.MPSRF;

if (R_hat_full <= 2)
    disp('The MCMC simulation has converged')
else
    disp('The MCMC simulation has not converged. Increase the number of samples or fine tune the algorithm.')
end

% store the MAP into myBayesianAnalysis.Results.PostProc.PointEstimate:
uq_postProcessInversion(BayesianAnalysis,'pointEstimate', 'MAP','burnIn',0.5,'posteriorPredictive',0)
X_MAP = BayesianAnalysis.Results.PostProc.PointEstimate.X;

%%
n_unc = length(X_MAP);
X_MAP_full = [X_MAP X_unperturbed(1,n_unc+1:end)];


for i = 1:n_runs
    
    %evaluate model at MAP
    if (Bayes_full == 0) % surrogate model used
        if (Surrogate_model_type == 0)
            Y_MAP = uq_evalModel(loaded_surrogate_models.mySurrogateModels(i).Model, X_MAP_full);
        else
            Y_MAP = uq_evalModel(mySurrogateModels(i).Model, X_MAP_full);
        end
    else
        Y_MAP = uq_evalModel(myForwardModels(i).Model, X_MAP_full);
    end
    % select mean
    mean_MAP = Y_MAP(1:n_coeffs:end);
    
    % put experimental data into a vector for plotting
    mean_exp_data = Data(i).y(1:n_coeffs:end);
    
    % get unperturbed model result from surrogate if not available yet
    if (test_run == 0)
        warning('unperturbed AeroModule results are obtained from the surrogate model');
        if (Surrogate_model_type == 0)
            Y_unperturbed = uq_evalModel(loaded_surrogate_models.mySurrogateModels(i).Model, X_unperturbed);
            mean_unperturbed = Y_unperturbed(1:n_coeffs:end);
        else
            Y_unperturbed = uq_evalModel(mySurrogateModels(i).Model, X_unperturbed);
            mean_unperturbed = Y_unperturbed(1:n_coeffs:end);
        end
    elseif (test_run == 1)
        mean_unperturbed = Y_unperturbed(i,1:n_coeffs:end);
    end
    
    figure
    plot(r_exp_data,mean_exp_data,'x');
    hold on
    plot(r_exp_data,mean_MAP,'o');
    plot(r_exp_data,mean_unperturbed,'s');
    
    grid on
    xlabel('r [m]');
    ylabel('Fn [N/m]');
    legend('Experimental data','Calibrated AeroModule (MAP)','Uncalibrated AeroModule');
    title(['Model vs. data for run ' num2str(select_runs(i))]);
end
