% ================================== PLOTTING and POSTPROCESSING =========================
% Print out a report of the results:
uq_print(BayesianAnalysis)

%% Display posterior (can be expensive to plot)
% uq_display(BayesianAnalysis,'scatterplot','all')
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

%% check convergence
% uq_display(BayesianAnalysis, 'meanConvergence', 'all')
 uq_display(BayesianAnalysis, 'trace', 'all')
%uq_display(BayesianAnalysis, 'acceptance', 'true')

% check convergence of MCMC
uq_postProcessInversion(BayesianAnalysis,'gelmanRubin', 'true','burnIn',0.5,'posteriorPredictive',0)
R_hat_full = BayesianAnalysis.Results.PostProc.MPSRF

if (R_hat_full <= 2)
    disp('The MCMC simulation has converged')
else
    disp('The MCMC simulation has not converged. Increase the number of samples or fine tune the algorithm.')
end

% store the MAP into myBayesianAnalysis.Results.PostProc.PointEstimate:
uq_postProcessInversion(BayesianAnalysis,'pointEstimate', 'MAP','burnIn',0.5,'posteriorPredictive',0)
% note that X_MAP contains both model parameters and the hyperparameters
% without the model constants (if present)
X_MAP = BayesianAnalysis.Results.PostProc.PointEstimate.X;
if (iscell(X_MAP)) %uqlab v1.4
 X_MAP = cell2mat(X_MAP);
end


%%
% number of model uncertainties, without hyperparameters
% as input to uq_eval, we need X without the hyperparameters, but with the
% constants
X_MAP_full = [X_MAP(1:nunc) X_unperturbed(1,nunc+1:end)];


for i = 1:n_runs
    
    % get posterior predictive
    nPred = 500;
    uq_postProcessInversion(BayesianAnalysis,'burnIn',0.5,'posteriorPredictive',nPred)

    if (isfield(BayesianAnalysis.Results.PostProc,'PostPred'))
        % uqlab v1.3
        postPred = BayesianAnalysis.Results.PostProc.PostPred.model.postPredRuns; 
    else
        % uqlab v1.4
        postPred = BayesianAnalysis.Results.PostProc.PostPredSample(i).PostPred;
    end


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
    
    % loop over selected Fourier coefficients
    for j=1:n_coeffs
        k = index_fourier(j);
        % select outputs at several sections corresponding to this Fourier coefficient
        % order is determined in NewMexico_readoutput
        QoI_MAP = Y_MAP(j:n_coeffs:end);

        % similarly put experimental data into a vector for plotting
        QoI_exp_data = Data(i).y(j:n_coeffs:end);

        % get unperturbed model result from surrogate if not available yet
        if (test_run == 0) % test run has not yet been performed
            warning('unperturbed AeroModule results are obtained from the surrogate model');
            if (Surrogate_model_type == 0)
                Y_unperturbed = uq_evalModel(loaded_surrogate_models.mySurrogateModels(i).Model, X_unperturbed);
                QoI_unperturbed = Y_unperturbed(j:n_coeffs:end);
            else
                Y_unperturbed = uq_evalModel(mySurrogateModels(i).Model, X_unperturbed);
                QoI_unperturbed = Y_unperturbed(j:n_coeffs:end);
            end
        elseif (test_run == 1)
            QoI_unperturbed = Y_unperturbed(i,j:n_coeffs:end);
        end

        figure
        h1 = plot(r_exp_data(r_index),QoI_exp_data,'x','markersize',16,'Linewidth', 2.5);
        t  = lines; %get(gca,'ColorOrder');
        
        hold on
        grey = [0.5 0.5 0.5];
        violin(postPred(:,j:n_coeffs:end),'x',r_exp_data(r_index),'edgecolor',grey-0.1,'facecolor',grey+0.1,'medc','','mc','','plotlegend','','facealpha',0.5);

        h2 = plot(r_exp_data(r_index),QoI_MAP,'o','markersize',16,'color',t(2,:),'Linewidth', 2.5);
        h3 = plot(r_exp_data(r_index),QoI_unperturbed,'s','markersize',16,'color',t(3,:),'Linewidth', 2.5);

        grid on
        xlabel('r [m]');
%         ylabel('Fn [N/m]');
        legend([h1 h2 h3],'Experimental data','Calibrated AeroModule (MAP)','Uncalibrated AeroModule');
        title(['Model vs. data for run ' num2str(select_runs(i)) ' and Fourier coefficient ' num2str(k)]);
    end
end

%% save surrogate model
if (Bayes_full == 0 && Surrogate_model_type == 1)

    filename = ['PCE_N' num2str(MetaOpts.ExpDesign.NSamples) '.mat'];
    filepath = fullfile('..','..','StoredSurrogates','NewMexico_calibrate');
    % prevent overwriting of file
    filename_new = avoidOverwrite(filename,filepath);
    disp(['saving surrogate model in ' filename_new]);

    save(fullfile(filepath,filename_new),'mySurrogateModels');
end
