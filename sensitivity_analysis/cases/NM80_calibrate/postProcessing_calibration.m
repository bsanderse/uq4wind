% ================================== PLOTTING and POSTPROCESSING =========================

% note: 
% the data in BayesianAnalysis.Results.PostPro are not updated but 
% (re)created during every call to uq_postProcessInversin using the 
% information in BayesianAnalysis.Results.Sample and in myBayesianAnalysis.Results.ForwardModel.

% UQLab defaults: 
% * burn-in is 50%, i.e. first 50% MCMC samples is discarded
% * 1000 samples of posterior predictive are drawn

% change this here to 50% burn in and 0 posterior predictive samples
uq_postProcessInversion(BayesianAnalysis,'burnIn',0.5,'posteriorPredictive',0)


%% Print out a report of the results:
uq_print(BayesianAnalysis)


%% Display marginals of posterior as scatter plot
uq_display(BayesianAnalysis,'scatterplot','all')
% selected parameters only:
% uq_display(BayesianAnalysis,'scatterplot',[1:4])

%% change fonts in scatterplot:
fontsize = 14;

ax = findall(gcf,'type','axes');
% set(ax,'TickLabelInterpreter','none');
set(ax,'TickLabelInterpreter','none','FontSize',fontsize);
%
xlabel(ax(1),'\theta_{E,4}','Interpreter','tex','FontSize',fontsize)
xlabel(ax(2),'\theta_{E,3}','Interpreter','tex','FontSize',fontsize)
xlabel(ax(3),'\theta_{E,2}','Interpreter','tex','FontSize',fontsize)
xlabel(ax(4),'\theta_{E,1}','Interpreter','tex','FontSize',fontsize)
xlabel(ax(5),'\Delta C_{l,4}','Interpreter','tex','FontSize',fontsize)
xlabel(ax(6),'\Delta C_{l,3}','Interpreter','tex','FontSize',fontsize)
xlabel(ax(7),'\Delta C_{l,2}','Interpreter','tex','FontSize',fontsize)
xlabel(ax(8),'\Delta C_{l,1}','Interpreter','tex','FontSize',fontsize)

ylabel(ax(8),'\theta_{E,4}','Interpreter','tex','FontSize',fontsize)
ylabel(ax(16),'\theta_{E,3}','Interpreter','tex','FontSize',fontsize)
ylabel(ax(24),'\theta_{E,2}','Interpreter','tex','FontSize',fontsize)
ylabel(ax(32),'\theta_{E,1}','Interpreter','tex','FontSize',fontsize)
ylabel(ax(40),'\Delta C_{l,4}','Interpreter','tex','FontSize',fontsize)
ylabel(ax(48),'\Delta C_{l,3}','Interpreter','tex','FontSize',fontsize)
ylabel(ax(56),'\Delta C_{l,2}','Interpreter','tex','FontSize',fontsize)
ylabel(ax(64),'\Delta C_{l,1}','Interpreter','tex','FontSize',fontsize)

% remove title
ax(65).Title.String='';
% change legend
ax(57).Legend.Interpreter='none';



%% trace plots
%uq_display(BayesianAnalysis, 'meanConvergence', 'all')
% plot trace plot of all parameters:
% uq_display(BayesianAnalysis, 'trace', 'all')
% trace plot of selected parameters:
uq_display(BayesianAnalysis, 'trace', [1;5])
% acceptance rate per chain
% uq_display(BayesianAnalysis, 'acceptance', 'true')

% plot convergence of mean of selected parameters:
% uq_display(BayesianAnalysis,'meanConvergence',[1;5]);



%% get Gelman-Rubin diagnostic 
uq_postProcessInversion(BayesianAnalysis,'gelmanRubin', 'true','burnIn',0.5,'posteriorPredictive',0); 
R_hat_full = BayesianAnalysis.Results.PostProc.MPSRF

if R_hat_full <= 2
    disp('The MCMC simulation has converged')
else
    disp('The MCMC simulation has not converged. Increase the number of samples or fine tune the algorithm.')
end

%% get the MAP
% store the MAP into BayesianAnalysis.Results.PostProc.PointEstimate:
uq_postProcessInversion(BayesianAnalysis,'pointEstimate', 'MAP','burnIn',0.5,'posteriorPredictive',0);
X_MAP = BayesianAnalysis.Results.PostProc.PointEstimate.X;


%% plot the posterior for a certain variable
% choose which variable to plot (give index):
X_choice  = 1;
% all posterior samples (except the ones disregarded with burn-in)
% posterior is an array of size Nsteps*Nvars*Nchains
posterior = BayesianAnalysis.Results.PostProc.PostSample;
% select the wanted variable, and reshape the rest into a 1D array
post_X_choice = reshape(squeeze(posterior(:,X_choice,:)),[],1);
% make density estimate:
[post_density, x_post_density] = ksdensity(post_X_choice);

% add the prior density:
% truncated gaussian:
x_prior = linspace(myPrior.Marginals(X_choice).Bounds(1),myPrior.Marginals(X_choice).Bounds(2),100);
[p,x] = truncnormpdf(x_prior,...
    myPrior.Marginals(X_choice).Moments(1),myPrior.Marginals(X_choice).Moments(2),...
    myPrior.Marginals(X_choice).Bounds(1),myPrior.Marginals(X_choice).Bounds(2));

figure
plot(x,p);
hold on
plot(x_post_density,post_density);


%% get posterior predictive
nPred = 500;
uq_postProcessInversion(BayesianAnalysis,'posteriorPredictive',nPred,'burnIn',0.5)
postPred = BayesianAnalysis.Results.PostProc.PostPred.model.postPredRuns;


%% plot uncalibrated, calibrated (MAP), and posterior predictive
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

% put experimental data into a vector for plotting
for i=1:4
    mean_exp_data(i) = mean(Data(i).y);   
end
% get unperturbed model result from surrogate if not available yet
if (test_run == 0)    
    warning('unperturbed AeroModule results are obtained from the surrogate model');
    if (Surrogate_model_type == 0)
        Y_unperturbed = uq_evalModel(loaded_surrogate_model.mySurrogateModel, X_unperturbed);  
    else
        Y_unperturbed = uq_evalModel(mySurrogateModel, X_unperturbed);  
    end     
end

figure
h1 = plot(r_exp_data,mean_exp_data,'x','markersize',16,'Linewidth', 2.5);
t  = lines; %get(gca,'ColorOrder');

hold on

grey = [0.5 0.5 0.5];
violin(postPred,'x',r_exp_data,'edgecolor',grey-0.1,'facecolor',grey+0.1,'medc','','mc','','plotlegend','','facealpha',0.5);

for i=1:4
    Ndata = length(Data(i).y);
    h4 = plot(r_exp_data(i)*ones(Ndata,1),Data(i).y,'x','markersize',6);
    set(h4,'color',t(1,:));
%     set(h4, 'color', get(h1, 'color')); % Use same color to fill in markers

end

h2 = plot(r_exp_data,Y_MAP,'o','markersize',16,'color',t(2,:),'Linewidth', 2.5);
% set(h2, 'markerfacecolor', get(h2, 'color')); % Use same color to fill in markers
h3 = plot(r_exp_data,Y_unperturbed,'s','markersize',16,'color',t(3,:),'Linewidth', 2.5);
% set(h3, 'markerfacecolor', get(h3, 'color')); % Use same color to fill in markers


% colors = get(gca,'colororder');



grid on
xlabel('r [m]');
ylabel('F_{N} [N/m]');
% legend('Experimental data','Calibrated AeroModule (MAP)','Uncalibrated AeroModule');
xlim([0 40])
legend([h1 h2 h3],{'Experimental data','Calibrated Aero-Module (MAP)','Uncalibrated Aero-Module'},'Location','northwest')
set(gca,'FontSize',14);
set(gcf,'Color',[1 1 1]);

% for saving as pdf:
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(gcf, '-dpdf', 'Sectional_force_calibrated_uncalibrate.pdf');

%% plot calibrated polars
P.FixedParameters.plot_polar = 1;
% run at unperturbed setting
writeAeroModuleInputReplacement(X_unperturbed(1:ndim),P);
% run at MAP
writeAeroModuleInputReplacement(X_MAP(1:ndim),P);

%% fix ranges
for i=1:ndim
    figure(100+i)
    xlim([-10 50])
    grid on 
    xlabel('alpha [deg]')
    ylabel('C_l [-]')
    legend('Uncalibrated','Calibrated (MAP)');
    box on
end

%% save surrogate model
if (Bayes_full == 0 && Surrogate_model_type == 1)

    disp('saving surrogate model');
    filename = ['PCE_LARS_N' num2str(MetaOpts.ExpDesign.NSamples) '_adapted_thickness_gaussian.mat'];
    filepath = 'StoredSurrogates/NM80_calibrate/';
    % prevent overwriting of file
    filename = avoidOverwrite(filename,filepath);
    save(fullfile(filepath,filename),'mySurrogateModel');
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


