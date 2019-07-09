%% Sobol indices computed with Monte Carlo and PCE-type methods

clc
close all
clearvars

input_file = 'airfoil_lift'; % specify directory which contains test case settings and model


%% Sobol options
SobolOpts.Type        = 'Sensitivity';
SobolOpts.Method      = 'Sobol';
SobolOpts.Sobol.Order = 1;


%% initialize UQlab
% add path
addpath(genpath('../../UQLabCore_Rel1.0.0/'));
% start uqlab
uqlab


%% start running
run(['cases/' input_file '/initialize.m']);

if (exist('mean_exact','var'))
    compare_mean = 1;
else
    compare_mean = 0;
end
if (exist('std_exact','var'))
    compare_std = 1;
else
    compare_std = 0;
end

myModel = uq_createModel(Model);      % create and add the model to UQLab
myInput = uq_createInput(Input) ;
uq_print(myInput);

uq_display(myInput);


%% MC
if (find(strcmp(methods,'MC')))
    
    disp('=========MC==========');
    
    for k = 1:MC_repeat
        for i=1:length(NsamplesMC)
            disp([NsamplesMC(i)]);
            X_ED2 = uq_getSample(NsamplesMC(i), 'MC') ;
            Y_ED2 = uq_evalModel(myModel,X_ED2) ;
            mean_MC(k,i) = mean(Y_ED2);
            std_MC(k,i)  = std(Y_ED2);
            
            
            SobolOpts.Sobol.SampleSize = NsamplesMC(i);
            SobolAnalysis_MC = uq_createAnalysis(SobolOpts);
            SobolResults_MC  = SobolAnalysis_MC.Results;
            Sobol_MC(k,i,1:ndim) = SobolResults_MC.FirstOrder;
            %             Sobol_MC(k,i,1:ndim) = SobolResults_MC.Total;
            
            
        end
    end
    % take average over first dimension
    Sobol_MC = squeeze(mean(Sobol_MC,1));
    if (compare_mean == 1)
        err_mean_MC = abs((mean(mean_MC,1)-mean_exact)/mean_ref);
    end
    if (compare_std == 1)
        err_std_MC = abs((mean(std_MC,1)-std_exact)/std_ref);
    end
    
end

%% PCE

if (find(strcmp(methods,'PCE_Quad')))
    
    disp('=========PCE==========');
    
    metamodelQuad.FullModel = myModel;
    metamodelQuad.Input     = myInput;
    metamodelQuad.Type      = 'Metamodel';
    metamodelQuad.MetaType  = 'PCE';
    
    metamodelQuad.Method          = 'Quadrature' ;
    metamodelQuad.Quadrature.Type = 'Full';
    
    for i = 1:length(DegreesPCE)
        
        metamodelQuad.Degree = DegreesPCE(i);
        myPCE_Quad           = uq_createModel(metamodelQuad);
        
        SobolAnalysis_PCE    = uq_createAnalysis(SobolOpts);
        SobolResults_PCE     = SobolAnalysis_PCE.Results;
        Sobol_PCE(i,1:ndim)  = SobolResults_PCE.FirstOrder;
        %         Sobol_PCE(i,1:ndim)  = SobolResults_PCE.Total;
        
        NsamplesPCE(i) = myPCE_Quad.ExpDesign.NSamples;
        mean_PCE(i)    = myPCE_Quad.PCE.Moments.Mean;
        std_PCE(i)     = sqrt(myPCE_Quad.PCE.Moments.Var);
        
    end
    
    if (compare_mean == 1)
        err_mean_PCE =  abs( (mean_PCE-mean_exact)/mean_ref);
    end
    if (compare_std == 1)
        err_std_PCE  =  abs( (std_PCE - std_exact)/std_ref);
    end
    
end

%% OLS
if (find(strcmp(methods,'PCE_OLS')))
    
    disp('=========OLS==========');
    
    metamodelOLS.FullModel = myModel;
    metamodelOLS.Input = myInput;
    metamodelOLS.Type = 'Metamodel';
    metamodelOLS.MetaType = 'PCE';
    metamodelOLS.Method = 'OLS' ;
    metamodelOLS.Degree = 1:6;
    
    % use this if issues with LOO: metamodelOLS.OLS.ModifiedLOO = 0;
    % note default sampling: LHS
    
    NsamplesOLS = NsamplesPCE;
    
    for k = 1:MC_repeat
        for i = 1:length(NsamplesOLS)
            
            metamodelOLS.ExpDesign.NSamples = NsamplesOLS(i);
            metamodelOLS.ExpDesign.Sampling = 'LHS'; % LHS is default
            myPCE_OLS = uq_createModel(metamodelOLS);
            mean_OLS(k,i) = myPCE_OLS.PCE.Moments.Mean;
            std_OLS(k,i) = sqrt(myPCE_OLS.PCE.Moments.Var);
        end
        
    end
    
    % take mean over first dimension (k)
    if (compare_mean == 1)
        err_mean_OLS =  abs((mean(mean_OLS,1)-mean_exact)/mean_ref);
    end
    if (compare_std == 1)
        err_std_OLS =  abs((mean(std_OLS,1)-std_exact)/std_ref);
    end
    
end

%% PCE LARS
if (find(strcmp(methods,'PCE_LARS')))
    
    disp('=========LARS==========');
    
    metamodelLARS.FullModel = myModel;
    metamodelLARS.Input     = myInput;
    metamodelLARS.Type      = 'Metamodel';
    metamodelLARS.MetaType  = 'PCE';
    metamodelLARS.Method    = 'LARS';
    metamodelLARS.Degree    = 1:6; % this automatically switches on degree adaptive PCE
    metamodelLARS.TruncOptions.qNorm = 0.75;
    
    % NsamplesLARS = [4 25 50 100 200];
    NsamplesLARS = NsamplesPCE;
    
    for k = 1:MC_repeat
        for i = 1:length(NsamplesLARS)
            
            % use manual experimental design:
            %         X_ED = uq_getSample(NsamplesLARS(i),'MC') ;
            %         Y_ED = uq_evalModel(myModel,X_ED);
            %         metamodelLARS.ExpDesign.X = X_ED;
            %         metamodelLARS.ExpDesign.Y = Y_ED;
            
            % use sampling strategy, note default MC!
            metamodelLARS.ExpDesign.Sampling = 'LHS'; % or 'LHS' or 'Sobol' or 'Halton'
            metamodelLARS.ExpDesign.NSamples = NsamplesLARS(i);
            
            myPCE_LARS     = uq_createModel(metamodelLARS);
            
            mean_LARS(k,i) = myPCE_LARS.PCE.Moments.Mean;
            std_LARS(k,i)  = sqrt(myPCE_LARS.PCE.Moments.Var);
        end
    end
    
    if (compare_mean == 1)
        err_mean_LARS =  abs((mean(mean_LARS,1)-mean_exact)/mean_ref);
    end
    if (compare_std == 1)
        err_std_LARS =  abs((mean(std_LARS,1)-std_exact)/std_ref);
    end
    
end

%% postprocessing
% uq_figure('Position', [50 50 500 400]);
% myColors = uq_cmap(2);

%% mean
if (compare_mean == 1)
    
    figure
    if (find(strcmp(methods,'MC')))
        loglog(NsamplesMC, err_mean_MC, 'x-','Linewidth', 2); %, 'Color', myColors(1,:));
        hold on
    end
    if (find(strcmp(methods,'PCE_Quad')))
        loglog(NsamplesPCE, err_mean_PCE, 's-','Linewidth', 2);%,'Color', myColors(2,:));
        hold on
    end
    if (find(strcmp(methods,'PCE_OLS')))
        loglog(NsamplesOLS, err_mean_OLS, 'o-','Linewidth', 2); %,'Color', 'y');
        hold on
    end
    if (find(strcmp(methods,'PCE_LARS')))
        loglog(NsamplesLARS, err_mean_LARS, 'd-','Linewidth', 2); %,'Color', 'r');
        hold on
    end
    xlabel('N') % Add proper labelling and a legend
    legend({'Monte Carlo','PCE - quadrature','PCE - least squares','PCE - LARS'})
    grid on;
    title('Error in mean')
end

%% standard deviation
if (compare_std == 1)
    
    figure
    if (find(strcmp(methods,'MC')))
        loglog(NsamplesMC, err_std_MC, 'x-','Linewidth', 2); %, 'Color', myColors(1,:));
        hold on
    end
    if (find(strcmp(methods,'PCE_Quad')))
        loglog(NsamplesPCE, err_std_PCE, 's-','Linewidth', 2);%,'Color', myColors(2,:));
        hold on
    end
    if (find(strcmp(methods,'PCE_OLS')))
        loglog(NsamplesOLS, err_std_OLS, 'o-','Linewidth', 2); %,'Color', 'y');
        hold on
    end
    if (find(strcmp(methods,'PCE_LARS')))
        loglog(NsamplesLARS, err_std_LARS, 'd-','Linewidth', 2); %,'Color', 'r');
        hold on
    end
    xlabel('N') % Add proper labelling and a legend
    legend({'Monte Carlo','PCE - quadrature','PCE - least squares','PCE - LARS'})
    grid on;
    title('Error in standard deviation')
    
end

%% Sobol indices

figure
cmap = get(gca,'ColorOrder');

if (find(strcmp(methods,'MC')))
    semilogx(NsamplesMC', Sobol_MC, 'x-','Linewidth', 2, 'Color', cmap(1,:),'EdgeColor', 'none');
    hold on
end
if (find(strcmp(methods,'PCE_Quad')))
    semilogx(NsamplesPCE', Sobol_PCE, 's-','Linewidth', 2,'Color', cmap(2,:));
    hold on
end
if (find(strcmp(methods,'PCE_OLS')))
    loglog(NsamplesOLS, Sobol_OLS, 'o-','Linewidth', 2,'Color', cmap(3,:));
    hold on
end
if (find(strcmp(methods,'PCE_LARS')))
    loglog(NsamplesLARS, Sobol_LARS, 'd-','Linewidth', 2,'Color', cmap(4,:));
    hold on
end
xlabel('N') % Add proper labelling and a legend
% legend({'Monte Carlo','PCE - quadrature','PCE - least squares','PCE - LARS'})
ylabel('Total index');
grid on;
title('Comparison of Sobol indices')

%% 
figure
uq_bar((1:ndim)-0.125, SobolResults_MC.FirstOrder, 0.25, 'FaceColor', cmap(1,:), 'EdgeColor', 'none')
hold on
uq_bar((1:ndim)+0.125, SobolResults_PCE.FirstOrder, 0.25, 'FaceColor', cmap(2,:), 'EdgeColor', 'none')
hold on
% uq_bar((1:ndim)+0.25, mySobolResultsLRA.Total, 0.25,...
%     'FaceColor', cm(64,:), 'EdgeColor', 'none')
% uq_setInterpreters(gca)
set(gca, 'XTick', 1:length(Input.Marginals),...
    'XTickLabel', SobolResults_PCE.VariableNames, 'FontSize', 14)
uq_legend({...
    sprintf('MC based (%.0e simulations)', NsamplesMC(end)),...
    sprintf('PCE-based (%d simulations)', myPCE_Quad.ExpDesign.NSamples)})
ylabel('Total order Sobol index');
ylim([0 1])


%% plot polynomial approximations for quadrature-based methods
if (find(strcmp(methods,'PCE_Quad')))
    
    n_inputs = length(Input.Marginals);
    
    if (n_inputs == 1)
        Xsamples  = myPCE_Quad.ExpDesign.X;
        Ysamples  = myPCE_Quad.ExpDesign.Y;
        
        N = 100;
        X = linspace(Input.Marginals(1).Bounds(1),Input.Marginals(1).Bounds(2),N)';
        
        Y     = uq_evalModel(myModel,X);
        Y_PCE = uq_evalModel(myPCE_Quad,X); % or myPCE_OLS, myPCE_LARS
        
        figure
        plot(X,Y,'-');
        hold on
        plot(X,Y_PCE,'--');
        plot(Xsamples,Ysamples,'s');
        
        
    elseif (n_inputs >= 2)
        
        Xsamples  = myPCE_Quad.ExpDesign.X;
        Ysamples  = myPCE_Quad.ExpDesign.Y;
        
        % plot two parameters
        N1 = 100;
        N2 = 100;
        
        % select which parameters to plot
        p1 = 1;
        p2 = 2;
        p3 = 3; % we will use mean for p3
        
        % create regular grid of points where surrogate model will be evaluated
        % (this should be cheap)
        if (strcmp(myInput.Marginals(p1).Type,'Gaussian'))
            X1 = linspace(myInput.Marginals(p1).Moments(1)-2*myInput.Marginals(p1).Moments(2),...
                myInput.Marginals(p1).Moments(1)+2*myInput.Marginals(p1).Moments(2),N1)';
        elseif (strcmp(myInput.Marginals(p1).Type,'Uniform'))
            X1 = linspace(myInput.Marginals(p1).Bounds(1),myInput.Marginals(p1).Bounds(2),N1)';
        elseif (strcmp(myInput.Marginals(p1).Type,'Weibull'))
            X1 = linspace(0,myInput.Marginals(p1).Moments(1) + 3*myInput.Marginals(p1).Moments(2),N1)';
        end
        if (strcmp(myInput.Marginals(p2).Type,'Gaussian'))
            X2 = linspace(myInput.Marginals(p2).Moments(1)-2*myInput.Marginals(p2).Moments(2),...
                myInput.Marginals(p2).Moments(1)+2*myInput.Marginals(p2).Moments(2),N2)';
        elseif (strcmp(myInput.Marginals(p2).Type,'Uniform'))
            X2 = linspace(myInput.Marginals(p2).Bounds(1),myInput.Marginals(p2).Bounds(2),N2)';
        elseif (strcmp(myInput.Marginals(p2).Type,'Weibull'))
            X2 = linspace(0,myInput.Marginals(p2).Moments(1) + 3*myInput.Marginals(p2).Moments(2),N2)';            
        end
        
        if (n_inputs==2)
            X1new=kron(X1,ones(N2,1));
            X2new=kron(ones(N1,1),X2);
            X = [X1new X2new];
            plotting = 1;
        elseif (n_inputs==3)
            % for n>2, it is difficult to plot something in 3D space
            % for quadrature methods based on tensor formulation,
            % number of samples in one direction = degree+1
            % take 'middle' (mean?) point of samples in X3
            samples1D = 1+DegreesPCE(end);
            start1D  = ceil(samples1D/2);
            X3    = Xsamples(start1D,p3);
            
            X1new = kron(X1,ones(N2,1));
            X2new = kron(ones(N1,1),X2);
            X3new = ones(size(X1new))*X3;
            X = [X1new X2new X3new];
            plotting = 1;
        else
            warning('plotting not correctly implemented for high dimensional spaces');
            plotting = 0;
        end
        
        if (plotting == 1)
            % evaluate original model (can be expensive)
            Y     = uq_evalModel(myModel,X);
            % evaluate surrogate model (should be cheap)
            Y_PCE = uq_evalModel(myPCE_Quad,X); % or myPCE_OLS, myPCE_LARS
            
            figure
            surf(X1,X2,reshape(Y_PCE,N2,N1));
            hold on
            %     surf(X1,X2,reshape(Y,N2,N1));
            if (n_inputs==2)
                plot3(Xsamples(:,1),Xsamples(:,2),Ysamples,'s','MarkerSize',16,'MarkerFaceColor','black');
            elseif (n_inputs==3)
                plot3(Xsamples(start1D:samples1D:end,1),Xsamples(start1D:samples1D:end,2),Ysamples(start1D:samples1D:end),'s','MarkerSize',16,'MarkerFaceColor','black');
            end
            
        end
        
    end
    
end

