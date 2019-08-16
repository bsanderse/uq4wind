%% Sobol indices computed with Monte Carlo and PCE-type methods

clc
close all
clearvars

input_file = 'aero_module'; % specify directory which contains test case settings and model


%% Sobol options

SobolOpts.Type        = 'Sensitivity';
SobolOpts.Method      = 'Sobol';
SobolOpts.Sobol.Order = 1;


%% initialize UQlab

% add path
addpath(genpath('../../UQLabCore_Rel1.0.0/'));
% start uqlab
uqlab


%% process input files
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

% create and add the model to UQLab
myModel = uq_createModel(Model);
% create input object with UQLab
myInput = uq_createInput(Input) ;

% display input properties
uq_print(myInput);
uq_display(myInput);


%% Monte-Carlo
if (find(strcmp(methods,'MC')))
    
    disp('=========MC==========');
    N_MC     = length(NsamplesMC);
    mean_MC  = zeros(MC_repeat,N_MC);
    std_MC   = zeros(MC_repeat,N_MC);
    Sobol_MC_FirstOrder = zeros(MC_repeat,N_MC,ndim);
    Sobol_MC_Total      = zeros(MC_repeat,N_MC,ndim);
    
    % perform multiple runs to decrease effect of random sampling
    for k = 1:MC_repeat
        
        
        for i=1:N_MC
            
            disp(NsamplesMC(i));
            % get random samples ('experimental design')
            X_ED = uq_getSample(NsamplesMC(i), 'MC') ;
            % evaluate model at sample
            Y_ED = uq_evalModel(myModel,X_ED) ;
            
            % moments of solution
            mean_MC(k,i) = mean(Y_ED);
            std_MC(k,i)  = std(Y_ED);
            
            % Sobol analysis; 
            SobolOpts.Sobol.SampleSize = NsamplesMC(i);
            SobolAnalysis_MC           = uq_createAnalysis(SobolOpts);
            SobolResults_MC            = SobolAnalysis_MC.Results;
            Sobol_MC_FirstOrder(k,i,1:ndim) = SobolResults_MC.FirstOrder;
            Sobol_MC_Total(k,i,1:ndim)      = SobolResults_MC.Total;
            
        end
    end
    
    % take average over first dimension (multiple MC runs)
    Sobol_MC_FirstOrder = squeeze(mean(Sobol_MC_FirstOrder,1));
    Sobol_MC_Total      = squeeze(mean(Sobol_MC_Total,1));
    
    if (compare_mean == 1)
        err_mean_MC = abs((mean(mean_MC,1)-mean_exact)/mean_ref);
    end
    if (compare_std == 1)
        err_std_MC = abs((mean(std_MC,1)-std_exact)/std_ref);
    end
 
end
%% Polynomial Chaos with quadrature

if (find(strcmp(methods,'PCE_Quad')))
    
    disp('=========PCE==========');
    
    N_Quad      = length(DegreesQuad);
    NsamplesQuad = zeros(N_Quad,1);
    mean_Quad    = zeros(N_Quad,1);
    std_Quad     = zeros(N_Quad,1);
    Sobol_Quad_FirstOrder = zeros(N_Quad,ndim);
    Sobol_Quad_Total = zeros(N_Quad,ndim);
    
    % set up PCE metamodel
    metamodelQuad.FullModel = myModel;
    metamodelQuad.Input     = myInput;
    metamodelQuad.Type      = 'Metamodel';
    metamodelQuad.MetaType  = 'PCE';
    
    metamodelQuad.Method          = 'Quadrature' ;
    metamodelQuad.Quadrature.Type = 'Full';
    
    
    for i = 1:N_Quad
        
        metamodelQuad.Degree = DegreesQuad(i);
        myPCE_Quad           = uq_createModel(metamodelQuad);
        
        % moments of solution
        NsamplesQuad(i) = myPCE_Quad.ExpDesign.NSamples;
        mean_Quad(i)    = myPCE_Quad.PCE.Moments.Mean;
        std_Quad(i)     = sqrt(myPCE_Quad.PCE.Moments.Var);
        
        % Sobol analysis
        % note the same options structure SobolOpts can be re-used to create a new analysis on the PCE model        
        SobolAnalysis_Quad    = uq_createAnalysis(SobolOpts);
        SobolResults_Quad     = SobolAnalysis_Quad.Results;
        Sobol_Quad_FirstOrder(i,1:ndim) = SobolResults_Quad.FirstOrder;
        Sobol_Quad_Total(i,1:ndim)      = SobolResults_Quad.Total;
        
        
    end
    
    if (compare_mean == 1)
        err_mean_Quad =  abs((mean_Quad - mean_exact)/mean_ref);
    end
    if (compare_std == 1)
        err_std_Quad  =  abs((std_Quad - std_exact)/std_ref);
    end
    
end


%% Polynomial Chaos with ordinary least squares (OLS)

if (find(strcmp(methods,'PCE_OLS')))
    
    disp('=========OLS==========');
    
    if (~exist('NsamplesOLS','var'))
        if (find(strcmp(methods,'PCE_Quad')))
            % we specify number of OLS samples based on samples used with
            % quadrature
            NsamplesOLS = NsamplesQuad;
            warning('number of OLS samples taken based on number of Quadrature samples');
        else
            error('please specify NsamplesOLS');
        end
    end
    
    
    N_OLS       = length(NsamplesOLS);
    mean_OLS    = zeros(N_OLS,1);
    std_OLS     = zeros(N_OLS,1);
    Sobol_OLS_FirstOrder = zeros(N_OLS,ndim);
    Sobol_OLS_Total      = zeros(N_OLS,ndim);
    
    metamodelOLS.FullModel = myModel;
    metamodelOLS.Input     = myInput;
    metamodelOLS.Type      = 'Metamodel';
    metamodelOLS.MetaType  = 'PCE';
    metamodelOLS.Method    = 'OLS';
    % specify array of possible degrees;
    % the degree with the lowest Leave-One-Out cross-validation error (LOO error)
    % is automatically selected:
    metamodelOLS.Degree    = 1:6;
    
    % if there are issues with LOO, try the following: metamodelOLS.OLS.ModifiedLOO = 0;
    % note that default sampling is LHS, this can be changed (see below)
    
    % as there is randomness in the experimental design, we can 
    % average over several runs
    for k = 1:OLS_repeat
        for i = 1:N_OLS
            
            metamodelOLS.ExpDesign.NSamples = NsamplesOLS(i);
            metamodelOLS.ExpDesign.Sampling = 'LHS'; % LHS is default
            myPCE_OLS = uq_createModel(metamodelOLS);
            
            % moments of solution            
            mean_OLS(k,i) = myPCE_OLS.PCE.Moments.Mean;
            std_OLS(k,i)  = sqrt(myPCE_OLS.PCE.Moments.Var);
            
            % Sobol analysis
            % note the same options structure SobolOpts can be re-used to create a new analysis on the PCE model
            SobolAnalysis_OLS    = uq_createAnalysis(SobolOpts);
            SobolResults_OLS     = SobolAnalysis_OLS.Results;
            Sobol_OLS_FirstOrder(i,1:ndim) = SobolResults_OLS.FirstOrder;
            Sobol_OLS_Total(i,1:ndim)      = SobolResults_OLS.Total;
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

%% Polynomial Chaos with LARS
if (find(strcmp(methods,'PCE_LARS')))
    
    disp('=========LARS==========');
    
    if (~exist('NsamplesLARS','var'))
        if (find(strcmp(methods,'PCE_Quad')))
            % we specify number of OLS samples based on samples used with
            % quadrature
            NsamplesLARS = NsamplesQuad;
            warning('number of LARS samples taken based on number of Quadrature samples');
        else
            error('please specify NsamplesLARS');
        end
    end
    
    N_LARS       = length(NsamplesLARS);
    mean_LARS    = zeros(N_LARS,1);
    std_LARS     = zeros(N_LARS,1);
    Sobol_LARS_FirstOrder = zeros(N_LARS,ndim);
    Sobol_LARS_Total      = zeros(N_LARS,ndim);
    
    metamodelLARS.FullModel = myModel;
    metamodelLARS.Input     = myInput;
    metamodelLARS.Type      = 'Metamodel';
    metamodelLARS.MetaType  = 'PCE';
    metamodelLARS.Method    = 'LARS';
    metamodelLARS.Degree    = 1:6; % this automatically switches on degree adaptive PCE
    metamodelLARS.TruncOptions.qNorm = 0.75;
    
    % as there is randomness in the experimental design, we can 
    % average over several runs    
    for k = 1:LARS_repeat
        for i = 1:N_LARS
            
            % use manual experimental design:
            %         X_ED = uq_getSample(NsamplesLARS(i),'MC') ;
            %         Y_ED = uq_evalModel(myModel,X_ED);
            %         metamodelLARS.ExpDesign.X = X_ED;
            %         metamodelLARS.ExpDesign.Y = Y_ED;
            
            % use sampling strategy, note that default is MC!
            metamodelLARS.ExpDesign.Sampling = 'LHS'; % or 'LHS' or 'Sobol' or 'Halton'
            metamodelLARS.ExpDesign.NSamples = NsamplesLARS(i);            
            myPCE_LARS     = uq_createModel(metamodelLARS);
            
            % moments of solution            
            mean_LARS(k,i) = myPCE_LARS.PCE.Moments.Mean;
            std_LARS(k,i)  = sqrt(myPCE_LARS.PCE.Moments.Var);

            % Sobol analysis
            % note the same options structure SobolOpts can be re-used to create a new analysis on the PCE model
            SobolAnalysis_LARS    = uq_createAnalysis(SobolOpts);
            SobolResults_LARS     = SobolAnalysis_LARS.Results;
            Sobol_LARS_FirstOrder(i,1:ndim) = SobolResults_LARS.FirstOrder;
            Sobol_LARS_Total(i,1:ndim)      = SobolResults_LARS.Total;
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
        loglog(NsamplesQuad, err_mean_Quad, 's-','Linewidth', 2);%,'Color', myColors(2,:));
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
        loglog(NsamplesQuad, err_std_Quad, 's-','Linewidth', 2);%,'Color', myColors(2,:));
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


%% Convergence of Sobol indices

figure
cmap = get(gca,'ColorOrder');

if (find(strcmp(methods,'MC')))
    semilogx(NsamplesMC', Sobol_MC_Total, 'x-','Linewidth', 2, 'Color', cmap(1,:));
    hold on
end
if (find(strcmp(methods,'PCE_Quad')))
    semilogx(NsamplesQuad', Sobol_Quad_Total, 's-','Linewidth', 2,'Color', cmap(2,:));
    hold on
end
if (find(strcmp(methods,'PCE_OLS')))
    loglog(NsamplesOLS, Sobol_OLS_Total, 'o-','Linewidth', 2,'Color', cmap(3,:));
    hold on
end
if (find(strcmp(methods,'PCE_LARS')))
    loglog(NsamplesLARS, Sobol_LARS_Total, 'd-','Linewidth', 2,'Color', cmap(4,:));
    hold on
end
xlabel('N') % Add proper labelling and a legend
% legend({'Monte Carlo','PCE - quadrature','PCE - least squares','PCE - LARS'})
ylabel('Total index');
grid on;
title('Comparison of Sobol indices')


%% bar chart of Sobol indices

figure
hold on

n_methods = length(methods);
bar_width = 0.5/n_methods;
bar_vec   = 1:n_methods;
coords    = (bar_vec - mean(bar_vec))*bar_width;
k         = 1;
if (find(strcmp(methods,'MC')))
    uq_bar((1:ndim)+coords(k), Sobol_MC_Total(end,:), bar_width, 'FaceColor', cmap(k,:), 'EdgeColor', 'none')
    k = k+1;
end
if (find(strcmp(methods,'PCE_Quad')))
    uq_bar((1:ndim)+coords(k), SobolResults_Quad.FirstOrder, bar_width, 'FaceColor', cmap(k,:), 'EdgeColor', 'none')
    k = k+1;
end
if (find(strcmp(methods,'PCE_OLS')))
    uq_bar((1:ndim)+coords(k), SobolResults_OLS.FirstOrder, bar_width, 'FaceColor', cmap(k,:), 'EdgeColor', 'none')
    k = k+1;
end
if (find(strcmp(methods,'PCE_LARS')))
    uq_bar((1:ndim)+coords(k), SobolResults_LARS.FirstOrder, bar_width, 'FaceColor', cmap(k,:), 'EdgeColor', 'none')
    k = k+1;
end

% uq_bar((1:ndim)+0.25, mySobolResultsLRA.Total, 0.25,...
%     'FaceColor', cm(64,:), 'EdgeColor', 'none')
% uq_setInterpreters(gca)
set(gca, 'XTick', 1:length(Input.Marginals),...
    'XTickLabel', SobolResults_Quad.VariableNames, 'FontSize', 14)
uq_legend({...
    sprintf('MC based (%.0e simulations)', NsamplesMC(end)),...
    sprintf('PCE-based (%d simulations)', myPCE_Quad.ExpDesign.NSamples)})
ylabel('Total order Sobol index');
ylim([0 1])



%% plot the polynomial response surface for PCE quadrature-based

if (find(strcmp(methods,'PCE_Quad')))
    
    n_inputs = length(Input.Marginals);
    
    if (n_inputs == 1)
        Xsamples  = myPCE_Quad.ExpDesign.X;
        Ysamples  = myPCE_Quad.ExpDesign.Y;
        
        N = 100;
        domain = getMarginalBounds(myInput.Marginals(1));
        X      = linspace(domain(1),domain(2),N)';
        
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
        domain1 = getMarginalBounds(myInput.Marginals(p1));
        X1      = linspace(domain1(1),domain1(2),N1)';
        domain2 = getMarginalBounds(myInput.Marginals(p2));
        X2      = linspace(domain2(1),domain2(2),N2)';
        
        
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
            samples1D = 1+DegreesQuad(end);
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
            title('polynomial response surface'); 
        end
        
    end
    
end

