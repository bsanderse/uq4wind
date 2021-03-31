%% test UQLab vector QoI behavior

clearvars
rng default
% rng('shuffle')

run('config.m');
addpath(genpath(UQLab_path));

uqlab;

%% model definition
Model.mString = '[sin(2*pi*X(:,1)) cos(2*pi*X(:,1))]';
% Model.mString = '[sin(2*pi*X(:,1))]';
myModel       = uq_createModel(Model);
Nout  = 2;

Nconst = 3;

%% list with number of samples to test
Samples_list = 2:5;
Nrun         = length(Samples_list);
error_LOO    = zeros(Nrun,1);

for k=1:Nrun
    
    %% input
    mu    = 0.;
    sigma = 0.25;  
    
%     Input.Marginals(1).Type = 'Uniform';
%     Input.Marginals(1).Parameters = [mu-2*sigma, mu+2*sigma];
%     Input.Marginals(1).Bounds = [mu-2*sigma, mu+2*sigma];
    Input.Marginals(1).Type = 'Gaussian';
    Input.Marginals(1).Parameters = [mu,sigma];
    Input.Marginals(1).Bounds = [mu-2*sigma, mu+2*sigma];
    for c=1:Nconst
        Input.Marginals(c+1).Type = 'Constant';
        Input.Marginals(c+1).Parameters = 1;
    end
    
    myInput = uq_createInput(Input);
    
    %% LARS
    metamodelLARS.FullModel = myModel;
    metamodelLARS.Input     = myInput;
    metamodelLARS.Type      = 'Metamodel';
    metamodelLARS.MetaType  = 'PCE';
    metamodelLARS.Method    = 'LARS';
    metamodelLARS.Degree    = 1:8; % this automatically switches on degree adaptive PCE
    metamodelLARS.TruncOptions.qNorm = 0.5:0.1:1.5;
    
    metamodelLARS.ExpDesign.Sampling = 'LHS'; % or 'LHS' or 'Sobol' or 'Halton'
    metamodelLARS.ExpDesign.NSamples = Samples_list(k);
    
    myPCE_LARS  = uq_createModel(metamodelLARS);
    
    
    %% plot results
    Xsamples  = myPCE_LARS.ExpDesign.X;
    Ysamples  = myPCE_LARS.ExpDesign.Y;
    
    % evaluate the PCE model at many points
    X_PCE = linspace(Input.Marginals(1).Bounds(1),Input.Marginals(1).Bounds(2),100)';
    % enrich X_PCE with constants
    X_PCE = [X_PCE ones(length(X_PCE),Nconst)];
    Y_PCE = uq_evalModel(myPCE_LARS,X_PCE);

    figure
    col = lines;
    for q=1:Nout
    plot(Xsamples(:,1),Ysamples(:,q),'o','Color',col(q,:));
    hold on
    plot(X_PCE(:,1),Y_PCE(:,q),'Color',col(q,:));
    end
        
%     error_LOO(k) = myPCE_LARS.Error.LOO;
    
end
