%% test UQLab LOO

%% example of issue:
% take LARS with 16 samples, rng default, and use degree 1:4. UQLab
% suggests a straight line (degree 1) as best fit

% maybe the issue is when we have periodic functions and the support of the
% pdf is larger than a typical period of the function

%% other issue:
% LHS sampling of 1D parameter is different from LHS sampling of 1D
% parameter plus constant -> reported to UQLab, see
% https://uqworld.org/t/surprising-results-with-lhs-sampling-and-constant-inputs/1080

clearvars
% rng(1000)
rng default
% rng('shuffle')

run('config.m');
addpath(genpath(UQLab_path));

uqlab;

%% model
% Model.mString = '3.87 - 0.03*sin(2*pi*(X(:,1)-100)/350)';
% Model.mString = 'sin(2*pi*(X(:,1)-50)/100)';
Model.mString = 'sin(2*pi*X(:,1))';%-0.5))';
myModel = uq_createModel(Model);

Samples_list = 5;
Nrun = length(Samples_list);
error_LOO = zeros(Nrun,1);

for k=1:Nrun
    
    %% input
    mu    = 0.;
    sigma = 0.25;  %0.25*350;
%     mu    = 296;  %0.5;
%     sigma = 120; %0.2;  %0.25*350;
%     Input.Marginals(1).Type = 'Gaussian';
%     Input.Marginals(1).Parameters = [mu, sigma];
%     Input.Marginals(1).Bounds = [mu-3*sigma mu+3*sigma];

    Input.Marginals(1).Type = 'Uniform';
    Input.Marginals(1).Parameters = [mu-2*sigma, mu+2*sigma];
    
    Input.Marginals(2).Type = 'Constant';
    Input.Marginals(2).Parameters = 1;
%     Input.Marginals(3).Type = 'Constant';
%     Input.Marginals(3).Parameters = 1;
%     Input.Marginals(4).Type = 'Constant';
%     Input.Marginals(4).Parameters = 1;
    
    myInput = uq_createInput(Input);
    
    %% LARS
    metamodelLARS.FullModel = myModel;
    metamodelLARS.Input     = myInput;
    metamodelLARS.Type      = 'Metamodel';
    metamodelLARS.MetaType  = 'PCE';
    metamodelLARS.Method    = 'LARS';
    metamodelLARS.Degree    = 1:4; % this automatically switches on degree adaptive PCE
%     metamodelLARS.TruncOptions.qNorm = 0.5:0.1:1.5;
    % metamodelLARS.LARS.ModifiedLOO = 0; % use standard LOO to decide on convergence
    
    metamodelLARS.ExpDesign.Sampling = 'LHS'; % or 'LHS' or 'Sobol' or 'Halton'
    metamodelLARS.ExpDesign.NSamples = Samples_list(k);
    
    myPCE_LARS  = uq_createModel(metamodelLARS);
    
    
    %% plot results
    Xsamples  = myPCE_LARS.ExpDesign.X;
    Ysamples  = myPCE_LARS.ExpDesign.Y;
    
    % evaluate the PCE model at many points
    bounds = getMarginalBounds(myInput.Marginals(1));
    X_PCE = linspace(bounds(1),bounds(2),100)';
%     X_PCE = linspace(Input.Marginals(1).Bounds(1),Input.Marginals(1).Bounds(2),100)';
    Y_PCE = uq_evalModel(myPCE_LARS,X_PCE);
    %
%     if  (Nrun==1)
        figure
        plot(Xsamples(:,1),Ysamples,'o');
        hold on
        plot(X_PCE(:,1),Y_PCE);
%     end
    error_LOO(k) = myPCE_LARS.Error.LOO;
    
    clear myInput myPCE_LARS
end

%%
if (Nrun>1)
    mean(error_LOO)
    std(error_LOO)
    figure
    semilogy(error_LOO,'s')
end
%%
% Xsamples =
%
%   119.8075
%   290.8403
%   419.0246
%   252.5712
%   503.9865
%   203.9508
%   313.7827
%   351.7113

% Input.Marginals
%
% ans =
%
%   struct with fields:
%
%           Type: 'Gaussian'
%     Parameters: [296 120]
%         Bounds: [-64 656]