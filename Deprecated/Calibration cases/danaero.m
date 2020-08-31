clc
close all
clearvars
format compact
format long

%% DANAERO experimental data
F_max = 1000; % Maximum azimuthally averaged spanwise force (N)

N_out = 4; % number of radial locations at which we have data
N_points = [20; 30; 10; 5]; % number of measurement data points per radial location

% Normal force (experimental data) in Newtons
% generate some (random) data points
F_n1 = 0.3 + 0.02*randn(N_points(1),1);
F_n2 = 0.55 + 0.02*randn(N_points(2),1);
F_n3 = 0.82 + 0.02*randn(N_points(3),1);
F_n4 = 0.87 + 0.02*randn(N_points(4),1);

% vector with all measurement data
F_n = [F_n1;F_n2;F_n3;F_n4]*F_max;
N_data = length(F_n); % total number of data points

% generate the indices of the data points
ind_dum      = cumsum(N_points);
ind_data_vec = [[1;ind_dum(1:end-1)+1] ind_dum];

% independent variables: at each radial location
% coordinate
x_r    = [11.87;17.82; 28.97; 35.53]; % in meters
x      = repelem(x_r,N_points);
% Beta angle distribution
beta_r = [65; 70; 75; 77]*(pi/180); % in radians
beta   = repelem(beta_r,N_points); % assuming constant across distribution
% Chord distribution
c_r    = [0.9; 0.6; 0.4; 0.3];
c      = repelem(c_r,N_points); % assuming constant across distribution
% Relative velocity distribution
W_r    = [22; 35; 50; 56];
W      = repelem(W_r,N_points); % assuming constant across distribution

% matrix with all the independent variables:
P      = [beta_r c_r W_r x_r];

figure
plot(x, F_n, 'kx')
xlabel('x [m]')
ylabel('F [N]')
title('DANAERO Experiment')
legend('Experiment');

%% Frequentist approach
p = 2*N_out; % Number of parameters to calibrate
% A_single = zeros(N_out,p);
A = zeros(N_data,p); % size of design matrix, normally N_data>p
q_test = eye(p,p);

% loop over possible parameter values
for j=1:p
    % matrix changes depending on the output of the model
    out = dan_model(q_test(j,:),P);
    % loop over different model outputs
    for i=1:N_out
        row = out(i,1);
        % same matrix for each data point:
        A(ind_data_vec(i,1):ind_data_vec(i,2),j) = row;
    end
end
% Ordinary least square
q_freq = pinv(A)*F_n;
%
cl_freq = q_freq(1:2:end);
cd_freq = q_freq(2:2:end);
F_n_freq = dan_model(q_freq,P);
figure
plot(x, F_n, 'kx')
hold on
plot(x_r,F_n_freq, 'g--x')
xlabel('x [m]')
ylabel('F [N]')
title('Experiment vs. Frequentist')
legend('Experiment','Frequentist');
%% Bayesian

% start uqlab
uqlab

Pt1 = 0; % Point1
Pt2 = 0; % Point 2
Pt3 = 0; % Point 3
Pt4 = 1; % Point 4

if (Pt1==1)
    % Prior
    PriorOpts.Marginals(1).Name = 'beta';
    PriorOpts.Marginals(1).Type = 'Constant';
    PriorOpts.Marginals(1).Parameters = beta_r(1);
    %
    PriorOpts.Marginals(2).Name = 'c';
    PriorOpts.Marginals(2).Type = 'Constant';
    PriorOpts.Marginals(2).Parameters = c_r(1);
    %
    PriorOpts.Marginals(3).Name = 'W';
    PriorOpts.Marginals(3).Type = 'Constant';
    PriorOpts.Marginals(3).Parameters = W_r(1);
    %
    PriorOpts.Marginals(4).Name = 'cL';
    PriorOpts.Marginals(4).Type = 'Gaussian';
    PriorOpts.Marginals(4).Parameters = [1 0.1];
    %
    PriorOpts.Marginals(5).Name = 'cD';
    PriorOpts.Marginals(5).Type = 'Gaussian';
    PriorOpts.Marginals(5).Parameters = [0.45 0.1];
    %
    myPriorDist = uq_createInput(PriorOpts);
    % Display input properties
    uq_print(myPriorDist);
    uq_display(myPriorDist);
    
    % Forward model
    ModelOpts1.Name = 'Axial force';
    ModelOpts1.mFile = 'dan_model_bayes';
    forwardModels(1).Model = uq_createModel(ModelOpts1);
    forwardModels(1).PMap = [1 2 3 4 5];
    
    % Experimental data
    myData(1).y = F_n1*F_max;
    myData(1).Name = 'Axial force 1';
    myData(1).MOMap = 1;
    
    % Discrepancy
    DiscrepancyOpts(1).Type = 'Gaussian';
    DiscrepancyOpts(1).Parameters = 1e-2; % known discrepancy variance
    BayesOpts.Discrepancy = DiscrepancyOpts;
end


if (Pt2==1)
    % Prior
    PriorOpts.Marginals(1).Name = 'beta';
    PriorOpts.Marginals(1).Type = 'Constant';
    PriorOpts.Marginals(1).Parameters = beta_r(2);
    %
    PriorOpts.Marginals(2).Name = 'c';
    PriorOpts.Marginals(2).Type = 'Constant';
    PriorOpts.Marginals(2).Parameters = c_r(2);
    %
    PriorOpts.Marginals(3).Name = 'W';
    PriorOpts.Marginals(3).Type = 'Constant';
    PriorOpts.Marginals(3).Parameters = W_r(2);
    %
    PriorOpts.Marginals(4).Name = 'cL';
    PriorOpts.Marginals(4).Type = 'Gaussian';
    PriorOpts.Marginals(4).Parameters = [1.1 0.1];
    %
    PriorOpts.Marginals(5).Name = 'cD';
    PriorOpts.Marginals(5).Type = 'Gaussian';
    PriorOpts.Marginals(5).Parameters = [0.41 0.1];
    %
    myPriorDist = uq_createInput(PriorOpts);
    % Display input properties
    uq_print(myPriorDist);
    uq_display(myPriorDist);
    
    % Forward model
    ModelOpts1.Name = 'Axial force';
    ModelOpts1.mFile = 'dan_model_bayes';
    forwardModels(1).Model = uq_createModel(ModelOpts1);
    forwardModels(1).PMap = [1 2 3 4 5];
    
    % Experimental data
    myData(1).y = F_n2*F_max;
    myData(1).Name = 'Axial force 2';
    myData(1).MOMap = 1;
    
    % Discrepancy
    DiscrepancyOpts(1).Type = 'Gaussian';
    DiscrepancyOpts(1).Parameters = 1e-2; % known discrepancy variance
    BayesOpts.Discrepancy = DiscrepancyOpts;
end

if (Pt3==1)
    % Prior
    PriorOpts.Marginals(1).Name = 'beta';
    PriorOpts.Marginals(1).Type = 'Constant';
    PriorOpts.Marginals(1).Parameters = beta_r(3);
    %
    PriorOpts.Marginals(2).Name = 'c';
    PriorOpts.Marginals(2).Type = 'Constant';
    PriorOpts.Marginals(2).Parameters = c_r(3);
    %
    PriorOpts.Marginals(3).Name = 'W';
    PriorOpts.Marginals(3).Type = 'Constant';
    PriorOpts.Marginals(3).Parameters = W_r(3);
    %
    PriorOpts.Marginals(4).Name = 'cL';
    PriorOpts.Marginals(4).Type = 'Gaussian';
    PriorOpts.Marginals(4).Parameters = [1.3 0.1];
    %
    PriorOpts.Marginals(5).Name = 'cD';
    PriorOpts.Marginals(5).Type = 'Gaussian';
    PriorOpts.Marginals(5).Parameters = [0.35 0.1];
    %
    myPriorDist = uq_createInput(PriorOpts);
    % Display input properties
    uq_print(myPriorDist);
    uq_display(myPriorDist);
    
    % Forward model
    ModelOpts1.Name = 'Axial force';
    ModelOpts1.mFile = 'dan_model_bayes';
    forwardModels(1).Model = uq_createModel(ModelOpts1);
    forwardModels(1).PMap = [1 2 3 4 5];
    
    % Experimental data
    myData(1).y = F_n3*F_max;
    myData(1).Name = 'Axial force 3';
    myData(1).MOMap = 1;
    
    % Discrepancy
    DiscrepancyOpts(1).Type = 'Gaussian';
    DiscrepancyOpts(1).Parameters = 1e-2; % known discrepancy variance
    BayesOpts.Discrepancy = DiscrepancyOpts;
end

if (Pt4==1)
    % Prior
    PriorOpts.Marginals(1).Name = 'beta';
    PriorOpts.Marginals(1).Type = 'Constant';
    PriorOpts.Marginals(1).Parameters = beta_r(4);
    %
    PriorOpts.Marginals(2).Name = 'c';
    PriorOpts.Marginals(2).Type = 'Constant';
    PriorOpts.Marginals(2).Parameters = c_r(4);
    %
    PriorOpts.Marginals(3).Name = 'W';
    PriorOpts.Marginals(3).Type = 'Constant';
    PriorOpts.Marginals(3).Parameters = W_r(4);
    %
    PriorOpts.Marginals(4).Name = 'cL';
    PriorOpts.Marginals(4).Type = 'Gaussian';
    PriorOpts.Marginals(4).Parameters = [1.45 0.1];
    %
    PriorOpts.Marginals(5).Name = 'cD';
    PriorOpts.Marginals(5).Type = 'Gaussian';
    PriorOpts.Marginals(5).Parameters = [0.33 0.1];
    %
    myPriorDist = uq_createInput(PriorOpts);
    % Display input properties
    uq_print(myPriorDist);
    uq_display(myPriorDist);

%     % Surrogate model
%     MetaOpts.Type = 'Metamodel';
%     MetaOpts.MetaType = 'PCE';
%     MetaOpts.ExpDesign.NSamples = 50;
%     mySurrogateModel = uq_createModel(MetaOpts);
    
    % Forward model
    ModelOpts1.Name = 'Axial force';
    ModelOpts1.mFile = 'dan_model_bayes';
    forwardModels(1).Model = uq_createModel(ModelOpts1);
    forwardModels(1).PMap = [1 2 3 4 5];

    % Experimental data
    myData(1).y = F_n4*F_max;
    myData(1).Name = 'Axial force 4';
    myData(1).MOMap = 1;
    
    % Discrepancy
    DiscrepancyOpts(1).Type = 'Gaussian';
    DiscrepancyOpts(1).Parameters = 1e-2; % known discrepancy variance
    BayesOpts.Discrepancy = DiscrepancyOpts;
end

%% Perform Bayesian analysis
BayesOpts.Data = myData;
BayesOpts.Type = 'Inversion';

% MCMC parameters
BayesOpts.Solver.Type = 'MCMC';
% MCMC algorithms available in UQLab
MH = 0; % Metropolis-Hastings
AM = 1; % Adaptive Metropolis
AIES = 0; % Affine invariant ensemble
HMC = 0; % Hamilton Monte Carlo

if (MH==1)
    BayesOpts.Solver.MCMC.Sampler = 'MH';
    BayesOpts.Solver.MCMC.Steps = 1e3;
    BayesOpts.Solver.MCMC.NChains = 1e2;
    BayesOpts.Solver.MCMC.T0 = 1e1;
end

if (AM==1)
    BayesOpts.Solver.MCMC.Sampler = 'AM';
    BayesOpts.Solver.MCMC.Steps = 1e3;
    BayesOpts.Solver.MCMC.NChains = 1e2;
    BayesOpts.Solver.MCMC.T0 = 1e1;
    BayesOpts.Solver.MCMC.Epsilon = 1e-4;
end

if (AIES==1)
    BayesOpts.Solver.MCMC.Sampler = 'AIES';
    BayesOpts.Solver.MCMC.Steps = 1e3;
    BayesOpts.Solver.MCMC.NChains = 1e2;
    BayesOpts.Solver.MCMC.a = 2;
end

if (HMC==1)
    BayesOpts.Solver.MCMC.Sampler = 'HMC';
    BayesOpts.Solver.MCMC.LeapfrogSteps = 1e3;
    BayesOpts.Solver.MCMC.LeapfrogSize = 0.01;
    BayesOpts.Solver.MCMC.Mass = 100;
end


%% Surrogate
MetaOpts.Type = 'Metamodel';
MetaOpts.MetaType = 'PCE';
MetaOpts.ExpDesign.NSamples = 3;
mySurrogateModel = uq_createModel(MetaOpts);

BayesOpts.ForwardModel.Model = mySurrogateModel;
myBayesianAnalysis_surrogateModel = uq_createAnalysis(BayesOpts);
uq_print(myBayesianAnalysis_surrogateModel)
uq_display(myBayesianAnalysis_surrogateModel)


%% Postprocessing

myBayesianAnalysis = uq_createAnalysis(BayesOpts);
uq_print(myBayesianAnalysis)
uq_display(myBayesianAnalysis)
uq_display(myBayesianAnalysis, 'meanConvergence', 'all')
uq_display(myBayesianAnalysis, 'trace', 'all')
uq_display(myBayesianAnalysis, 'acceptance', 'true')
uq_postProcessInversion(myBayesianAnalysis,'pointEstimate', 'MAP')
uq_postProcessInversion(myBayesianAnalysis,'gelmanRubin', 'true')
R_hat = myBayesianAnalysis.Results.PostProc.MPSRF;

if R_hat <= 1.5
    disp('The MCMC simulation has converged')
else
    disp('The MCMC simulation has not converged. Increase the number of samples or fine tune the algorithm.')
end

%% get the Maximum a Posteriori value (MAP)



q_Bayes = myBayesianAnalysis.Results.PostProc.PointEstimate.X;



if (Pt1==1)
    cL_freq = cl_freq(1)
    cD_freq = cd_freq(1)
    cL_bayes = q_Bayes(1)
    cD_bayes = q_Bayes(2)
    F_n_Bayes = 0.5*1.225.*c_r(1)*(W_r(1)^2)*(q_Bayes(1)*sin(beta_r(1))...
        + q_Bayes(2)*cos(beta_r(1)));
    figure
    plot(x, F_n, 'kx')
    hold on
    plot(x_r,F_n_freq, 'g--x')
    hold on
    plot(x_r(1),F_n_Bayes, 'ro', 'MarkerSize',10)
    xlabel('x [m]')
    ylabel('F [N]')
    title('Experiment vs. Frequentist vs. Bayesian')
    legend('Experiment','Frequentist','Bayesian');
end

if (Pt2==1)
    cL_freq = cl_freq(2)
    cD_freq = cd_freq(2)
    cL_bayes = q_Bayes(1)
    cD_bayes = q_Bayes(2)
    F_n_Bayes = 0.5*1.225.*c_r(2)*(W_r(2)^2)*(q_Bayes(1)*sin(beta_r(2))...
        + q_Bayes(2)*cos(beta_r(2)));
    figure
    plot(x, F_n, 'kx')
    hold on
    plot(x_r,F_n_freq, 'g--x')
    hold on
    plot(x_r(2),F_n_Bayes, 'ro', 'MarkerSize',10)
    xlabel('x [m]')
    ylabel('F [N]')
    title('Experiment vs. Frequentist vs. Bayesian')
    legend('Experiment','Frequentist','Bayesian');
end

if (Pt3==1)
    cL_freq = cl_freq(3)
    cD_freq = cd_freq(3)
    cL_bayes = q_Bayes(1)
    cD_bayes = q_Bayes(2)
    F_n_Bayes = 0.5*1.225.*c_r(3)*(W_r(3)^2)*(q_Bayes(1)*sin(beta_r(3))...
        + q_Bayes(2)*cos(beta_r(3)));
    figure
    plot(x, F_n, 'kx')
    hold on
    plot(x_r,F_n_freq, 'g--x')
    hold on
    plot(x_r(3),F_n_Bayes, 'ro', 'MarkerSize',10)
    xlabel('x [m]')
    ylabel('F [N]')
    title('Experiment vs. Frequentist vs. Bayesian')
    legend('Experiment','Frequentist','Bayesian');
end

if (Pt4==1)
    cL_freq = cl_freq(4)
    cD_freq = cd_freq(4)
    cL_bayes = q_Bayes(1)
    cD_bayes = q_Bayes(2)
    F_n_Bayes = 0.5*1.225.*c_r(4)*(W_r(4)^2)*(q_Bayes(1)*sin(beta_r(4))...
        + q_Bayes(2)*cos(beta_r(4)));
    figure
    plot(x, F_n, 'kx')
    hold on
    plot(x_r,F_n_freq, 'g--x')
    hold on
    plot(x_r(4),F_n_Bayes, 'ro', 'MarkerSize',10)
    xlabel('x [m]')
    ylabel('F [N]')
    title('Experiment vs. Frequentist vs. Bayesian')
    legend('Experiment','Frequentist','Bayesian');
end



