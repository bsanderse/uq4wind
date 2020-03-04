clc
close all
clearvars
format compact
format long

do_Bayes = 0;

%% DANAERO experimental data

F_max = 0.01; % Maximum azimuthally averaged spanwise force (N)

N_out = 4; % number of radial locations at which we have data
N_points = [20; 30; 10; 5]; % number of measurement data points per radial location

% Normal force (experimental data) in Newtons
% generate some (random) data points
F_n1 = 0.3 + 0.02*randn(N_points(1),1);
F_n2 = 0.55 + 0.02*randn(N_points(2),1);
F_n3 = 0.82 + 0.02*randn(N_points(3),1);
F_n4 = 0.87 + 0.02*randn(N_points(4),1);

% F_n1 = [0.27; 0.28; 0.29; 0.3; 0.31 ; 0.32; 0.33];
% F_n2 = [ 0.53; 0.54; 0.55; 0.56; 0.57];
% F_n3 = [0.80; 0.81; 0.82; 0.83; 0.84];
% F_n4 = [ 0.86; 0.87; 0.88];

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
beta_r = [5; 10; 15; 17]*(pi/180); % in radians
beta   = repelem(beta_r,N_points); % assuming constant across distribution
% Chord distribution
c_r    = [0.9; 0.6; 0.4; 0.3];
c      = repelem(c_r,N_points); % assuming constant across distribution
% Relative velocity distribution
v_r    = [22; 35; 50; 56];
v      = repelem(v_r,N_points); % assuming constant across distribution

% matrix with all the independent variables:
P      = [beta_r c_r v_r x_r];

% DANAERO plot (modified)
% figure
% plot(x,F_n,'.')
% xlabel('x [m]')
% ylabel('F_{n} [-]')
% title('Modified experimental data')

%% Frequentist approach
p = 2*N_out; % Number of parameters to calibrate
% A_single = zeros(N_out,p);
A = zeros(N_data,p); % size of design matrix, normally N_data>p
q_test = eye(p,p);

% loop over possible parameter values
for j=1:p
    % matrix changes depending on the output of the model
    out = dan_model_ext(q_test(j,:),P);
    % loop over different model outputs
    for i=1:N_out
        row = out(i,1);
        % same matrix for each data point:
        A(ind_data_vec(i,1):ind_data_vec(i,2),j) = row;
    end
end

q_freq = regress(F_n,A);
F_n_OLS = dan_model_ext(q_freq,P);

cl_freq = (q_freq(1:2:end)./0.5*1.225.*v_r.*c_r)
cd_freq = (q_freq(2:2:end)./0.5*1.225.*v_r.*c_r)


figure(10)
plot(x, F_n, 'kx')
hold on
plot(x_r,F_n_OLS, 'g--x')
xlabel('x [m]')
ylabel('F_n [-]');

%% Bayesian approach
if (do_Bayes == 1)
    
    % start uqlab
    uqlab
    
    % Model settings
    ModelOpts.mFile = 'dan_model_vectorized';
    % pass design matrix as parameter to the M-file
    ModelOpts.Parameters   = A;
    ModelOpts.isVectorized = true;
    myForwardModel = uq_createModel(ModelOpts);
    
    % Prior
    % Note that the range needs to be a good estimate
    PriorOpts.Marginals(1).Name = 'dL*';
    PriorOpts.Marginals(1).Type = 'Uniform';
    PriorOpts.Marginals(1).Parameters = [0,4e-2];
    %
    PriorOpts.Marginals(2).Name = 'dD*';
    PriorOpts.Marginals(2).Type = 'Gaussian';
    PriorOpts.Marginals(2).Parameters = [4e-4,1e-3];
    %
    myPriorDist = uq_createInput(PriorOpts);
    
    
    % Display input properties
    uq_print(myPriorDist);
    uq_display(myPriorDist);
    
    % Discrepancy
    DiscrepancyOpts.Type = 'Gaussian';
    DiscrepancyOpts.Parameters = 1e-6;
    BayesOpts.Discrepancy = DiscrepancyOpts;
    
    % Experimental data
    myData.y       = F_n'; % note: UQLab uses a row vector here
    myData.Name    = 'measurement data';
    BayesOpts.Data = myData;
    BayesOpts.Type = 'Inversion';
    
    
    % MCMC parameters
    BayesOpts.Solver.Type = 'MCMC';
    % MCMC algorithms available in UQLab
    MH = 1; % Metropolis-Hastings
    AM = 0; % Adaptive Metropolis
    AIES = 0; % Affine invariant ensemble
    HMC = 0; % Hamilton Monte Carlo
    
    if (MH==1)
        BayesOpts.Solver.MCMC.Sampler = 'MH';
        BayesOpts.Solver.MCMC.Steps = 1e4;
        BayesOpts.Solver.MCMC.NChains = 1e2;
        BayesOpts.Solver.MCMC.T0 = 1e1;
    end
    
    if (AM==1)
        BayesOpts.Solver.MCMC.Sampler = 'AM';
        BayesOpts.Solver.MCMC.Steps = 1e4;
        BayesOpts.Solver.MCMC.NChains = 1e2;
        BayesOpts.Solver.MCMC.T0 = 1e1;
        BayesOpts.Solver.MCMC.Epsilon = 1e-4;
    end
    
    if (AIES==1)
        BayesOpts.Solver.MCMC.Sampler = 'AIES';
        BayesOpts.Solver.MCMC.Steps = 1e3;
        BayesOpts.Solver.MCMC.NChains = 1e2;
        BayesOpts.Solver.MCMC.T0 = 1e2;
    end
    
    if (HMC==1)
        BayesOpts.Solver.MCMC.Sampler = 'HMC';
        BayesOpts.Solver.MCMC.LeapfrogSteps = 1e3;
        BayesOpts.Solver.MCMC.LeapfrogSize = 0.01;
        BayesOpts.Solver.MCMC.Mass = 1;
    end
    
    % Post-processing
    myBayesianAnalysis = uq_createAnalysis(BayesOpts);
    uq_print(myBayesianAnalysis)
    uq_display(myBayesianAnalysis)
    uq_display(myBayesianAnalysis, 'meanConvergence', 'all')
    uq_display(myBayesianAnalysis, 'trace', 'all')
    uq_display(myBayesianAnalysis, 'acceptance', 'true')
    
    %Point value for Q[q1, q2]
    uq_postProcessInversion(myBayesianAnalysis,'pointEstimate', 'MAP')
    uq_postProcessInversion(myBayesianAnalysis,'gelmanRubin', 'true')
    % get the Maximum a Posteriori value (MAP)
    q_Bayes = myBayesianAnalysis.Results.PostProc.PointEstimate.X;
    F_n_Bayes = dan_model(q_Bayes,beta);
    
    
    
    cl_bayes = (myBayesianAnalysis.Results.PostProc.PointEstimate.X(1)./0.5*1.225.*v_r.*c_r) % MAP
    cd_bayes = (myBayesianAnalysis.Results.PostProc.PointEstimate.X(2)./0.5*1.225.*v_r.*c_r) % MAP
    
    figure(10)
    hold on
    plot(beta,F_n_Bayes, 'r--o')
    xlabel('\beta [radians]')
    ylabel('F_{n} [-]')
    title('Experiment vs. Frequentist vs. Bayesian')
    legend('measurement data','OLS','Bayesian');
    
end
