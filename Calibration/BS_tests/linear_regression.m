%% Sobol indices computed with Monte Carlo and PCE-type methods

clc
close all
clearvars

rng default

% select which methods to run
runOLS = 1; % Ordinary Least Squares
runUQLab = 0; % UQLab
runIQR = 0; % Implicit quadrature rule
runSampling = 1; % Brute-force sampling of prior 

%% generate data and exact solution

% generate artificial measurement data by using the linear model
% with Gaussian noise on top of it

% number of data points
N_data = 10;
% interval [0,1]
x_data = linspace(0,1,N_data)';

% stdev in measurement error - also used for likelihood in Bayes
sigma = 0.1;

% exact value for beta used to generate artificial measurement data, and to
% plot exact solution
beta_exact = [0.5;0.25];
% measurement data
y_data  = linearmodel(beta_exact,x_data) + sigma*randn(N_data,1);

% generate exact solution
N_exact = 100;
x_exact = linspace(0,1,N_exact)';
y_exact = linearmodel(beta_exact,x_exact);


% exact solution for Bayes:
% - multivariate normal likelihood, with mean = 0 and covariance sigma^2 * I
% - uniform prior
% the posterior MAP estimate is then the same as the OLS estimate

figure(1)
plot(x_data,y_data,'x');
hold on
plot(x_exact,y_exact,'-');

%% ordinary least squares (OLS) solution
if (runOLS==1)
    % we now create the design matrix, built with the x_data points
    % we can create it automatically by calling linearmodel several times,
    % each time with only one of the betas active
    p = length(beta_exact); % number of terms to calibrate
    A = zeros(N_data,p); % size of design matrix, normally N_data>p
    beta_test = eye(p,p);
    for i=1:p
        A(:,i) = linearmodel(beta_test(i,:),x_data);
    end
    % solve the (overdetermined) least-squares problem with the regress command
    % this gives the estimate for the beta parameters
    beta_OLS = regress(y_data,A)
    % note that this is the same as
    % beta_OLS = (X' * X)\(X' * y_data)
    % and also the same as
    % beta_OLS = pinv(X)*y_data
    
    % with this estimate of beta we get the following solution
    y_OLS = linearmodel(beta_OLS,x_exact);
    % alternatively, we can use linearmodel_vectorized, which will return the
    % OLS solution at the data points
    % y_OLS = linearmodel_vectorized(beta_OLS,A);
    figure(1)
    plot(x_exact,y_OLS,'--');
end

%% Bayesian solution with UQLab
if (runUQLab == 1)
    
    %% initialize UQlab
    
    % add path
    addpath(genpath('../../../UQLabCore_Rel1.3.0/'));
    % start uqlab
    uqlab
    
    % model settings
    ModelOpts.mFile = 'linearmodel_vectorized';
    % pass design matrix as parameter to the M-file
    ModelOpts.Parameters   = A;
    ModelOpts.isVectorized = true;
    myForwardModel = uq_createModel(ModelOpts);
    
    % prior
    PriorOpts.Marginals(1).Name = 'beta0';
    PriorOpts.Marginals(1).Type = 'Uniform';
    PriorOpts.Marginals(1).Parameters = [0,1];
    %
    PriorOpts.Marginals(2).Name = 'beta1';
    PriorOpts.Marginals(2).Type = 'Uniform';
    PriorOpts.Marginals(2).Parameters = [0,1];
    %
    myPriorDist = uq_createInput(PriorOpts);
    
    
    % display input properties
    uq_print(myPriorDist);
    uq_display(myPriorDist);
    
    
    %% likelihood
    DiscrepancyOptsKnown.Type = 'Gaussian';
    DiscrepancyOptsKnown.Parameters = sigma^2; % this is sigma^2
    
    %% Bayes options
    Solver.Type = 'MCMC';
    Solver.MCMC.Sampler = 'AM';
    Solver.MCMC.Steps = 1e3;
    Solver.MCMC.NChains = 1e2;
    Solver.MCMC.T0 = 1e2;
    % Solver.MCMC.Proposal.PriorScale = 0.1;
    
    myData.y       = y_data'; % note: UQLab uses a row vector here
    myData.Name    = 'measurement data';
    BayesOpts.Data = myData;
    BayesOpts.Type = 'Inversion';
    BayesOpts.Discrepancy = DiscrepancyOptsKnown;
    BayesOpts.Solver = Solver;
    
    %% perform MCMC
    myBayesianAnalysis = uq_createAnalysis(BayesOpts);
    
    %% postprocessing
    % text output of Bayesian analysis
    uq_print(myBayesianAnalysis)
    % graphical display of posterior
    uq_display(myBayesianAnalysis)
    
    % uq_postProcessInversion(myBayesianAnalysis,'pointEstimate', 'MAP')
    
    % get the Maximum a Posteriori value
    beta_MAP = myBayesianAnalysis.Results.PostProc.PointEstimate.X
    
    y_MAP = linearmodel(beta_MAP,x_exact);
    
    figure(1)
    plot(x_exact,y_MAP,'-.');
    legend('measurement data','exact solution','least-squares fit','Bayesian MAP');
end

%% Implicit Quadrature rule
if (runIQR == 1)
    addpath(genpath('/Users/sanderse/Dropbox/work/students/Laurent van den Bos/code/software/'));
    [x, w, err, n] = calibrate_linear(@linearmodel_vectorized,'nz',N_data,'d',p,'sigma',sigma,'A',A,'N',25,'S',1000,'progress',true);
    
    f = quadmarginals(x, w);
end


%% Chebyshev sampling of prior
if (runSampling == 1)
    if (runOLS ==0)
        error('first run OLS');
    end
    addpath(genpath('/Users/sanderse/Dropbox/work/Programming/libs/chebfun-master/'));
    
    % uniform sampling of prior parameter space:
    % assume domain [0,1] for both parameters
    Nbeta1 = 100;
    Nbeta2 = 100;
    beta1_cheb = linspace(0,1,Nbeta1)'; %chebpts(Nbeta1,[0,1]);
    beta2_cheb = linspace(0,1,Nbeta2)'; %chebpts(Nbeta2,[0,1]);
    L = zeros(Nbeta1,Nbeta2); % likelihood
    for i=1:Nbeta1
        for j=1:Nbeta2
            beta_ij = [beta1_cheb(i);beta2_cheb(j)];
            L(i,j) = exp( - (y_data - linearmodel(beta_ij,x_data))'*(y_data - linearmodel(beta_ij,x_data)) / (2*sigma^2) );
        end
    end
%     contourf(beta1_cheb,beta2_cheb,L');
    
    % posterior is directly proportional to likelihood because we use a uniform prior
    % MAP is found by looking for maximum of L
    
    
    % alternatively, we can construct a chebfun on the likelihood directly
    L_cheb = chebfun2( @(beta1,beta2) exp( - (y_data - A*[beta1;beta2])'*(y_data - A*[beta1;beta2]) / (2*sigma^2) ), [0 1 0 1]);
    % and find its maximum:
    [L_max, beta_max] = max2(L_cheb);
    beta_max
    figure
    contourf(L_cheb')
end
