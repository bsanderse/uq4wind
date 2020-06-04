% Name of Matlab file representing the turbine data
turbineName = 'NM80_calibrate'; % 'NM80', 'AVATAR'
% check NM80.m or AVATAR.m or (turbine_name).m for turbine-specific
% settings and definition of uncertainties

%% model description 
% Name of Matlab file representing the model
Model.mHandle = @aero_module_axial;
% Optionally, one can pass parameters to model stored in the cell array P
P = getParameterAeroModule(turbineName);
Model.Parameters = P;
Model.isVectorized = false;

%% Experimental data
% Marco's script reading data in N/m
filename_exp = ('../../../Experimental/WINDTRUE/raw.dat');
output_raw = importfile3(filename_exp, 2);
Data.y = [output_raw.Fy03, output_raw.Fy05, output_raw.Fy08, output_raw.Fy10]; % Raw data
% Data.y = [mean(output_raw.Fy03), mean(output_raw.Fy05), mean(output_raw.Fy08),...
%           mean(output_raw.Fy10)]; % Mean data
Data.Name = 'Axial force';
%% Surrogate model
load_surrogate = 0;

if (load_surrogate == 1)
    load surrogate.mat;
    
elseif (load_surrogate == 0)

    MetaOpts.Type = 'Metamodel';
    MetaOpts.MetaType = 'PCE';
    MetaOpts.Method = 'LARS'; % Quadrature, OLS, LARS

    MetaOpts.ExpDesign.Sampling = 'LHS';
    MetaOpts.ExpDesign.NSamples = 5;
    MetaOpts.Degree = 1:4;
    MetaOpts.TruncOptions.qNorm = 0.75;
    
end

%% likelihood description
% for simplicity, assume a value for sigma, i.e. the standard deviation
% between model output and data
% this needs to be changed! sigma should be part of the calibration
sigma = 0.1;
DiscrepancyOptsKnown.Type = 'Gaussian';
DiscrepancyOptsKnown.Parameters = sigma^2; % this is sigma^2


% SigmaOpts.Name = 'Prior of sigma2';
% SigmaOpts.Marginals.Name = 'sigma2';
% SigmaOpts.Marginals.Type = 'Uniform';
% SigmaOpts.Marginals.Parameters = [0 20];
% 
% mySigmaDist = uq_createInput(SigmaOpts);
% 
% DiscrepancyOpts.Type = 'Gaussian';
% DiscrepancyOpts.Prior = mySigmaDist;




%% Bayes options
% MCMC parameters
Solver.Type = 'MCMC';
% MCMC algorithms available in UQLab
MH = 1; % Metropolis-Hastings
AM = 0; % Adaptive Metropolis
AIES = 0; % Affine invariant ensemble
HMC = 0; % Hamilton Monte Carlo

if (MH==1)
    Solver.MCMC.Sampler = 'MH';
    Solver.MCMC.Steps = 1e2;
    Solver.MCMC.NChains = 1e2;
    Solver.MCMC.T0 = 1e1;
end

if (AM==1)
    Solver.MCMC.Sampler = 'AM';
    Solver.MCMC.Steps = 1e1;
    Solver.MCMC.NChains = 1e1;
    Solver.MCMC.T0 = 1e1;
    Solver.MCMC.Epsilon = 1e-2;
end

if (AIES==1)
    Solver.MCMC.Sampler = 'AIES';
    Solver.MCMC.Steps = 1e2;
    Solver.MCMC.NChains = 1e2;
    Solver.MCMC.a = 2;
end

if (HMC==1)
    Solver.MCMC.Sampler = 'HMC';
    Solver.MCMC.LeapfrogSteps = 1e3;
    Solver.MCMC.LeapfrogSize = 0.01;
    Solver.MCMC.Mass = 100;
end

BayesOpts.Data = Data;
BayesOpts.Type = 'Inversion';
BayesOpts.Discrepancy = DiscrepancyOptsKnown;
BayesOpts.Solver = Solver;

%% Assemble the Input.Marginal for Bayesian calibration through text comparison
% NOTE: check getParameterAeroModule.m to see the definition of the P array
% P{26} contains the uncertain parameters for which we will do calibration
ndim = length(P{26}); 
% P{25} contains all possible parameters, deterministic and uncertain, of
% which a subset is used in the sensitivity study (as defined in P{25})
ntot = length(P{25}.Marginals); 
discrete_index = [];
cont_index = [];
discrete_param_vals = [];

% loop over P{26} and for each uncertain parameter get the distribution as
% stored in P{25}
for i=1:ndim    
    for j = 1:ntot
        % find which index we need by looking in struct P{25}
        % store the required information in Input.Marginals(i), which will
        % be used by UQLab
        if(strcmp([P{25}.Marginals(j).Name,num2str(P{25}.Marginals(j).Index)],[P{26}{i}{1},num2str(P{26}{i}{2})]))
            Prior.Marginals(i).Name =  [P{26}{i}{1},num2str(P{26}{i}{2})];
            Prior.Marginals(i).Type = P{25}.Marginals(j).Type; 
            Prior.Marginals(i).Parameters = P{25}.Marginals(j).Parameters;
            Prior.Marginals(i).Bounds = P{25}.Marginals(j).Bounds;
            
            if(P{26}{i}{3} ==1) % Get the index and parameter of discrete variable
                discrete_index = [discrete_index i];
                discrete_param_vals = [discrete_param_vals Input.Marginals(i).Parameters(2)];
            else
                cont_index = [cont_index i];
            end
            break;
        end
    end
end


