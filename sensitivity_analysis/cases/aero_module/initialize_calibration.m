%% Turbine input parameters
% Name of Matlab file representing the turbine data for calibration
turbineName = 'NM80_calibrate'; % 'NM80', 'AVATAR'
% check |NM80_calibrate.m| or |(TurbineName_calibrate).m| for turbine-specific
% settings and definition of uncertainties

%% Forward model description 
% Name of Matlab file representing the model
Model.mHandle = @aero_module_calibration;
% Optionally, one can pass parameters to model stored in the cell array P
P = getParameterAeroModule(turbineName);
Model.Parameters = P;
Model.isVectorized = false;

%% Experimental data
% Marco's (ECN) script for reading the data in N/m
filename_exp = ('../../../Experimental/WINDTRUE/raw.dat');
output_raw = read_exp_data(filename_exp, 2);
% Because the model has different discrepancy options at different radial locations, 
% the measurement data is stored in four different data structures:
Data(1).y = [mean(output_raw.Fy03)]; % [N/m]
Data(1).Name = 'Fy03';
Data(1).MOMap = 1; % Model Output Map 1

Data(2).y = [mean(output_raw.Fy05)]; % [N/m]
Data(2).Name = 'Fy05';
Data(2).MOMap = 2; % Model Output Map 2

Data(3).y = [mean(output_raw.Fy08)]; % [N/m]
Data(3).Name = 'Fy08';
Data(3).MOMap = 3; % Model Output Map 3

Data(4).y = [mean(output_raw.Fy10)]; % [N/m]
Data(4).Name = 'Fy10';
Data(4).MOMap = 4; % Model Output Map 4

%% Surrogate model options
% Switch for Bayesian analysis with the AeroModule or with the surrogate model
Bayes_full = 0; % 0: use surrogate model (PCE); 1: run full model for Bayes (Computationally expensive!)

% If Bayes_full = 0, we need to specify options for loading a surrogate model
Surrogate_model_type = 0; % 0: Uses a stored PCE surrogate model, 1: create surrogate model

% Options for loading a surrogate model
Surrogate_model_filename = 'surrogate/PCE_LARS.mat'; % Specify the surrogate model file to be used

% Options for creating a surrogate model
% These are used if Bayes_full = 0 and Surrogate_model_type = 1
MetaOpts.Type = 'Metamodel';
MetaOpts.MetaType = 'PCE';
MetaOpts.Method = 'LARS'; % Quadrature, OLS, LARS

MetaOpts.ExpDesign.Sampling = 'LHS';
MetaOpts.ExpDesign.NSamples = 100;
MetaOpts.Degree = 1:4;
MetaOpts.TruncOptions.qNorm = 0.75;


%% Likelihood description 
% Here, the discrepancy for each data structure |y| are chosen to be
% independent and identically distributed Gaussian random variables.
% For the current case, 2*standard deviations of the experimental
% measurements is chosen as the prior.
DiscrepancyPriorOpts1.Name = 'Prior of sigma 1';
DiscrepancyPriorOpts1.Marginals(1).Name = 'Sigma1';
DiscrepancyPriorOpts1.Marginals(1).Type = 'Uniform';
DiscrepancyPriorOpts1.Marginals(1).Parameters = [0, 2*std(output_raw.Fy03)];
DiscrepancyPrior1 = uq_createInput(DiscrepancyPriorOpts1);

DiscrepancyOpts(1).Type = 'Gaussian';
DiscrepancyOpts(1).Prior = DiscrepancyPrior1;

DiscrepancyPriorOpts2.Name = 'Prior of sigma 2';
DiscrepancyPriorOpts2.Marginals(1).Name = 'Sigma2';
DiscrepancyPriorOpts2.Marginals(1).Type = 'Uniform';
DiscrepancyPriorOpts2.Marginals(1).Parameters = [0, 2*std(output_raw.Fy05)];
DiscrepancyPrior2 = uq_createInput(DiscrepancyPriorOpts2);

DiscrepancyOpts(2).Type = 'Gaussian';
DiscrepancyOpts(2).Prior = DiscrepancyPrior2;

DiscrepancyPriorOpts3.Name = 'Prior of sigma 3';
DiscrepancyPriorOpts3.Marginals(1).Name = 'Sigma3';
DiscrepancyPriorOpts3.Marginals(1).Type = 'Uniform';
DiscrepancyPriorOpts3.Marginals(1).Parameters = [0, 2*std(output_raw.Fy08)];
DiscrepancyPrior3 = uq_createInput(DiscrepancyPriorOpts3);

DiscrepancyOpts(3).Type = 'Gaussian';
DiscrepancyOpts(3).Prior = DiscrepancyPrior3;

DiscrepancyPriorOpts4.Name = 'Prior of sigma 4';
DiscrepancyPriorOpts4.Marginals(1).Name = 'Sigma4';
DiscrepancyPriorOpts4.Marginals(1).Type = 'Uniform';
DiscrepancyPriorOpts4.Marginals(1).Parameters = [0, 2*std(output_raw.Fy10)];
DiscrepancyPrior4 = uq_createInput(DiscrepancyPriorOpts4);

DiscrepancyOpts(4).Type = 'Gaussian';
DiscrepancyOpts(4).Prior = DiscrepancyPrior4;



%% Bayes options
% MCMC parameters
Solver.Type = 'MCMC';
% MCMC algorithms available in UQLab
MH = 0; % Metropolis-Hastings
AM = 0; % Adaptive Metropolis
AIES = 1; % Affine invariant ensemble
HMC = 0; % Hamilton Monte Carlo

if (MH==1)
    Solver.MCMC.Sampler = 'MH';
    Solver.MCMC.Steps = 1e2;
    Solver.MCMC.NChains = 1e2;
    Solver.MCMC.T0 = 1e1;
end

if (AM==1)
    Solver.MCMC.Sampler = 'AM';
    Solver.MCMC.Steps = 1e2;
    Solver.MCMC.NChains = 1e2;
    Solver.MCMC.T0 = 1e1;
    Solver.MCMC.Epsilon = 1e-2;
end

if (AIES==1)
    Solver.MCMC.Sampler = 'AIES';
    Solver.MCMC.Steps = 1e2;
    Solver.MCMC.NChains = 1e2;
    Solver.MCMC.a = 5;
end

if (HMC==1)
    Solver.MCMC.Sampler = 'HMC';
    Solver.MCMC.LeapfrogSteps = 1e3;
    Solver.MCMC.LeapfrogSize = 0.01;
    Solver.MCMC.Mass = 100;
end


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
%             Prior.Marginals(i).Bounds = P{25}.Marginals(j).Bounds;
            
            if(P{26}{i}{3} ==1) % Get the index and parameter of discrete 2*stdiable
                discrete_index = [discrete_index i];
                discrete_param_vals = [discrete_param_vals Input.Marginals(i).Parameters(2)];
            else
                cont_index = [cont_index i];
            end
            break;
        end
    end
end


