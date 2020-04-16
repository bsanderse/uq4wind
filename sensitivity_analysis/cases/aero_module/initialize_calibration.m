% Name of Matlab file representing the turbine data
turbineName = 'NM80_calibrate'; % 'NM80', 'AVATAR'
% check NM80.m or AVATAR.m or (turbine_name).m for turbine-specific
% settings and definition of uncertainties

%% model description 
% Name of Matlab file representing the model
Model.mHandle = @aero_module;
% Optionally, one can pass parameters to the model stored in the cell
% array P
P = getParameterAeroModule(turbineName);
Model.Parameters = P;
Model.isVectorized = false;

%% data description
Data.y = 96832; % need to put in data (N/m)
Data.Name = 'Axial force';

%% likelihood description
% for simplicity, assume a value for sigma, i.e. the standard deviation
% between model output and data
% this needs to be changed! sigma should be part of the calibration
sigma = 0.1;
DiscrepancyOptsKnown.Type = 'Gaussian';
DiscrepancyOptsKnown.Parameters = sigma^2; % this is sigma^2

%% Bayes options
Solver.Type = 'MCMC';
% Adaptive Metropolis:
Solver.MCMC.Sampler = 'MH';
Solver.MCMC.Steps = 1e2;
Solver.MCMC.NChains = 1e1;
%Solver.MCMC.T0 = 1e2;
%     Solver.MCMC.Proposal.PriorScale = 0.1;
% AIES:
%     Solver.MCMC.Sampler = 'AIES';

% show MCMC chain convergence:
%     Solver.MCMC.Visualize.Parameters = 1;
%     Solver.MCMC.Visualize.Interval = 100;

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


