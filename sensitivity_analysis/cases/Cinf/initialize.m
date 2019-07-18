%% model description 

% infinitely continuous function
% see Patterson,
% ﻿On the construction of a practical Ermakov-Zolotukhin multiple integrator
% integral 5

% name of Matlab file representing the model
Model.mHandle = @Cinf;



%% list of UQ methods to be used for analysis

% specify a list of options from the following list:
% methods = {'MC','PCE_Quad','PCE_OLS','PCE_LARS'};
methods = {'MC','PCE_Quad'};

% for Monte Carlo, specify number of times to repeat MC-based methods to obtain 'nice' convergence
% graphs
MC_repeat = 1;
NsamplesMC = [1e1 1e2 1e3 1e4 1e5];

% for PCE-Quad, specify the polynomial degrees to be tested
DegreesQuad = 1:6; %[1 2 3 4 5 6];


%% input description

ndim = 2;

for ii = 1 : ndim
    Input.Marginals(ii).Type = 'Uniform'; 
    Input.Marginals(ii).Parameters = [0, 1];
    Input.Marginals(ii).Bounds = [0, 1]; % not necessary for uniform but useful for plotting
end

% exact mean:
mean_exact = (2^ndim)*cos(0.6*pi + 0.5*(ndim*0.5))*(sin(0.25)^ndim)/(0.5^ndim);
mean_ref   = mean_exact; % used to scale error