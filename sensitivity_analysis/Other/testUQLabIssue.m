%% test UQLab bug
clear all

uqlab;

%% model
Model.mString = 'X(:,1)';
myModel = uq_createModel(Model);

%% prior
PriorOpts.Marginals(1).Type = 'Gaussian';
PriorOpts.Marginals(1).Parameters = [0 1];
PriorOpts.Marginals(1).Bounds = [-3 3];
myPrior = uq_createInput(PriorOpts);


%% data 
myData(1).y = [0; 1; 2];

%% likelihood and hyperparameters
DiscrepancyPriorOpts1.Name = 'Prior of sigma 1';
DiscrepancyPriorOpts1.Marginals(1).Name = 'Sigma1';
DiscrepancyPriorOpts1.Marginals(1).Type = 'Uniform';
DiscrepancyPriorOpts1.Marginals(1).Parameters = [0 2];
DiscrepancyPrior1 = uq_createInput(DiscrepancyPriorOpts1);

DiscrepancyOpts(1).Type = 'Gaussian';
DiscrepancyOpts(1).Prior = DiscrepancyPrior1;


%% Bayes
BayesOpts.Prior = myPrior;                  % Prior
BayesOpts.Data  = myData;                     % Measurement data
BayesOpts.Type  = 'Inversion';              % Calibration
BayesOpts.Discrepancy = DiscrepancyOpts;    % Likelihood


%
% This is a bug fix for UQLab
% UQLab gives an error if the number of field names in the prior are not
% the same as in the discrepancy options
% see https://uqworld.org/t/bayesian-inference-error-when-fields-of-prior-are-different-from-field-of-discrepancy-options/937

% Here we set the fields of the discrepancy to be same as those of the prior
for i=1:length(myPrior.Marginals)
    fieldnames_prior = fieldnames(myPrior.Marginals(i));
    fieldnames_disc  = fieldnames(DiscrepancyOpts(i).Prior.Marginals);
    match_fieldnames = isfield(DiscrepancyOpts(i).Prior.Marginals,fieldnames_prior);
    ind = find(match_fieldnames==0);
    for j=1:length(ind)
        add_field = fieldnames_prior(ind(j));
        DiscrepancyOpts(i).Prior.Marginals.(add_field{1}) = [];
    end
end

BayesianAnalysis = uq_createAnalysis(BayesOpts);
