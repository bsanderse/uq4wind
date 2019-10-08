function KLDiv = uq_inversion_test_MCMC_computeKL(BayesianModel,Conj)
% UQ_INVERSION_TEST_MCMC_COMPUTEKL estimates the KL-divergence between 
%   the analytical Gaussian posterior and samples generated by the 
%   inversion module. Used in many selftests of the Bayesian inversion
%   module.
%
%   See also: UQ_SELFTEST_UQ_INVERSION, UQ_INVERSION_TEST_MCMC_SETUP

%extract information from input
nDim = 1;
CombinedSamples = reshape(permute(BayesianModel.Results.Sample(floor(end/2):end,:,:),[2 1 3]),nDim,[])';
Approx.Mean = mean(CombinedSamples);
Approx.Cov = cov(CombinedSamples);

%compute divergence
cov1 = Conj.posteriorVariance;
cov2 = Approx.Cov;
mu1 = Conj.posteriorMean;
mu2 = Approx.Mean;
KLDiv = 1/2*(log(det(cov2)/det(cov1)) - nDim + ...
    trace(cov2\cov1) + (mu2-mu1)*(cov2\(mu2'-mu1')));
