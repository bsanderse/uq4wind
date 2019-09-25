function [x_cond,output_samples] = uq_getKucherenkoSamples(myInput,N,Method,CondingIdx,samples,Estimator)
% UQ_GETKUCHERENKOSAMPLES produces cross-conditioned samples made from
% two samples
%   [x_cond,output_samples] = uq_getKucherenkoSamples(myInput,N,Method,CondingIdx,samples,Estimator)
%   - CondingIdx: contains the indices of the conditioning variables
%   (numeric or logical, function works with logical)
%   - samples: the input contains provided samples as well as sampling
%   strategies and sample size. The output additionally contains two
%   created uniform, independent samples (u1, u2), uniform dependent
%   samples (u1corr, u2corr) and their transformations into the physical
%   space (x1, x2), that are needed for the Kucherenko indices.
%
% See also: UQ_CLOSED_SENS_INDEX, UQ_TOTAL_SENS_INDEX

%% Setup
% Amount of variables
M = length(myInput.Marginals);

%%
% Check the format of the indices. We want them to be logical.
if islogical(CondingIdx)
    % there should be M entries
    if length(CondingIdx) ~= M
        fprintf('\n\nError: The conditioning indices are provided as logicals, \n but the length of the array is not equal to M!\n');
        error('While initiating the conditional sampling')
    end
elseif isnumeric(CondingIdx)
    % maybe it's numeric 1's and 0's and meant to be logical
    if all(ismember(CondingIdx,[0 1])) && length(CondingIdx)==M
        CondingIdx = logical(CondingIdx);
    
    % but maybe it's variable indices $\subset (1,2,...,M)$. Turn into logical
    elseif all(CondingIdx < M+1) && length(unique(CondingIdx)) == length(CondingIdx)
        logidx = false(1,M);
        logidx(CondingIdx) = true;
        CondingIdx = logidx;
        
    else
        fprintf('\n\nError: The provided conditioning indices are neither logical nor numeric!\n');
        error('While initiating the conditional sampling')
    end
else
    fprintf('\n\nError: The provided conditioning indices are neither logical nor numeric!\n');
    error('While initiating the conditional sampling')
end

%% Get two samples

% setup which samples are still needed
need1all = false;
need1u = false;
need2all = false;
need2u = false;

% Some info on U to use the isoprobabilistic transform. Is always useful.
[U_marginals(1:M).Type] = deal('uniform');
[U_marginals(1:M).Parameters] = deal([0 1]);
U_copula.Type = 'Independent';

% check for available samples. ONLY EXPECTED CASES ARE TREATED / INVESTIGATED.
if ~isfield(samples,'u1') && ~isfield(samples,'x1') % nothing provided
    need1all = true;
    need2all = true;
    
elseif isfield(samples,'x1') % a sample in physical space is provided
    if ~isfield(samples,'u1')
        need1u = true;
    end
    if isfield(samples,'x2') % second sample in physical space is provided
        if ~isfield(samples,'u2')
            need2u = true;
        end
    else
        need2all = true;
    end
end

% Produce the needed samples:
% 'u' uniform indep, 'ucorr' uniform dep, 'x' physical
% sampling options
samplingopt.Method = Method;
if need1all
    samples.u1 = uq_sampleU(N,M,samplingopt);
    samples.u1corr = uq_GeneralIsopTransform(samples.u1, U_marginals, U_copula, U_marginals, myInput.Copula);
    samples.x1 = uq_IsopTransform(samples.u1corr, U_marginals, myInput.Marginals);
end
if need2all
    samples.u2 = uq_sampleU(N,M,samplingopt);
    samples.u2corr = uq_GeneralIsopTransform(samples.u2, U_marginals, U_copula, U_marginals, myInput.Copula);
    samples.x2 = uq_IsopTransform(samples.u2corr, U_marginals, myInput.Marginals);
end
if need1u
    samples.u1corr = uq_IsopTransform(samples.x1, myInput.Marginals, U_marginals);
    samples.u1 = uq_GeneralIsopTransform(samples.u1corr, U_marginals, myInput.Copula, U_marginals, U_copula);
end
if need2u
    samples.u2corr = uq_IsopTransform(samples.x2, myInput.Marginals, U_marginals);
    samples.u2 = uq_GeneralIsopTransform(samples.u2corr, U_marginals, myInput.Copula, U_marginals, U_copula);
end


%% Conditioning
% mix the samples
u_mix = zeros(N,M);
switch Estimator
    case {'standard', 'modified'}
        u_mix(:,CondingIdx) = samples.u1corr(:,CondingIdx);
        u_mix(:,~CondingIdx) = samples.u2(:,~CondingIdx);
        
    case 'alternative'
        u_mix(:,CondingIdx) = samples.u2corr(:,CondingIdx);
        u_mix(:,~CondingIdx) = samples.u1(:,~CondingIdx);
        
end
% get the needed cross conditioning
x_cond = uq_getCondSample(myInput,N,Method,CondingIdx,u_mix,true);

% also return the extended samples structure
output_samples = samples;

end
