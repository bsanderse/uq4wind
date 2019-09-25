function varargout = uq_postProcessInversion(module, varargin)
% UQ_POSTPROCESSINVERSION post-processes an inverse analysis carried out
%    with the Bayesian inversion module of UQLab.
%
%    UQ_DISPLAY_UQ_INVERSION(MODULE, NAME, VALUE) allows to choose
%    more post processing options by specifying Name/Value pairs:
%
%       Name                  VALUE
%       'burnIn'              Removes the burn in from the sample,
%                             specified as a fraction (between 0 and 1) or
%                             an Integer smaller than the number of
%                             MCMC iterations
%                             - Double or Integer
%                             default : 0.5
%       'badChains'           Removes chains specified by their id from the
%                             sample
%                             - Integer
%                             default : 0
%       'pointEstimate'       Computes a point estimate based on the
%                             supplied sample.
%                             - String : ('Mean','MAP')
%                             default : 'mean'
%       'percentiles'         Computes dimension-wise percentiles of 
%                             supplied sample.
%                             - Double
%                             default : [0.05, 0.95]
%       'gelmanRubin'         Computes the Rubin-Gelman scale reduction
%                             factor
%                             - Logical
%                             default : false
%       'prior'               Samples specified number of sample points 
%                             from the prior distribution.
%                             - Integer
%                             default : 10,000
%       'priorPredictive'     Samples specified number of sample points 
%                             from the prior predictive distribution. This 
%                             requires additional calls to the forward 
%                             models
%                             - Integer
%                             default : 1,000
%       'posteriorPredictive' Samples specified number of sample points 
%                             from the posterior predictive distribution.
%                             - Integer
%                             default : 1,000
%                          
% See also: UQ_PRINT_UQ_INVERSION, UQ_DISPLAY_UQ_INVERSION

%% CONSISTENCY CHECKS
if ~strcmp(module.Type, 'uq_inversion')
    error('uq_inversion_postProc only operates on objects of type ''Inversion''') 
end

%check if MCMC Solver
if ~strcmp(module.Internal.Solver.Type, 'MCMC')
    error('No results to post-process')
end

% switch if custom likelihood
if module.Internal.customLikeli
    CUSTOM_LIKELI = true;
else
    CUSTOM_LIKELI = false;
end

%% INITIALIZE
% get Sample
Sample = module.Results.Sample;

% get loglikelihood
LogLikeliEval = module.Results.LogLikeliEval;

% get forward model evaluations
ForwardModelRuns = module.Results.ForwardModel;
nForwardModels = length(ForwardModelRuns);

% get length variables
[nIter,nDim,nChains] = size(Sample);

% get indices of discrepancy Sample excluding constants
nonConstId = module.Internal.ConstInfo.idNonConst;
discrepancyIndices = find(module.Internal.paramDiscrepancyID(nonConstId));
modelIndices = find(module.Internal.paramDiscrepancyID(nonConstId) == 0);

%% Default behavior
% Burn in
Default.burnIn_flag = true;
Default.burnIn = floor(nIter/2); %discard half the chain
% Point estimate
Default.pointEstimate_flag = true;
Default.pointEstimate = 'Mean';
% Prior sample
Default.prior_flag = true;
Default.nPriorSamples = 10000;
% Percentiles
Default.percentiles_flag = true;
Default.perc_probabilities = [0.05, 0.95];
if ~CUSTOM_LIKELI
    % only do predictive distributions if no user-specified likelihood
    % Posterior predictive
    Default.posteriorPredictive_flag = true;
    Default.nPostPredSamples = 1000;
    % Prior predictive
    Default.priorPredictive_flag = true;
    Default.nPriorPredSamples = 1000;
else
    % Posterior predictive
    Default.posteriorPredictive_flag = false;
    Default.nPostPredSamples = 0;
    % Prior predictive
    Default.priorPredictive_flag = false;
    Default.nPriorPredSamples = 0;
end
% Bad chains
Default.badChains_flag = false;
% Gelman Rubin
Default.gelmanRubin_flag = false;

%% Check for input arguments
%set optional arguments

if nargin > 1
    % vargin given
    parse_keys = {'burnIn', 'badChains', 'pointEstimate','percentiles',...
        'gelmanRubin', 'prior','priorPredictive', 'posteriorPredictive'};
    parse_types = {'p','p','p','p','p','p','p','p'};
    [uq_cline, ~] = uq_simple_parser(varargin, parse_keys, parse_types);
else
    % no vargin, use default options
    nOpts = 8;
    uq_cline = cell(nOpts,1);
    for ii = 1:nOpts
        uq_cline{ii} = 'false';
    end
end

% 'burnIn' option removes burnin from the sample before
% plotting
if ~strcmp(uq_cline{1}, 'false')
    burnIn_flag = true;
    if floor(uq_cline{1})==uq_cline{1}
        %integer passed
        burnIn = uq_cline{1};
    elseif and(uq_cline{1} >= 0, uq_cline{1} <= 1)
        %fraction passed
        burnIn = floor(uq_cline{1}*size(module.Results.Sample,1));
    else
        error('Argument of burnIn not valid.')    
    end
else
    burnIn_flag = Default.burnIn_flag;
    burnIn = Default.burnIn;
end

% 'badChains' option removes chains specified by badChains from the
% sample before plotting
if ~strcmp(uq_cline{2}, 'false')
    badChains = uq_cline{2};
    badChains_flag = true;
    if max(badChains > nChains)
        error('Argument of badChains not valid.')    
    end
else
    badChains_flag = Default.badChains_flag;
end

% 'pointEstimate' option adds a point estimate to the plot
if and(ischar(uq_cline{3}),~strcmp(uq_cline{3}, 'false'))
    pointEstimate_flag = true;
    pointEstimate = uq_cline{3};
elseif isnumeric(uq_cline{3})
    pointEstimate_flag = true;
    pointEstimate = 'custom';
    pointParamIn = uq_cline{3};
else
    pointEstimate_flag = Default.pointEstimate_flag;
    pointEstimate = Default.pointEstimate;
end

% 'percentiles' computes the dimensionwise percentiles of the sample
if ~strcmp(uq_cline{4}, 'false')
    percentiles_flag = true;
    perc_probabilities = uq_cline{4};
    if or(max(perc_probabilities) > 1, min(perc_probabilities) < 0)
        % supplied percentiles are out of bounds
        error('The percentiles are not between [0,1].')
    end
else
    percentiles_flag = Default.percentiles_flag;
    perc_probabilities = Default.perc_probabilities;
end

% 'gelmanRubin' computes the Gelman-Rubin scale reduction factor 
if ~strcmp(uq_cline{5}, 'false')
    gelmanRubin_flag = uq_cline{5};
else
    gelmanRubin_flag = Default.gelmanRubin_flag;
end

% 'prior' draws samples from the prior distribution
if ~strcmp(uq_cline{6}, 'false')
    prior_flag = true;
    nPriorSamples = uq_cline{6};
else
    prior_flag = Default.prior_flag;
    nPriorSamples = Default.nPriorSamples;
end

% 'priorPredictive' draws samples from the prior predictive
% distribution
if ~strcmp(uq_cline{7}, 'false')
    if CUSTOM_LIKELI
        error('Predictive distributions are not supported with user-specified likelihood functions.')
    end
    priorPredictive_flag = true;
    nPriorPredSamples = uq_cline{7};
else
    priorPredictive_flag = Default.priorPredictive_flag;
    nPriorPredSamples = Default.nPriorPredSamples;
end

% 'postPredictive' draws samples from the posterior predictive
% distribution
if ~strcmp(uq_cline{8}, 'false')
    if CUSTOM_LIKELI
        error('Predictive distributions are not supported with user-specified likelihood functions.')
    end
    posteriorPredictive_flag = true;
    nPostPredSamples = uq_cline{8};
else
    posteriorPredictive_flag = Default.posteriorPredictive_flag;
    nPostPredSamples = Default.nPostPredSamples;
end

%% PREPROCESSING - MODIFICATIONS OF SAMPLE
if burnIn_flag
    % preprocessing was chosen
    % remove the burn in
    % ensure that burnIn is smaller than the number of Sample
    if burnIn >= size(Sample,1)
        error('Burn in larger than number of sample points per chain')
    end
    Sample = Sample(burnIn:end,:,:);
    
    % logliklihood evaluations
    LogLikeliEval = LogLikeliEval(burnIn:end,:);
    
    % model evaluations
    for mm = 1:nForwardModels
        ForwardModelRuns(mm).evaluation = ...
            ForwardModelRuns(mm).evaluation(burnIn:end,:,:);
    end
    
    %Return to analysis object
    module.Results.PostProc.PostSample = Sample;
    module.Results.PostProc.PostLogLikeliEval = LogLikeliEval;
    module.Results.PostProc.PostModel = ForwardModelRuns;
    
end

if badChains_flag
    % preprocessing was chosen
    % remove the bad chains
    goodChainsIndex = 1:nChains; goodChainsIndex(badChains) = [];
    Sample = Sample(:,:,goodChainsIndex);
    
    % logliklihood evaluations
    LogLikeliEval = LogLikeliEval(:,goodChainsIndex);
    
    % model evaluations
    for mm = 1:nForwardModels
        ForwardModelRuns(mm).evaluation = ...
            ForwardModelRuns(mm).evaluation(:,:,goodChainsIndex);
    end
    
    %Return to analysis object
    module.Results.PostProc.PostSample = Sample;
    module.Results.PostProc.PostLogLikeliEval = LogLikeliEval;
    module.Results.PostProc.PostModel = ForwardModelRuns;
end

% Combine chains
% Sample
Sample3D = Sample;
Sample2D = reshape(permute(Sample,[2 1 3]),nDim,[])';

% LogLikelihood
LogLikelihood3D = LogLikeliEval;
LogLikelihood2D = reshape(LogLikeliEval,[],1);

% Forward model
ForwardModelRuns3D = ForwardModelRuns;
for mm = 1:nForwardModels
    nOut = size(ForwardModelRuns(mm).evaluation,2);
    ForwardModelRuns2D(mm).evaluation = ...
        reshape(permute(ForwardModelRuns(mm).evaluation,[2 1 3]),nOut,[])';
end

if pointEstimate_flag
    switch lower(pointEstimate)
        case 'mean'
            % set to posterior mean
            pointParam = uq_inversion_mean(Sample3D);
        case 'map'
            % take the sample with the maximum (unnormalized) posterior
            % density as the MAP
            LogLikelihood = LogLikelihood2D;
            LogPrior = uq_eval_logpdf(Sample2D,module.Internal.FullPrior);
            [~,maxIndex] = max(LogLikelihood + LogPrior);
            % maximum sample is the MAP
            pointParam = Sample2D(maxIndex,:);            
        case 'custom'
            % do nothing
            pointParam = pointParamIn;
    end
    % Return to analysis object
    module.Results.PostProc.PointEstimate.Parameter = pointParam;
    module.Results.PostProc.PointEstimate.Type = pointEstimate;
end

if percentiles_flag
   % compute percentiles from supplied sample
   percentiles = uq_inversion_percentiles(Sample3D, perc_probabilities);
   
   % return percentiles to analysis object
   module.Results.PostProc.Percentiles.Values = percentiles;
   module.Results.PostProc.Percentiles.Probabilities = perc_probabilities;
end

if gelmanRubin_flag
    % compute potential scale reduction factor (if at least two 
    % chains are available)
    [nSamples,nDim,nChains] = size(Sample3D);
    if nChains > 1
        %within-sequence variance 
        W = zeros(nDim);
        chainMeans = zeros(nChains,nDim);
        for ii = 1:nChains
            samplesCurr = squeeze(Sample3D(:,:,ii));
            W = W+cov(samplesCurr);
            chainMeans(ii,:) = mean(samplesCurr);
        end
        W = W/nChains;

        %between-sequence variance
        B = cov(chainMeans);

        %scale reduction factor
        tempMat = W\B;
        MPSRF = nSamples/(nSamples+1)+(nChains+1)/nChains*eigs(tempMat,1);
    else
        error('Need at least two chains to compute Rubin-Gelman scale reduction factor!')
    end
    %Return to analysis object
    module.Results.PostProc.MPSRF = MPSRF;
end

%% SAMPLING FROM DISTRIBUTIONS

if prior_flag
    % draw prior samples
    nSamples = nPriorSamples;
    PriorSample = uq_getSample(module.Internal.FullPrior ,nSamples);
    
    % Return to analysis object
    module.Results.PostProc.PriorSample = PriorSample;
end

if priorPredictive_flag
    % compute the prior predictive
    % n plot Sample
    nPlotSamples = nPriorPredSamples;
    nSamples = size(Sample2D,1);
    % get reduced number of Sample for plotting
    plotIndices = ceil(rand(nPlotSamples,1)*nSamples);
    % prior Sample
    priorSample = uq_getSample(module.PriorDist,nPlotSamples);
    priorSampleModel = priorSample(:,modelIndices);
    
    %loop over forward models
    for mm = 1:nForwardModels
        %point estimate
        if pointEstimate_flag
            %remove discrepancy terms
            pointParamModel = pointParam;
            pointParamModel(:,discrepancyIndices) = [];
            model(mm).pointEstimateRun = uq_evalModel(module.ForwardModel(mm).Model,pointParamModel);
        end
        %prior Runs
        model(mm).priorRuns = uq_evalModel(module.ForwardModel(mm).Model,priorSampleModel);
        model(mm).priorPredRuns = zeros(size(model(mm).priorRuns));
        %data
        for jj = 1:size(model(mm).priorRuns,2)
            model(mm).data(jj).y = [];
        end
    end
    
    %loop over data groups
    for ii = 1:module.Internal.nDataGroups
        % get currents
        yCurr = module.Internal.Data(ii).y;
        discrepancyCurr = module.Internal.Discrepancy(ii);
        %get model output index
        MOMap = module.Internal.Data(ii).MOMap;
        %use MOMap to extract prior runs relevant for the current
        %data group
        priorRunsCurr = zeros(nPlotSamples,size(yCurr,2));
        priorPredRunsCurr = zeros(nPlotSamples,size(yCurr,2));
        for jj = 1:module.Internal.nForwardModels
            if any(MOMap(1,:)==jj)
                priorRunsCurr = model(jj).priorRuns(:,MOMap(2,MOMap(1,:)==jj));
            end
        end
        %add discrepancy from likelihood Sample to the current prior/priorerior runs 
        if strcmp(discrepancyCurr.Type,'Gaussian')
            if discrepancyCurr.ParamKnown
                % Known discrepancy
                % draw Sample from the conditional distribution
                % loop over plot Sample
                for kk = 1:nPlotSamples
                    priorPredRunsCurr(kk,:) = ...
                        discrepancyCurr.likelihoodSamples_g(priorRunsCurr(kk,:),discrepancyCurr.Parameters);
                end
            else
                % Unknown discrepancy
                %get index
                currParamIndex = module.Internal.paramDiscrepancyID == ii;
                %remove constants
                currParamIndex = currParamIndex(nonConstId);
                priorSamplesDiscrepancyCurr = priorSample(:,currParamIndex);
                %loop over sample points
                for kk = 1:nPlotSamples
                    priorPredRunsCurr(kk,:) = ...
                        discrepancyCurr.likelihoodSamples_g(priorRunsCurr(kk,:),priorSamplesDiscrepancyCurr(kk,:));
                end
            end
        else
            error('This discrepancy model is not supported for this plot')
        end
        
        %assign the current prior runs and the data to the model(i) structure
        for mm = 1:nForwardModels
            if any(MOMap(1,:)==mm)
                MIndex = MOMap(1,:)==mm;
                OMapCurr = MOMap(2,MIndex);
                %predictive runs
                model(mm).priorPredRuns(:,OMapCurr) = priorPredRunsCurr;
                %data related to current model
                for kk = 1:length(OMapCurr)
                    model(mm).data(OMapCurr(kk)).y = ...
                        [model(mm).data(OMapCurr(kk)).y;yCurr(:,kk)];
                end
            end
        end
    end
    
    % Return to analysis object
    module.Results.PostProc.PriorPred.model = model;
end

%clear model variable
clear model

if posteriorPredictive_flag
    % compute the post predictive
    % n plot Sample
    nPlotSamples = nPostPredSamples;
    nSamples = size(Sample2D,1);
    % get reduced number of Sample for plotting
    plotIndices = ceil(rand(nPlotSamples,1)*nSamples);
    % posterior samples
    postSample = Sample2D(plotIndices,:);
    
    %loop over forward models
    for mm = 1:nForwardModels
        %point estimate
        if pointEstimate_flag
            %remove discrepancy terms
            pointParamModel = pointParam;
            pointParamModel(:,discrepancyIndices) = [];
            model(mm).pointEstimateRun = uq_evalModel(module.ForwardModel(mm).Model,pointParamModel);
        end
        %post Runs
        model(mm).postRuns = ForwardModelRuns2D(mm).evaluation(plotIndices,:);
        model(mm).postPredRuns = zeros(size(model(mm).postRuns));
        %data
        for jj = 1:size(model(mm).postRuns,2)
            model(mm).data(jj).y = [];
        end
    end
    
    %loop over data groups
    for ii = 1:module.Internal.nDataGroups
        % get currents
        yCurr = module.Internal.Data(ii).y;
        discrepancyCurr = module.Internal.Discrepancy(ii);
        %get model output index
        MOMap = module.Internal.Data(ii).MOMap;
        %use MOMap to extract posterior runs relevant for the current
        %data group
        postRunsCurr = zeros(nPlotSamples,size(yCurr,2));
        postPredRunsCurr = zeros(nPlotSamples,size(yCurr,2));
        for jj = 1:module.Internal.nForwardModels
            if any(MOMap(1,:)==jj)
                postRunsCurr = model(jj).postRuns(:,MOMap(2,MOMap(1,:)==jj));
            end
        end
        %add discrepancy from likelihood Sample to the current posterior runs 
        if strcmp(discrepancyCurr.Type,'Gaussian')
            if discrepancyCurr.ParamKnown
                % Known discrepancy
                % draw Sample from the conditional distribution
                % loop over plot Sample
                for kk = 1:nPlotSamples
                    postPredRunsCurr(kk,:) = ...
                        discrepancyCurr.likelihoodSamples_g(postRunsCurr(kk,:),discrepancyCurr.Parameters);
                end
            else
                % Unknown discrepancy
                %get index
                currParamIndex = module.Internal.paramDiscrepancyID == ii;
                %remove constants
                currParamIndex = currParamIndex(nonConstId);
                postSamplesDiscrepancyCurr = postSample(:,currParamIndex);
                %loop over sample points
                for kk = 1:nPlotSamples
                    postPredRunsCurr(kk,:) = ...
                        discrepancyCurr.likelihoodSamples_g(postRunsCurr(kk,:),postSamplesDiscrepancyCurr(kk,:));
                end
            end
        else
            error('This discrepancy model is not supported for this computation')
        end
        
        %assign the current post/posterior runs and the data to the model(i) structure
        for mm = 1:nForwardModels
            if any(MOMap(1,:)==mm)
                MIndex = MOMap(1,:)==mm;
                OMapCurr = MOMap(2,MIndex);
                %predictive runs
                model(mm).postPredRuns(:,OMapCurr) = postPredRunsCurr;
                %data related to current model
                for kk = 1:length(OMapCurr)
                    model(mm).data(OMapCurr(kk)).y = ...
                        [model(mm).data(OMapCurr(kk)).y;yCurr(:,kk)];
                end
            end
        end
    end
    
    % Return to analysis object
    module.Results.PostProc.PostPred.model = model;
end

if nargout 
    varargout{1} = module;
end